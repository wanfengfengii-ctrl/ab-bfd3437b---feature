"""环状 DNA 条码拼接核心算法。

把每条长度 L = k+1 的读数看成德布鲁因多重图中的一条有向边
    s[:-1]  ->  s[1:]
（序列每次出现都是一条独立平行边，即一份证据）。
"每条序列恰好使用一次、相邻重叠 k = L-1 的闭环拼接" 等价于该多重图的
欧拉回路。

单株模式（默认 / barcode_count=1）：全部读数必须属于同一条欧拉回路。

双株模式（barcode_count=2）：全部读数来自恰好两条环状条码。算法
**联合枚举完整划分**——把每条读数（实例）分到组 A 或组 B（读数 0 固定归 A
以破除两群无先后的对称），仅当两个组的边子图**同时**度数平衡且弱连通、
各自至少 3 条边时才接受；不是先拼好一条再拿余料兜底。每组内部仍按原欧拉
回路规则枚举闭环。

等价归一化（用于判断是否为同一个可信结论）：
  * 环的循环移位；
  * 整条环的反向互补；
  * 两条条码不分先后（组交换）；
  * 完全相同的重复读数（平行边）互换不算新答案。

因此一个"完整划分的规范等价类"就是两组规范条码组成的无序对。枚举找到 2
个不同等价类即判定为歧义。
"""
from __future__ import annotations

from dataclasses import dataclass

ALPHABET = frozenset("ACGT")
_COMP = str.maketrans("ACGT", "TGCA")

MIN_READS, MAX_READS = 6, 24
DUAL_MAX_READS = 18
MIN_LEN, MAX_LEN = 3, 8
MIN_GROUP_READS = 3


class ValidationError(ValueError):
    """输入不满足 6~24 条、等长 3~8、仅含 ACGT 等约束。"""


@dataclass(frozen=True)
class Edge:
    u: str       # 起点 k-mer
    v: str       # 终点 k-mer
    label: str   # 原始读数
    idx: int     # 从 0 开始的输入序号


def reverse_complement(s: str) -> str:
    return s.translate(_COMP)[::-1]


def canonical_circle(seq: str) -> str:
    """环序列在"循环移位 + 整条反向互补"下的规范代表（字典序最小者）。"""
    n = len(seq)
    best = seq
    s = seq
    for _ in range(n - 1):
        s = s[1:] + s[0]
        if s < best:
            best = s
    s = reverse_complement(seq)
    for _ in range(n):
        if s < best:
            best = s
        s = s[1:] + s[0]
    return best


def validate(sequences: list[str], barcode_count: int = 1) -> int:
    """校验输入，返回 k = L-1；不合法抛 ValidationError。"""
    if not isinstance(sequences, list) or not sequences:
        raise ValidationError("sequences 必须是非空数组")
    n = len(sequences)
    if barcode_count == 1:
        if not (MIN_READS <= n <= MAX_READS):
            raise ValidationError(
                f"序列条数必须在 {MIN_READS}~{MAX_READS} 之间，实际 {n} 条"
            )
    elif barcode_count == 2:
        if not (MIN_READS <= n <= DUAL_MAX_READS):
            raise ValidationError(
                f"双株模式序列条数必须在 {MIN_READS}~{DUAL_MAX_READS} 之间，"
                f"实际 {n} 条"
            )
    else:
        raise ValidationError(f"barcode_count 只支持 1 或 2，实际 {barcode_count!r}")
    norm: list[str] = []
    for i, s in enumerate(sequences):
        if not isinstance(s, str):
            raise ValidationError(f"第 {i + 1} 条序列不是字符串")
        t = s.strip().upper()
        if not (MIN_LEN <= len(t) <= MAX_LEN):
            raise ValidationError(
                f"第 {i + 1} 条序列长度必须在 {MIN_LEN}~{MAX_LEN}，实际 {len(t)}"
            )
        bad = sorted(set(t) - ALPHABET)
        if bad:
            raise ValidationError(f"第 {i + 1} 条序列含非法碱基: {''.join(bad)}")
        norm.append(t)
    lengths = {len(s) for s in norm}
    if len(lengths) != 1:
        raise ValidationError(f"所有序列必须等长，实际长度集合 {sorted(lengths)}")
    # 规范化后写回原列表（大写、去空白）
    sequences[:] = norm
    return len(norm[0]) - 1


def _components(edges: list[Edge]) -> list[list[str]]:
    """非零度顶点的弱连通分量（并查集）。"""
    parent: dict[str, str] = {}

    def find(x: str) -> str:
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: str, b: str) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for e in edges:
        union(e.u, e.v)
    groups: dict[str, list[str]] = {}
    for v in sorted(parent):
        groups.setdefault(find(v), []).append(v)
    return [sorted(g) for g in groups.values()]


def _remaining_connected(mask: int, edges: list[Edge], cur: str, start: str) -> bool:
    """剩余边的所有端点（连同 cur、start）是否处于同一个弱连通分量。

    不满足说明当前走法跨过了桥，余下的边不可能再被一条连续轨迹走完并回到
    起点，可安全剪枝。
    """
    parent: dict[str, str] = {}

    def find(x: str) -> str:
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    m = mask
    while m:
        b = m & -m
        i = b.bit_length() - 1
        e = edges[i]
        ra, rb = find(e.u), find(e.v)
        if ra != rb:
            parent[rb] = ra
        m ^= b

    root = find(cur)
    if find(start) != root:
        return False
    for v in parent:
        if find(v) != root:
            return False
    return True


def _barcode_from_path(path: list[int], edges: list[Edge], k: int) -> str:
    """由边序列拼出长度为 n 的环状条码（闭环后末尾 k-mer 与起点重合）。"""
    s = edges[path[0]].label
    for eid in path[1:]:
        s += edges[eid].label[-1]
    return s[: len(path)]  # 丢掉与起点重合的末尾 k 个字符


def _align_order(
    path: list[int], circle: str, canon: str
) -> tuple[list[int], str]:
    """尽量用纯循环移位把使用次序对齐到规范条码；否则只能经反向互补归一。"""
    n = len(circle)
    s = circle
    for r in range(n):
        if s == canon:
            return path[r:] + path[:r], "cyclic_rotation" if r else "identical"
        s = s[1:] + s[0]
    return path, "reverse_complement"


def _enumerate_classes(
    edges: list[Edge], start: str, k: int, limit: int = 2
) -> list[tuple[str, list[int], str]]:
    """枚举至多 limit 个不同规范等价类，返回 (规范条码, 边位置路径, 归一关系)。

    路径中的数字是 ``edges`` 列表中的**位置**（不一定等于 ``Edge.idx``：
    双株模式会为每个分组构造只含本组边、但保留原输入序号的子列表）。
    """
    m = len(edges)
    full = (1 << m) - 1

    outgoing: dict[str, list[int]] = {}
    for pos, e in enumerate(edges):
        outgoing.setdefault(e.u, []).append(pos)

    found: list[tuple[str, list[int], str]] = []

    def dfs(cur: str, used: int, path: list[int]) -> None:
        if len(found) >= limit:
            return
        if used == full:
            if cur != start:
                return
            circle = _barcode_from_path(path, edges, k)
            canon = canonical_circle(circle)
            if not any(c == canon for c, _, _ in found):
                order, rel = _align_order(path, circle, canon)
                found.append((canon, order, rel))
            return

        # 同标签平行边只尝试当前未用的最小位置 -> 重复读数互换不产生歧义，
        # 且使用次序按输入序号稳定分配（子列表按原序号升序构造）。
        choices: dict[str, int] = {}
        for pos in outgoing.get(cur, ()):
            if not (used >> pos) & 1:
                choices.setdefault(edges[pos].label, pos)

        can_finish = (full ^ used).bit_count() == 1
        for pos in choices.values():
            e = edges[pos]
            nused = used | (1 << pos)
            if not can_finish and not _remaining_connected(
                full ^ nused, edges, e.v, start
            ):
                continue
            dfs(e.v, nused, path + [pos])
            if len(found) >= limit:
                return

    dfs(start, 0, [])
    return found


def _build_degrees(edges: list[Edge]) -> tuple[dict[str, int], dict[str, int], list[str]]:
    indeg: dict[str, int] = {}
    outdeg: dict[str, int] = {}
    for e in edges:
        outdeg[e.u] = outdeg.get(e.u, 0) + 1
        indeg[e.v] = indeg.get(e.v, 0) + 1
    vertices = sorted(set(indeg) | set(outdeg))
    return indeg, outdeg, vertices


def _imbalance(
    vertices: list[str], indeg: dict[str, int], outdeg: dict[str, int]
) -> list[dict]:
    return [
        {"vertex": v, "in_degree": indeg.get(v, 0), "out_degree": outdeg.get(v, 0)}
        for v in vertices
        if indeg.get(v, 0) != outdeg.get(v, 0)
    ]


def _single_witness(
    canon: str, order_pos: list[int], edges: list[Edge], k: int, relation: str
) -> dict:
    """由"边列表位置"路径组装单条闭环见证（序号一律换算回原输入序号）。"""
    return {
        "canonical_barcode": canon,
        "barcode": _barcode_from_path(order_pos, edges, k),
        "canonical_relation": relation,
        "order": [edges[p].idx + 1 for p in order_pos],
        "evidence": _evidence(order_pos, edges, k),
    }


def _evidence(order: list[int], edges: list[Edge], k: int) -> list[dict]:
    """相邻（含首尾相接）重叠证据；order 为 edges 列表位置，输出原输入序号。"""
    out = []
    m = len(order)
    for pos, epos in enumerate(order):
        a = edges[epos]
        b = edges[order[(pos + 1) % m]]
        overlap = a.label[1:]
        assert overlap == b.label[:-1]
        out.append(
            {
                "position": pos,
                "prev": a.idx + 1,
                "next": b.idx + 1,
                "overlap": overlap,
                "overlap_length": k,
                "appended_base": b.label[-1],
            }
        )
    return out


def _assemble_single(seq: list[str], k: int, edges: list[Edge]) -> dict:
    indeg, outdeg, vertices = _build_degrees(edges)

    reasons: list[dict] = []

    imbalance = _imbalance(vertices, indeg, outdeg)
    if imbalance:
        reasons.append(
            {
                "code": "degree_imbalance",
                "message": "存在入度不等于出度的 k-mer，无法形成闭环",
                "vertices": imbalance,
            }
        )

    comps = _components(edges)
    if len(comps) > 1:
        reasons.append(
            {
                "code": "fragmented_graph",
                "message": "非零度顶点落在多个互不连通的分量中，无法拼成单一闭环",
                "component_count": len(comps),
                "components": comps,
            }
        )

    if reasons:
        return {
            "status": "no_solution",
            "overlap_length": k,
            "read_count": len(edges),
            "reasons": reasons,
        }

    start = vertices[0]  # 回路经过所有顶点，固定最小顶点作为枚举起点
    classes = _enumerate_classes(edges, start, k, limit=2)

    if len(classes) == 1:
        canon, order_pos, rel = classes[0]
        result = {
            "status": "unique",
            "overlap_length": k,
            "read_count": len(edges),
        }
        result.update(_single_witness(canon, order_pos, edges, k, rel))
        result["normalizations"] = ["cyclic_rotation", "reverse_complement"]
        return result

    if not classes:
        # 度数平衡且连通时必然存在欧拉回路；到不了这里只是防御性兜底。
        return {
            "status": "no_solution",
            "overlap_length": k,
            "read_count": len(edges),
            "reasons": [
                {
                    "code": "euler_search_failed",
                    "message": "图满足必要条件但未找到欧拉回路（内部错误）",
                }
            ],
        }

    witnesses = [_single_witness(c, o, edges, k, r) for c, o, r in classes]
    return {
        "status": "ambiguous",
        "overlap_length": k,
        "read_count": len(edges),
        "message": "存在多个不同规范等价类的闭环拼法，给出两条不同规范见证",
        "witnesses": witnesses,
        "normalizations": ["cyclic_rotation", "reverse_complement"],
    }


def _enumerate_partitions(n: int) -> "object":
    """枚举把 n 个带序号实例分成 A/B 两组的完整划分（位掩码表示 A 组）。

    读数 0 固定归 A，破除"两条条码不分先后"的组交换对称。
    采用 Gray 码顺序走查 2^(n-1) 种自由分配（边 1..n-1 的归属），相邻两次
    只翻转一条边，三元组依次产出 (A 组掩码, 本次翻转的边序号, 是否新进入 A)，
    配合调用方做 O(1) 增量度数维护。
    """
    free = n - 1
    total = 1 << free
    prev = 0
    for t in range(total):
        gray = t ^ (t >> 1)
        if t == 0:
            # 初始状态：只有边 0 在 A，无翻转（翻转边序号用 0 作哨兵）。
            yield 1, 0, True
            continue
        diff = gray ^ prev
        bit = diff.bit_length() - 1
        yield (gray << 1) | 1, bit + 1, bool((gray >> bit) & 1)
        prev = gray


def _group_state(edges: list[Edge], mask: int):
    """构造一组边（按原输入序号升序）并返回其度数、顶点、连通分量信息。"""
    group = [edges[i] for i in range(len(edges)) if (mask >> i) & 1]
    indeg, outdeg, vertices = _build_degrees(group)
    return group, indeg, outdeg, vertices


def _assemble_dual(seq: list[str], k: int, edges: list[Edge]) -> dict:
    n = len(edges)
    full = (1 << n) - 1
    # 完整划分规范等价类 -> 见证；逐划分增量构建，最多保留 2 个不同类。
    found: dict[tuple[str, str], dict] = {}

    # 顶点编号化，配合 Gray 码翻转做 O(1) 增量度数平衡维护。
    gin, gout, all_vertices = _build_degrees(edges)
    vid = {v: i for i, v in enumerate(all_vertices)}
    uids = [vid[e.u] for e in edges]
    vids = [vid[e.v] for e in edges]
    gbal = [0] * len(all_vertices)  # 全体（A∪B）每点出度-入度
    gdeg = [0] * len(all_vertices)
    for j in range(n):
        gbal[uids[j]] += 1
        gbal[vids[j]] -= 1
        gdeg[uids[j]] += 1
        gdeg[vids[j]] += 1

    # 全体已不平衡时，任何划分都不可能让两组同时平衡（两组平衡之和必为 0）。
    globally_balanced = all(b == 0 for b in gbal)

    # 组 A 的增量状态：初始 Gray 状态只有边 0 在 A。
    bal_a = [0] * len(all_vertices)
    deg_a = [0] * len(all_vertices)
    bal_a[uids[0]] += 1
    bal_a[vids[0]] -= 1
    deg_a[uids[0]] += 1
    deg_a[vids[0]] += 1

    def flip(j: int, into_a: bool) -> None:
        d = 1 if into_a else -1
        u, v = uids[j], vids[j]
        bal_a[u] += d
        bal_a[v] -= d
        deg_a[u] += d
        deg_a[v] += d

    def groups_balanced() -> bool:
        """组 A 增量已知；组 B 度数 = 全体 - A，逐点比较。"""
        for i, b in enumerate(bal_a):
            if b != 0:
                return False
            if gdeg[i] - deg_a[i] > 0 and gbal[i] - b != 0:
                return False
        return True

    # 标签多重集对相同的划分只可能产生同一批规范条码对（重复实例互换等价），
    # 记忆化跳过；所有标签互异时签名必然不同，直接省去这层开销。
    labels = [e.label for e in edges]
    memoize = len(set(labels)) < n
    seen_signatures: set[tuple[tuple[str, int], ...]] = set()
    label_counts: dict[str, int] = {labels[0]: 1}

    size_a = 1
    examined = 0
    eligible = 0
    first_fail: dict | None = None

    for mask, fj, added in _enumerate_partitions(n):
        examined += 1
        if fj:
            flip(fj, added)
            size_a += 1 if added else -1
            if memoize:
                lab = labels[fj]
                if added:
                    label_counts[lab] = label_counts.get(lab, 0) + 1
                else:
                    c = label_counts[lab] - 1
                    if c:
                        label_counts[lab] = c
                    else:
                        del label_counts[lab]

        if not (MIN_GROUP_READS <= size_a <= n - MIN_GROUP_READS):
            continue
        if not globally_balanced:
            continue
        eligible += 1

        if memoize:
            signature = tuple(sorted(label_counts.items()))
            if signature in seen_signatures:
                continue
            seen_signatures.add(signature)

        if not groups_balanced():
            continue

        other = full ^ mask
        ga, ia, oa, va = _group_state(edges, mask)
        gb, ib, ob, vb = _group_state(edges, other)
        # 平衡是欧拉回路的必要条件；再确认两组各自弱连通（单一闭环）。
        comps_a = _components(ga)
        comps_b = _components(gb)
        if len(comps_a) > 1 or len(comps_b) > 1:
            if first_fail is None:
                first_fail = {
                    "group_a_reads": [e.idx + 1 for e in ga],
                    "group_b_reads": [e.idx + 1 for e in gb],
                    "imbalance_a": [],
                    "imbalance_b": [],
                    "components_a": len(comps_a),
                    "components_b": len(comps_b),
                }
            continue

        # 两组各自独立枚举闭环（同株内仍按原重叠规则与规范归一裁决）。
        classes_a = _enumerate_classes(ga, va[0], k, limit=2)
        classes_b = _enumerate_classes(gb, vb[0], k, limit=2)
        if not classes_a or not classes_b:
            continue  # 防御性：平衡连通时理论上不会发生

        # 组内若本身存在多个规范类，所有跨组配对都要参与全局等价类判定
        #（不同配对即不同的完整双组结论）。
        stop = False
        for ca, oa0, ra in classes_a:
            for cb, ob0, rb in classes_b:
                key = tuple(sorted((ca, cb)))
                if key in found:
                    continue
                wa = _single_witness(ca, oa0, ga, k, ra)
                wb = _single_witness(cb, ob0, gb, k, rb)
                # 两条条码不分先后：按规范条码排序两组；条码相同则按组内
                # 最小输入序号稳定排序（重复实例互换不制造歧义）。
                groups = [wa, wb]
                groups.sort(key=lambda w: (w["canonical_barcode"], min(w["order"])))
                assignment = [0] * n
                for rank, w in enumerate(groups):
                    for one in w["order"]:
                        assignment[one - 1] = rank + 1
                found[key] = {
                    "barcodes": [w["canonical_barcode"] for w in groups],
                    "groups": groups,
                    "assignment": assignment,
                }
                if len(found) >= 2:
                    stop = True
                    break
            if stop:
                break
        if stop:
            break

    base = {
        "overlap_length": k,
        "read_count": n,
        "barcode_count": 2,
        "normalizations": [
            "cyclic_rotation",
            "reverse_complement",
            "group_swap",
            "duplicate_read_interchange",
        ],
    }

    if len(found) == 1:
        witness = next(iter(found.values()))
        return {
            "status": "unique",
            **base,
            "message": "全部读数可唯一地划分为两组，各自形成符合重叠规则的闭环",
            "canonical_barcodes": witness["barcodes"],
            "groups": witness["groups"],
            "assignment": witness["assignment"],
        }

    if len(found) >= 2:
        witnesses = []
        for wit in found.values():
            witnesses.append(
                {
                    "canonical_barcodes": wit["barcodes"],
                    "groups": wit["groups"],
                    "assignment": wit["assignment"],
                }
            )
        return {
            "status": "ambiguous",
            **base,
            "message": "存在两份不同的完整双组划分（规范条码对不同），各给一份见证",
            "witnesses": witnesses,
        }

    reason = {
        "code": "no_dual_partition",
        "message": (
            "联合枚举了每组至少 "
            f"{MIN_GROUP_READS} 条、至多 {n - MIN_GROUP_READS} 条的完整划分"
            "（读数 1 固定归第一组以消除组交换对称，相邻划分按 Gray 码只改派一条"
            "读数并增量核对度数）；没有任何划分能让两组边子图同时平衡且弱连通，"
            "即无法把全部读数拆成两个各用一次的合格闭环"
        ),
        "partitions_examined": examined,
        "size_eligible_partitions": eligible,
        "min_reads_per_barcode": MIN_GROUP_READS,
    }
    if not globally_balanced:
        reason["global_degree_imbalance"] = _imbalance(all_vertices, gin, gout)
    if first_fail is not None:
        reason["example_failed_partition"] = first_fail
    return {
        "status": "no_solution",
        **base,
        "reasons": [reason],
    }


def assemble(sequences: list[str], barcode_count: int = 1) -> dict:
    """主入口：返回唯一 / 歧义 / 无解三类结果之一。

    barcode_count=1（默认，兼容原行为）：全部读数拼成单一闭环。
    barcode_count=2：全部读数恰好分属两条环状条码，联合划分、分别成环。
    """
    seq = list(sequences)
    k = validate(seq, barcode_count)

    edges = [
        Edge(s[:-1], s[1:], s, i)
        for i, s in enumerate(seq)
    ]

    if barcode_count == 2:
        return _assemble_dual(seq, k, edges)
    return _assemble_single(seq, k, edges)