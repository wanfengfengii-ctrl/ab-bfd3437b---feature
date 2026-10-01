"""环状 DNA 条码拼接核心算法。

把每条长度 L = k+1 的读数看成德布鲁因多重图中的一条有向边
    s[:-1]  ->  s[1:]
（序列每次出现都是一条独立平行边，即一份证据）。
"每条序列恰好使用一次、相邻重叠 k = L-1 的闭环拼接" 等价于该多重图的
欧拉回路。

等价归一化（用于判断是否为同一个可信条码）：
  * 环的循环移位；
  * 整条环的反向互补；
  * 完全相同的重复读数（平行边）互换不算新答案。

枚举在固定起点上进行带剪枝的欧拉回路 DFS，同一顶点上标签相同的候选边
只取当前未用的最小输入序号一条（稳定分配），每一步都丢掉"跨过桥"的走法，
找到 2 个不同规范等价类即可判定为歧义。

双株模式（assemble(..., barcode_count=2)，见本模块下半部分）：联合枚举
全部读数的二分划分，两组各自必须是平衡且弱连通的欧拉子多重图（各至少
3 条），以排序后的规范条码对为完整划分的规范等价类键。
"""
from __future__ import annotations

from dataclasses import dataclass

ALPHABET = frozenset("ACGT")
_COMP = str.maketrans("ACGT", "TGCA")

MIN_READS, MAX_READS = 6, 24
MIN_READS_DUAL, MAX_READS_DUAL = 6, 18
MIN_LEN, MAX_LEN = 3, 8
MIN_GROUP_READS = 3  # 双株模式下每条条码至少使用的读数条数


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
    """校验输入，返回 k = L-1；不合法抛 ValidationError。

    barcode_count=2 为双株模式：6~18 条读数，联合划分为两个各至少 3 条的
    闭环。
    """
    if not isinstance(sequences, list) or not sequences:
        raise ValidationError("sequences 必须是非空数组")
    lo, hi = (MIN_READS_DUAL, MAX_READS_DUAL) if barcode_count == 2 else (MIN_READS, MAX_READS)
    if not (lo <= len(sequences) <= hi):
        raise ValidationError(
            f"序列条数必须在 {lo}~{hi} 之间，实际 {len(sequences)} 条"
        )
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
    """枚举至多 limit 个不同规范等价类，返回 (规范条码, 边序号路径, 归一关系)。"""
    m = len(edges)
    full = (1 << m) - 1

    outgoing: dict[str, list[int]] = {}
    for e in edges:
        outgoing.setdefault(e.u, []).append(e.idx)

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

        # 同标签平行边只尝试当前未用的最小序号 -> 重复读数互换不产生歧义，
        # 且使用次序按输入序号稳定分配。
        choices: dict[str, int] = {}
        for eid in outgoing.get(cur, ()):
            if not (used >> eid) & 1:
                choices.setdefault(edges[eid].label, eid)

        can_finish = (full ^ used).bit_count() == 1
        for eid in choices.values():
            e = edges[eid]
            nused = used | (1 << eid)
            if not can_finish and not _remaining_connected(
                full ^ nused, edges, e.v, start
            ):
                continue
            dfs(e.v, nused, path + [eid])
            if len(found) >= limit:
                return

    dfs(start, 0, [])
    return found


def _evidence(order: list[int], edges: list[Edge], k: int) -> list[dict]:
    """相邻（含首尾相接）重叠证据。"""
    out = []
    m = len(order)
    for pos, eid in enumerate(order):
        a = edges[eid]
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


# ---------------------------------------------------------------------------
# 双株模式：全部读数恰好划分为两个子多重图，各自构成一条独立闭环。
#
# 一个划分 (A, B) 合格的必要条件（欧拉定理）：
#   1. |A|, |B| >= MIN_GROUP_READS；
#   2. 每个子图内入度 = 出度（平衡）；
#   3. 每个子图的非零度顶点弱连通（恰为一个欧拉分量）。
# 划分的联合枚举保证"先拼一条再处理余料"与一次性联合判定同解集，
# 这里直接联合枚举、同时验证两组。
#
# 划分等价类归一化：
#   * 组不分先后（{A,B} 与 {B,A} 同一类）——以排序后的规范条码对为键；
#   * 组内环的循环移位 / 整条反向互补——canonical_circle；
#   * 相同读数（平行边）实例互换——组内取最小输入序号的实例，
#     实例的不同分配不产生新的规范等价类。
# ---------------------------------------------------------------------------


def _dual_balanced(mask: int, edges: list[Edge]) -> bool:
    """mask 选出的子图是否每个顶点入度 = 出度（空集视为平衡）。"""
    delta: dict[str, int] = {}
    m = mask
    while m:
        b = m & -m
        e = edges[b.bit_length() - 1]
        delta[e.u] = delta.get(e.u, 0) + 1
        delta[e.v] = delta.get(e.v, 0) - 1
        m ^= b
    return all(d == 0 for d in delta.values())


def _dual_connected(mask: int, edges: list[Edge]) -> bool:
    """mask 选出的子图非零度顶点是否落在同一个弱连通分量内。"""
    parent: dict[str, str] = {}

    def find(x: str) -> str:
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    anchor = None
    m = mask
    while m:
        b = m & -m
        e = edges[b.bit_length() - 1]
        ru, rv = find(e.u), find(e.v)
        if ru != rv:
            parent[rv] = ru
        if anchor is None:
            anchor = e.u
        m ^= b
    if anchor is None:
        return False
    root = find(anchor)
    return all(find(v) == root for v in parent)


def _canonicalize_mask(mask: int, edges: list[Edge]) -> int:
    """把划分中的重复读数实例按输入序号稳定归一。

    对每个读数标签，组 A 若需要 a 个实例，则取该标签序号最小的 a 条；
    其余归组 B。标签多重集不变，故平衡/连通性与所有组内闭环不变，
    而不同实例分配（重复平行边互换）收敛到同一个规范掩码。
    """
    by_label: dict[str, list[int]] = {}
    need: dict[str, int] = {}
    m = mask
    while m:
        b = m & -m
        i = b.bit_length() - 1
        lab = edges[i].label
        by_label.setdefault(lab, []).append(i)
        need[lab] = need.get(lab, 0) + 1
        m ^= b
    out = 0
    for lab, idxs in by_label.items():  # dict 按首次出现有序，结果确定
        idxs.sort()
        for i in idxs[: need[lab]]:
            out |= 1 << i
    return out


def _split_subgraphs(
    mask: int, edges: list[Edge]
) -> tuple[list[Edge], list[Edge]]:
    """按掩码拆成两张子多重图，子图内边重新局部编号。"""
    sub_a: list[Edge] = []
    sub_b: list[Edge] = []
    for e in edges:
        if (mask >> e.idx) & 1:
            sub_a.append(Edge(e.u, e.v, e.label, len(sub_a)))
        else:
            sub_b.append(Edge(e.u, e.v, e.label, len(sub_b)))
    return sub_a, sub_b


def _subgraph_class_paths(
    sub: list[Edge], k: int, limit: int = 2
) -> list[tuple[str, list[str], str]]:
    """子多重图的至多 limit 个规范等价类。

    返回 (规范条码, 回路读数标签序列, 归一关系)；以标签序列而非实例序号
    表达，天然对重复实例互换不变。
    """
    start = min({e.u for e in sub} | {e.v for e in sub})
    return [
        (canon, [sub[j].label for j in order], rel)
        for canon, order, rel in _enumerate_classes(sub, start, k, limit=limit)
    ]


def _enumerate_dual_partitions(
    edges: list[Edge], k: int, limit: int = 2
) -> tuple[dict, dict]:
    """枚举至多 limit 个不同的规范双组划分等价类。

    键为排序后的 (规范条码 A, 规范条码 B)；值为
    (规范掩码, 组0见证三元组, 组1见证三元组)。联合枚举所有划分，
    两组的平衡与连通同时判定，不存在"先拼一条再处理余料"。
    """
    n = len(edges)
    full = (1 << n) - 1
    min_size = MIN_GROUP_READS
    max_size = n - min_size

    # 键 -> (规范掩码, (canon, labels, rel) 组0, (canon, labels, rel) 组1)
    found: dict[tuple[str, str], tuple] = {}
    seen_norms: set[int] = set()
    stats = {"partitions_checked": 0, "valid_partitions": 0}

    def record(mask: int) -> None:
        norm = _canonicalize_mask(mask, edges)
        if norm in seen_norms:  # 同标签多重集的实例置换已归一，不重复取证
            return
        seen_norms.add(norm)
        sub_a, sub_b = _split_subgraphs(norm, edges)
        classes_a = _subgraph_class_paths(sub_a, k, limit=2)
        classes_b = _subgraph_class_paths(sub_b, k, limit=2)
        if not classes_a or not classes_b:
            return
        for ca, la, ra in classes_a:
            for cb, lb, rb in classes_b:
                key = tuple(sorted((ca, cb)))
                if key not in found:
                    wa = (ca, la, ra)
                    wb = (cb, lb, rb)
                    found[key] = (norm, wa, wb)

    # 固定边 0 属于组 A，消去组交换对称（组不分先后）。
    # n <= 18，叶节点至多 2^17，配合尺寸剪枝毫秒级完成。
    def grow(i: int, mask: int, size: int) -> None:
        if len(found) >= limit:
            return
        if size > max_size:
            return
        if size + (n - i) < min_size:  # 剩余边全给 A 也不够最小尺寸
            return
        if i == n:
            stats["partitions_checked"] += 1
            other = full ^ mask
            if not (_dual_balanced(mask, edges) and _dual_connected(mask, edges)):
                return
            if not (_dual_balanced(other, edges) and _dual_connected(other, edges)):
                return
            stats["valid_partitions"] += 1
            record(mask)
            return

        grow(i + 1, mask | (1 << i), size + 1)  # 边 i 归 A
        if len(found) >= limit:
            return
        grow(i + 1, mask, size)  # 边 i 归 B

    grow(1, 1, 1)
    return found, stats


def _circle_from_labels(labels: list[str]) -> str:
    s = labels[0] + "".join(t[-1] for t in labels[1:])
    return s[: len(labels)]


def _assign_instances(
    labels: list[str], mask: int, edges: list[Edge], in_group_a: bool
) -> list[int]:
    """把回路标签序列稳定映射到全局输入序号：同标签取组内最小的可用序号。"""
    pool: dict[str, list[int]] = {}
    for e in edges:
        in_a = (mask >> e.idx) & 1
        if in_a == in_group_a:
            pool.setdefault(e.label, []).append(e.idx)
    for idxs in pool.values():
        idxs.sort()  # 本身即有序，显式排序表明按输入序号稳定分配
    cursor: dict[str, int] = {}
    order: list[int] = []
    for lab in labels:
        j = cursor.get(lab, 0)
        order.append(pool[lab][j])
        cursor[lab] = j + 1
    return order


def _dual_witness(entry: tuple, edges: list[Edge], k: int) -> dict:
    """把一条规范划分记录展开为完整双组见证。"""
    mask, wa, wb = entry

    def build(
        trip: tuple[str, list[str], str], in_group_a: bool
    ) -> tuple[dict, str, list[int]]:
        canon, labels, rel = trip
        order0 = _assign_instances(labels, mask, edges, in_group_a)
        return (
            {
                "canonical_barcode": canon,
                "barcode": _circle_from_labels(labels),
                "canonical_relation": rel,
                "order": [i + 1 for i in order0],
                "evidence": _evidence(order0, edges, k),
            },
            canon,
            order0,
        )

    g0, c0, o0 = build(wa, True)
    g1, c1, o1 = build(wb, False)
    # 按规范条码排序两组；同码时按稳定输入序号序列排序，输出确定
    if (c0, o0) > (c1, o1):
        g0, g1 = g1, g0
        c0, c1 = c1, c0
        o0, o1 = o1, o0

    group_of: dict[int, int] = {}
    for eid in o0:
        group_of[eid] = 0
    for eid in o1:
        group_of[eid] = 1
    return {
        "barcodes": [g0, g1],
        "canonical_barcodes": [c0, c1],
        "assignment": [
            {"read_index": i + 1, "group": group_of[i]} for i in range(len(edges))
        ],
    }


def _global_imbalance(edges: list[Edge]) -> list[dict]:
    indeg: dict[str, int] = {}
    outdeg: dict[str, int] = {}
    for e in edges:
        outdeg[e.u] = outdeg.get(e.u, 0) + 1
        indeg[e.v] = indeg.get(e.v, 0) + 1
    return [
        {"vertex": v, "in_degree": indeg.get(v, 0), "out_degree": outdeg.get(v, 0)}
        for v in sorted(set(indeg) | set(outdeg))
        if indeg.get(v, 0) != outdeg.get(v, 0)
    ]


_DUAL_NORMALIZATIONS = [
    "cyclic_rotation",
    "reverse_complement",
    "group_swap",
    "duplicate_read_interchange",
]


def assemble_dual(sequences: list[str]) -> dict:
    """双株模式入口：把全部读数联合划分为两个各自合格的闭环。"""
    seq = list(sequences)
    k = validate(seq, barcode_count=2)
    edges = [Edge(s[:-1], s[1:], s, i) for i, s in enumerate(seq)]
    n = len(edges)

    imbalance = _global_imbalance(edges)
    if imbalance:
        return {
            "status": "no_solution",
            "barcode_count": 2,
            "overlap_length": k,
            "read_count": n,
            "reasons": [
                {
                    "code": "degree_imbalance",
                    "message": "全体读数在某些 k-mer 上入度不等于出度，"
                               "任意划分都无法让两组同时形成闭环",
                    "vertices": imbalance,
                }
            ],
        }

    found, stats = _enumerate_dual_partitions(edges, k, limit=2)

    if not found:
        comps = _components(edges)
        reason = {
            "code": "no_dual_partition",
            "message": (
                "不存在把全部读数恰好分成两组的完整划分，使两组各自"
                "构成度数平衡且弱连通的闭合多重图（每个合格闭环"
                f"至少 {MIN_GROUP_READS} 条读数，每条读数恰好使用一次）"
            ),
            "constraints": {
                "reads_used_exactly_once": True,
                "min_reads_per_barcode": MIN_GROUP_READS,
                "group_eulerian": ["balanced_degree", "weakly_connected"],
                "group_size_range": [MIN_GROUP_READS, n - MIN_GROUP_READS],
            },
            "diagnostics": {
                "weak_component_count": len(comps),
                "weak_component_sizes": [len(c) for c in comps],
                "weak_components": comps,
                "note": (
                    "弱连通分量多于两个时无法只分成两个连通组；恰有两个时"
                    "每个分量必须独立成环且各含足够读数"
                ),
            },
            "enumeration": stats,
        }
        return {
            "status": "no_solution",
            "barcode_count": 2,
            "overlap_length": k,
            "read_count": n,
            "reasons": [reason],
        }

    if len(found) == 1:
        entry = next(iter(found.values()))
        w = _dual_witness(entry, edges, k)
        return {
            "status": "unique",
            "barcode_count": 2,
            "overlap_length": k,
            "read_count": n,
            "canonical_barcodes": w["canonical_barcodes"],
            "barcodes": w["barcodes"],
            "assignment": w["assignment"],
            "normalizations": list(_DUAL_NORMALIZATIONS),
        }

    witnesses = [_dual_witness(entry, edges, k) for entry in found.values()]
    return {
        "status": "ambiguous",
        "barcode_count": 2,
        "overlap_length": k,
        "read_count": n,
        "message": "存在两个以上规范等价类的完整双组划分，给出两份不同的完整双组见证",
        "witnesses": witnesses,
        "normalizations": list(_DUAL_NORMALIZATIONS),
    }


def assemble(sequences: list[str], barcode_count: int = 1) -> dict:
    """主入口：单株（默认，原行为）或双株（barcode_count=2）。"""
    if barcode_count == 2:
        return assemble_dual(sequences)
    if barcode_count != 1:
        raise ValidationError("barcode_count 只支持 1（默认）或 2")
    seq = list(sequences)
    k = validate(seq)

    edges = [
        Edge(s[:-1], s[1:], s, i)
        for i, s in enumerate(seq)
    ]

    indeg: dict[str, int] = {}
    outdeg: dict[str, int] = {}
    for e in edges:
        outdeg[e.u] = outdeg.get(e.u, 0) + 1
        indeg[e.v] = indeg.get(e.v, 0) + 1
    vertices = sorted(set(indeg) | set(outdeg))

    reasons: list[dict] = []

    imbalance = [
        {"vertex": v, "in_degree": indeg.get(v, 0), "out_degree": outdeg.get(v, 0)}
        for v in vertices
        if indeg.get(v, 0) != outdeg.get(v, 0)
    ]
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

    def witness(canon: str, order0: list[int], relation: str) -> dict:
        return {
            "canonical_barcode": canon,
            "barcode": _barcode_from_path(order0, edges, k),
            "canonical_relation": relation,
            "order": [i + 1 for i in order0],
            "evidence": _evidence(order0, edges, k),
        }

    if len(classes) == 1:
        canon, order0, rel = classes[0]
        return {
            "status": "unique",
            "overlap_length": k,
            "read_count": len(edges),
            "canonical_barcode": canon,
            "barcode": _barcode_from_path(order0, edges, k),
            "canonical_relation": rel,
            "order": [i + 1 for i in order0],
            "evidence": _evidence(order0, edges, k),
            "normalizations": ["cyclic_rotation", "reverse_complement"],
        }

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

    witnesses = [witness(c, o, r) for c, o, r in classes]
    return {
        "status": "ambiguous",
        "overlap_length": k,
        "read_count": len(edges),
        "message": "存在多个不同规范等价类的闭环拼法，给出两条不同规范见证",
        "witnesses": witnesses,
        "normalizations": ["cyclic_rotation", "reverse_complement"],
    }