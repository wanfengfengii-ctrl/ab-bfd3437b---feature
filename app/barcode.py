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
"""
from __future__ import annotations

from dataclasses import dataclass

ALPHABET = frozenset("ACGT")
_COMP = str.maketrans("ACGT", "TGCA")

MIN_READS, MAX_READS = 6, 24
MIN_LEN, MAX_LEN = 3, 8


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


def validate(sequences: list[str]) -> int:
    """校验输入，返回 k = L-1；不合法抛 ValidationError。"""
    if not isinstance(sequences, list) or not sequences:
        raise ValidationError("sequences 必须是非空数组")
    if not (MIN_READS <= len(sequences) <= MAX_READS):
        raise ValidationError(
            f"序列条数必须在 {MIN_READS}~{MAX_READS} 之间，实际 {len(sequences)} 条"
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


def assemble(sequences: list[str]) -> dict:
    """主入口：返回唯一 / 歧义 / 无解三类结果之一。"""
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