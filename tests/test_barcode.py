"""barcode 核心算法单元测试 + 独立暴力参照交叉验证。"""
import itertools
import os
import random
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "app"))

from barcode import (  # noqa: E402
    canonical_circle,
    reverse_complement,
    ValidationError,
    assemble,
)


def reads_of(circle: str, length: int) -> list[str]:
    """从环状串生成每个位置的定长读数（允许 length 大于环长）。"""
    n = len(circle)
    doubled = circle * (length // n + 2)
    return [doubled[i:i + length] for i in range(n)]


def brute_classes(reads: list[str]) -> set[str]:
    """独立参照：枚举边实例的全部排列，收集不同规范环。"""
    n = len(reads)
    classes = set()
    for perm in itertools.permutations(range(n)):
        ok = True
        for a, b in zip(perm, perm[1:]):
            if reads[a][1:] != reads[b][:-1]:
                ok = False
                break
        if ok and reads[perm[-1]][1:] != reads[perm[0]][:-1]:
            ok = False
        if ok:
            circle = reads[perm[0]] + "".join(reads[i][-1] for i in perm[1:])
            classes.add(canonical_circle(circle[:n]))
    return classes


def check_witness(test: unittest.TestCase, w: dict, reads: list[str], k: int) -> None:
    """见证内部一致性：次序、重叠证据、条码窗口、规范代表。"""
    n = len(reads)
    length = k + 1
    test.assertEqual(sorted(w["order"]), list(range(1, n + 1)))
    test.assertEqual(len(w["evidence"]), n)
    reps = (n + length) // n + 1  # 读数可能比环长，需要足够多圈
    doubled = w["barcode"] * reps
    for pos, ev in enumerate(w["evidence"]):
        test.assertEqual(ev["position"], pos)
        a = reads[ev["prev"] - 1]
        b = reads[ev["next"] - 1]
        test.assertEqual(ev["prev"], w["order"][pos])
        test.assertEqual(ev["next"], w["order"][(pos + 1) % n])
        test.assertEqual(a[1:], b[:-1])
        test.assertEqual(ev["overlap"], a[1:])
        test.assertEqual(ev["overlap_length"], k)
        test.assertEqual(ev["appended_base"], b[-1])
        test.assertEqual(doubled[pos:pos + length], reads[w["order"][pos] - 1])
    test.assertEqual(canonical_circle(w["barcode"]), w["canonical_barcode"])


class CanonicalTests(unittest.TestCase):
    def test_rotation_and_reverse_complement(self):
        s = "AATCGCAGT"
        for r in range(len(s)):
            rot = s[r:] + s[:r]
            self.assertEqual(canonical_circle(rot), canonical_circle(s))
        self.assertEqual(canonical_circle(reverse_complement(s)), canonical_circle(s))
        self.assertEqual(reverse_complement(reverse_complement(s)), s)
        self.assertEqual(canonical_circle("AACAAC"), "AACAAC")


class ValidationTests(unittest.TestCase):
    def test_count_bounds(self):
        with self.assertRaises(ValidationError):
            assemble(["AAA"] * 5)
        with self.assertRaises(ValidationError):
            assemble(["AAA"] * 25)

    def test_length_bounds(self):
        with self.assertRaises(ValidationError):
            assemble(["AA"] * 6)
        with self.assertRaises(ValidationError):
            assemble(["AAAAAAAAA"] * 6)

    def test_equal_length_and_alphabet(self):
        with self.assertRaises(ValidationError):
            assemble(["AAA", "AAC", "ACA", "CAA", "ATT", "TT"])  # 不等长
        with self.assertRaises(ValidationError):
            assemble(["AAA", "AAC", "ACA", "CAN", "AAA", "AAT"])  # 非法碱基 N

    def test_lowercase_is_normalized(self):
        r = assemble(["aaa", "aac", "aca", "caa", "aat", "ata"])
        self.assertNotEqual(r["status"], "invalid_input")


class UniqueTests(unittest.TestCase):
    READS = ["AAT", "ATC", "TCG", "CGC", "GCA", "CAA"]  # 环 AATCGC

    def test_unique_basic(self):
        r = assemble(list(self.READS))
        self.assertEqual(r["status"], "unique")
        self.assertEqual(r["canonical_barcode"], "AATCGC")
        self.assertEqual(r["overlap_length"], 2)
        check_witness(self, r, self.READS, 2)

    def test_reverse_complement_input_same_class(self):
        rc_reads = [reverse_complement(s) for s in self.READS]
        r1 = assemble(list(self.READS))
        r2 = assemble(rc_reads)
        self.assertEqual(r2["status"], "unique")
        self.assertEqual(r1["canonical_barcode"], r2["canonical_barcode"])
        check_witness(self, r2, rc_reads, 2)

    def test_shuffled_input_same_barcode_stable_order(self):
        shuffled = ["TCG", "CAA", "AAT", "GCA", "ATC", "CGC"]
        r = assemble(list(shuffled))
        self.assertEqual(r["status"], "unique")
        self.assertEqual(r["canonical_barcode"], "AATCGC")
        # 按输入序号稳定分配：起点边 AAT 在打乱后的序号为 3
        self.assertEqual(r["order"][0], 3)
        check_witness(self, r, shuffled, 2)
        # 重复运行结果一致
        self.assertEqual(assemble(list(shuffled))["order"], r["order"])

    def test_duplicate_reads_interchange_not_ambiguous(self):
        # 周期环 AAC：每个读数出现两次，平行边互换只能得到一个规范类
        reads = ["AAC", "ACA", "CAA", "AAC", "ACA", "CAA"]
        r = assemble(list(reads))
        self.assertEqual(r["status"], "unique")
        self.assertEqual(r["canonical_barcode"], "AACAAC")
        self.assertEqual(r["order"], [1, 2, 3, 4, 5, 6])  # 最小可用序号
        check_witness(self, r, reads, 2)

        shuffled = ["ACA", "CAA", "AAC", "AAC", "ACA", "CAA"]
        r2 = assemble(list(shuffled))
        self.assertEqual(r2["status"], "unique")
        self.assertEqual(r2["canonical_barcode"], "AACAAC")
        self.assertEqual(r2["order"], [3, 1, 2, 4, 5, 6])
        check_witness(self, r2, shuffled, 2)

    def test_length_8_small_circle(self):
        # 环长 6 小于读数长度 8：重复结构仍可闭环
        reads = reads_of("AAAAAC", 8)
        self.assertEqual(len(reads), 6)
        r = assemble(list(reads))
        self.assertEqual(r["status"], "unique")
        self.assertEqual(r["canonical_barcode"], "AAAAAC")
        check_witness(self, r, reads, 7)


class AmbiguousTests(unittest.TestCase):
    # θ 图：关节点 AA 与 CC 之间各有两条不同路径，去程/回程的交错方式
    # （P+R、Q+S 与 P+S、Q+R）产生不能经旋转或反向互补归一的环。
    #   P: AA→AC→CC   AAC, ACC        R: CC→CA→AA  CCA, CAA
    #   Q: AA→AG→GC→CC AAG, AGC, GCC  S: CC→CG→GA→AA CCG, CGA, GAA
    READS = ["AAC", "ACC", "AAG", "AGC", "GCC",
             "CCA", "CAA", "CCG", "CGA", "GAA"]

    def test_two_distinct_canonical_witnesses(self):
        r = assemble(list(self.READS))
        self.assertEqual(r["status"], "ambiguous")
        self.assertEqual(len(r["witnesses"]), 2)
        barcodes = {w["canonical_barcode"] for w in r["witnesses"]}
        self.assertEqual(len(barcodes), 2)
        for w in r["witnesses"]:
            check_witness(self, w, self.READS, 2)

    def test_matches_brute_force(self):
        expected = brute_classes(self.READS)
        self.assertGreaterEqual(len(expected), 2)
        got = {w["canonical_barcode"] for w in assemble(list(self.READS))["witnesses"]}
        self.assertTrue(got <= expected)


class NoSolutionTests(unittest.TestCase):
    def test_degree_imbalance(self):
        reads = ["AAA", "AAA", "AAC", "ACA", "CAA", "CAA"]
        r = assemble(list(reads))
        self.assertEqual(r["status"], "no_solution")
        codes = {x["code"] for x in r["reasons"]}
        self.assertIn("degree_imbalance", codes)
        self.assertNotIn("fragmented_graph", codes)
        imba = next(x for x in r["reasons"] if x["code"] == "degree_imbalance")
        bad = {v["vertex"] for v in imba["vertices"]}
        self.assertEqual(bad, {"AA", "CA"})

    def test_fragmented_graph(self):
        reads = ["AAA", "AAA", "CCC", "CCC", "GGG", "GGG"]
        r = assemble(list(reads))
        self.assertEqual(r["status"], "no_solution")
        codes = {x["code"] for x in r["reasons"]}
        self.assertIn("fragmented_graph", codes)
        self.assertNotIn("degree_imbalance", codes)
        frag = next(x for x in r["reasons"] if x["code"] == "fragmented_graph")
        self.assertEqual(frag["component_count"], 3)
        self.assertEqual(
            sorted(c[0] for c in frag["components"]), ["AA", "CC", "GG"]
        )
        self.assertTrue(all(len(c) == 1 for c in frag["components"]))

    def test_imbalance_and_fragmentation_reported_together(self):
        reads = ["AAA", "AAA", "AAC", "CCC", "GGG", "TTT"]
        r = assemble(list(reads))
        self.assertEqual(r["status"], "no_solution")
        codes = {x["code"] for x in r["reasons"]}
        self.assertEqual(codes, {"degree_imbalance", "fragmented_graph"})
        frag = next(x for x in r["reasons"] if x["code"] == "fragmented_graph")
        # 分量 {AA, AC}（含两个顶点）以及 {CC}、{GG}、{TT}
        self.assertEqual(frag["component_count"], 4)
        self.assertEqual(
            sorted(frag["components"]),
            [["AA", "AC"], ["CC"], ["GG"], ["TT"]],
        )


class BruteForceCrossCheck(unittest.TestCase):
    """随机环状游走生成的小图：求解器结论必须与全排列暴力枚举一致。"""

    def test_random_circles(self):
        rng = random.Random(20260929)
        for _ in range(60):
            n = rng.randint(6, 7)
            circle = "".join(rng.choice("AC") for _ in range(n))
            reads = reads_of(circle, 3)
            expected = brute_classes(reads)
            r = assemble(list(reads))
            if len(expected) == 1:
                self.assertEqual(r["status"], "unique", reads)
                self.assertEqual(r["canonical_barcode"], next(iter(expected)))
                check_witness(self, r, reads, 2)
            else:
                self.assertEqual(r["status"], "ambiguous", reads)
                got = {w["canonical_barcode"] for w in r["witnesses"]}
                self.assertEqual(len(got), 2)
                self.assertTrue(got <= expected)


# ---------------------------------------------------------------------------
# 双株模式
# ---------------------------------------------------------------------------


def _circles_for(sub: list[str]) -> set[str]:
    """一组读数能拼出的全部不同规范环（重叠约束下的 DFS 暴力）。"""
    m = len(sub)
    classes: set[str] = set()

    def dfs(path: list[int], used: int) -> None:
        if len(path) == m:
            first, last = sub[path[0]], sub[path[-1]]
            if last[1:] == first[:-1]:
                circle = first + "".join(sub[i][-1] for i in path[1:])
                classes.add(canonical_circle(circle[:m]))
            return
        last = sub[path[-1]]
        for j in range(m):
            if not (used >> j) & 1 and last[1:] == sub[j][:-1]:
                dfs(path + [j], used | (1 << j))

    for s in range(m):
        dfs([s], 1 << s)
    return classes


def brute_dual_classes(reads: list[str], min_group: int = 3) -> set[tuple[str, str]]:
    """独立参照：枚举所有二分划分（固定读数 0 归 A 消去组交换），
    两组各自暴力求环，收集排序后的规范条码对。"""
    n = len(reads)
    classes: set[tuple[str, str]] = set()
    for mask in range(1, 1 << n):
        if not (mask & 1):
            continue
        if not (min_group <= mask.bit_count() <= n - min_group):
            continue
        a = [reads[i] for i in range(n) if (mask >> i) & 1]
        b = [reads[i] for i in range(n) if not ((mask >> i) & 1)]
        ca, cb = _circles_for(a), _circles_for(b)
        if ca and cb:
            for x in ca:
                for y in cb:
                    classes.add(tuple(sorted((x, y))))
    return classes


def check_dual_witness(
    test: unittest.TestCase, w: dict, reads: list[str], k: int
) -> None:
    """一份完整双组见证的内部一致性。"""
    n = len(reads)
    length = k + 1
    test.assertEqual(len(w["barcodes"]), 2)
    all_orders: list[int] = []
    for g in w["barcodes"]:
        all_orders += g["order"]
        size = len(g["order"])
        test.assertGreaterEqual(size, 3)
        test.assertEqual(len(g["evidence"]), size)
        reps = (size + length) // size + 1
        doubled = g["barcode"] * reps
        for pos, ev in enumerate(g["evidence"]):
            test.assertEqual(ev["position"], pos)
            test.assertEqual(ev["prev"], g["order"][pos])
            test.assertEqual(ev["next"], g["order"][(pos + 1) % size])
            a = reads[ev["prev"] - 1]
            b = reads[ev["next"] - 1]
            test.assertEqual(a[1:], b[:-1])
            test.assertEqual(ev["overlap"], a[1:])
            test.assertEqual(ev["overlap_length"], k)
            test.assertEqual(ev["appended_base"], b[-1])
            test.assertEqual(doubled[pos:pos + length], reads[ev["prev"] - 1])
        test.assertEqual(canonical_circle(g["barcode"]), g["canonical_barcode"])
    test.assertEqual(sorted(all_orders), list(range(1, n + 1)))
    test.assertEqual(
        [a["read_index"] for a in w["assignment"]], list(range(1, n + 1))
    )
    group_of = {}
    for gi, g in enumerate(w["barcodes"]):
        for idx in g["order"]:
            group_of[idx] = gi
    test.assertEqual(
        [a["group"] for a in w["assignment"]],
        [group_of[i] for i in range(1, n + 1)],
    )
    test.assertEqual(w["canonical_barcodes"], sorted(w["canonical_barcodes"]))
    test.assertEqual(
        w["canonical_barcodes"],
        [g["canonical_barcode"] for g in w["barcodes"]],
    )


class DualValidationTests(unittest.TestCase):
    def test_count_bounds_dual(self):
        with self.assertRaises(ValidationError):
            assemble(["AAA"] * 5, barcode_count=2)
        with self.assertRaises(ValidationError):
            assemble(["AAA"] * 19, barcode_count=2)

    def test_bad_barcode_count(self):
        with self.assertRaises(ValidationError):
            assemble(["AAA"] * 6, barcode_count=3)
        with self.assertRaises(ValidationError):
            assemble(["AAA"] * 6, barcode_count=0)

    def test_default_is_single_strain_unchanged(self):
        r = assemble(reads_of("AATCGC", 3))
        self.assertEqual(r["status"], "unique")
        self.assertNotIn("barcode_count", r)


class DualUniqueTests(unittest.TestCase):
    C1, C2 = "AATCGC", "ACCTTG"  # 两环 k-mer 完全不相交 -> 划分唯一

    def setUp(self):
        self.reads = reads_of(self.C1, 3) + reads_of(self.C2, 3)

    def test_unique_two_barcodes(self):
        r = assemble(list(self.reads), barcode_count=2)
        self.assertEqual(r["status"], "unique")
        self.assertEqual(r["barcode_count"], 2)
        self.assertEqual(
            r["canonical_barcodes"],
            sorted([canonical_circle(self.C1), canonical_circle(self.C2)]),
        )
        check_dual_witness(self, r, self.reads, 2)

    def test_group_swap_equivalent(self):
        swapped = reads_of(self.C2, 3) + reads_of(self.C1, 3)
        r1 = assemble(list(self.reads), barcode_count=2)
        r2 = assemble(swapped, barcode_count=2)
        self.assertEqual(r2["status"], "unique")
        self.assertEqual(r1["canonical_barcodes"], r2["canonical_barcodes"])
        check_dual_witness(self, r2, swapped, 2)

    def test_shuffled_input_stable(self):
        rng = random.Random(7)
        sh = list(self.reads)
        rng.shuffle(sh)
        r = assemble(sh, barcode_count=2)
        self.assertEqual(r["status"], "unique")
        self.assertEqual(
            r["canonical_barcodes"],
            sorted([canonical_circle(self.C1), canonical_circle(self.C2)]),
        )
        check_dual_witness(self, r, sh, 2)
        self.assertEqual(
            assemble(sh, barcode_count=2)["assignment"], r["assignment"]
        )

    def test_reverse_complement_equivalent(self):
        rc = [reverse_complement(s) for s in self.reads]
        r = assemble(rc, barcode_count=2)
        self.assertEqual(r["status"], "unique")
        self.assertEqual(
            r["canonical_barcodes"],
            assemble(list(self.reads), barcode_count=2)["canonical_barcodes"],
        )

    def test_identical_barcode_sequences_unique_not_ambiguous(self):
        # 两条条码序列相同、但证据（各一份读数）足以拆成两份：
        # 组标签互换与重复实例分配都不得制造歧义。
        reads = reads_of("AATCGC", 3) * 2
        r = assemble(list(reads), barcode_count=2)
        self.assertEqual(r["status"], "unique", r.get("witnesses"))
        self.assertEqual(r["canonical_barcodes"], ["AATCGC", "AATCGC"])
        check_dual_witness(self, r, reads, 2)
        sizes = sorted(
            sum(1 for a in r["assignment"] if a["group"] == g) for g in (0, 1)
        )
        self.assertEqual(sizes, [6, 6])

    def test_identical_barcodes_with_internal_duplicate_reads(self):
        # 组标签互换与重复实例分配不制造歧义：本原环两份、每种读数各两条，
        # 跨组实例分配方式极多，但唯一可行的完整划分仍是 6+6 两份同码。
        reads = reads_of("AATCGC", 3) * 2
        rng = random.Random(11)
        sh = list(reads)
        rng.shuffle(sh)
        r = assemble(sh, barcode_count=2)
        self.assertEqual(r["status"], "unique", r.get("witnesses"))
        self.assertEqual(r["canonical_barcodes"], ["AATCGC", "AATCGC"])
        check_dual_witness(self, r, sh, 2)
        self.assertEqual(
            assemble(sh, barcode_count=2)["assignment"], r["assignment"]
        )

    def test_periodic_circle_duplicates_can_be_genuinely_ambiguous(self):
        # 周期 3 的环 AACAAC 两份：除 6+6（两条 AACAAC）外还存在 3+9
        # （AAC 与 AACAACAAC）的完整划分，这是真实歧义而非实例分配噪声。
        reads = reads_of("AACAAC", 3) * 2
        r = assemble(list(reads), barcode_count=2)
        self.assertEqual(r["status"], "ambiguous")
        pairs = {tuple(w["canonical_barcodes"]) for w in r["witnesses"]}
        self.assertIn(("AACAAC", "AACAAC"), pairs)
        self.assertTrue(any(len(set(p)) == 2 for p in pairs))
        expected = brute_dual_classes(reads)
        self.assertTrue(pairs <= expected)

    def test_minimum_3_plus_3_reads(self):
        reads = reads_of("AAC", 3) + reads_of("CGT", 3)
        r = assemble(list(reads), barcode_count=2)
        self.assertEqual(r["status"], "unique")
        check_dual_witness(self, r, reads, 2)


class DualAmbiguousTests(unittest.TestCase):
    # 两环共享读数 ATC：除天然 6+6 划分外，还能交错成 8+4 两个闭环，
    # 规范条码对不同 -> 歧义。
    READS = reads_of("AATCGC", 3) + reads_of("ACGATC", 3)

    def test_two_distinct_complete_witnesses(self):
        r = assemble(list(self.READS), barcode_count=2)
        self.assertEqual(r["status"], "ambiguous")
        self.assertEqual(r["barcode_count"], 2)
        self.assertEqual(len(r["witnesses"]), 2)
        pairs = [tuple(w["canonical_barcodes"]) for w in r["witnesses"]]
        self.assertEqual(len(set(pairs)), 2)
        for w in r["witnesses"]:
            check_dual_witness(self, w, self.READS, 2)

    def test_matches_brute_force(self):
        expected = brute_dual_classes(self.READS)
        self.assertGreaterEqual(len(expected), 2)
        got = {
            tuple(w["canonical_barcodes"])
            for w in assemble(list(self.READS), barcode_count=2)["witnesses"]
        }
        self.assertTrue(got <= expected)


class DualNoSolutionTests(unittest.TestCase):
    def test_single_circle_not_splittable(self):
        reads = reads_of("AATCGC", 3)
        r = assemble(list(reads), barcode_count=2)
        self.assertEqual(r["status"], "no_solution")
        codes = {x["code"] for x in r["reasons"]}
        self.assertEqual(codes, {"no_dual_partition"})
        reason = r["reasons"][0]
        self.assertEqual(
            reason["constraints"]["min_reads_per_barcode"], 3
        )
        self.assertGreaterEqual(reason["enumeration"]["partitions_checked"], 1)

    def test_global_degree_imbalance(self):
        reads = ["AAA", "AAC", "ACC", "CCC", "GGG", "TTT"]
        r = assemble(list(reads), barcode_count=2)
        self.assertEqual(r["status"], "no_solution")
        codes = {x["code"] for x in r["reasons"]}
        self.assertEqual(codes, {"degree_imbalance"})

    def test_no_solution_matches_brute_force(self):
        # 6 条只成一个环：任何 3+3 划分都不可能各自平衡连通
        reads = reads_of("AATCGC", 3)
        self.assertEqual(brute_dual_classes(reads), set())
        self.assertEqual(
            assemble(list(reads), barcode_count=2)["status"], "no_solution"
        )


class DualBruteForceCrossCheck(unittest.TestCase):
    """小实例上双株结论与独立暴力划分枚举一致（含唯一/歧义/无解）。"""

    def test_random_dual_batches(self):
        rng = random.Random(20261001)
        checked = 0
        for _ in range(24):
            n1 = rng.randint(3, 5)
            n2 = rng.randint(3, 5)
            alpha = "AC" if rng.random() < 0.5 else "ACGT"
            c1 = "".join(rng.choice(alpha) for _ in range(n1))
            c2 = "".join(rng.choice(alpha) for _ in range(n2))
            reads = reads_of(c1, 3) + reads_of(c2, 3)
            if len(reads) > 10:  # 暴力代价可控
                continue
            expected = brute_dual_classes(reads)
            r = assemble(list(reads), barcode_count=2)
            checked += 1
            if len(expected) == 1:
                self.assertEqual(r["status"], "unique", reads)
                self.assertEqual(
                    r["canonical_barcodes"], sorted(next(iter(expected))), reads
                )
                check_dual_witness(self, r, reads, 2)
            elif len(expected) > 1:
                self.assertEqual(r["status"], "ambiguous", reads)
                got = {tuple(w["canonical_barcodes"]) for w in r["witnesses"]}
                self.assertEqual(len(got), 2)
                self.assertTrue(got <= expected, reads)
            else:
                self.assertEqual(r["status"], "no_solution", reads)
        self.assertGreaterEqual(checked, 10)

    def test_random_single_circles_unsplittable_or_matched(self):
        rng = random.Random(4242)
        for _ in range(30):
            n = rng.randint(6, 8)
            reads = reads_of(
                "".join(rng.choice("ACGT") for _ in range(n)), 3
            )
            expected = brute_dual_classes(reads)
            r = assemble(list(reads), barcode_count=2)
            if not expected:
                self.assertEqual(r["status"], "no_solution", reads)
            elif len(expected) == 1:
                self.assertEqual(r["status"], "unique", reads)
                check_dual_witness(self, r, reads, 2)
            else:
                self.assertEqual(r["status"], "ambiguous", reads)
                got = {tuple(w["canonical_barcodes"]) for w in r["witnesses"]}
                self.assertTrue(got <= expected)


if __name__ == "__main__":
    unittest.main(verbosity=2)
