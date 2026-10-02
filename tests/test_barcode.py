"""barcode 核心算法单元测试 + 独立暴力参照交叉验证。"""
import itertools
import os
import random
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "app"))

from barcode import (  # noqa: E402
    Edge,
    MIN_GROUP_READS,
    canonical_circle,
    reverse_complement,
    ValidationError,
    _build_degrees,
    _components,
    _enumerate_classes,
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


def brute_dual_classes(reads: list[str], k: int) -> set[tuple[str, str]]:
    """双株独立参照：枚举所有完整划分（读 0 固定归 A），组内枚举全部规范类，
    收集"两组规范条码无序对"的集合。"""
    n = len(reads)
    edges = [Edge(s[:-1], s[1:], s, i) for i, s in enumerate(reads)]
    pairs: set[tuple[str, str]] = set()
    for bits in range(1 << (n - 1)):
        mask = 1 | (bits << 1)
        am = [i for i in range(n) if (mask >> i) & 1]
        bm = [i for i in range(n) if not (mask >> i) & 1]
        if not (MIN_GROUP_READS <= len(am) <= n - MIN_GROUP_READS):
            continue
        ga = [edges[i] for i in am]
        gb = [edges[i] for i in bm]
        ia, oa, va = _build_degrees(ga)
        ib, ob, vb = _build_degrees(gb)
        if any(ia.get(v, 0) != oa.get(v, 0) for v in va):
            continue
        if any(ib.get(v, 0) != ob.get(v, 0) for v in vb):
            continue
        if len(_components(ga)) > 1 or len(_components(gb)) > 1:
            continue
        ca = {c for c, _, _ in _enumerate_classes(ga, va[0], k, limit=99)}
        cb = {c for c, _, _ in _enumerate_classes(gb, vb[0], k, limit=99)}
        for a in ca:
            for b in cb:
                pairs.add(tuple(sorted((a, b))))
    return pairs


def check_dual_witness(test: unittest.TestCase, result: dict,
                       reads: list[str], k: int) -> None:
    """双株见证一致性：每组闭环、每条读数恰好一次、归属与证据自洽。"""
    n = len(reads)
    if "status" in result:
        witnesses = [result] if result["status"] == "unique" else result["witnesses"]
    else:
        witnesses = [result]  # 单个 ambiguous 见证
    for wit in witnesses:
        groups = wit["groups"]
        test.assertEqual(len(groups), 2)
        all_orders: list[int] = []
        for rank, g in enumerate(groups, start=1):
            m = len(g["order"])
            test.assertGreaterEqual(m, MIN_GROUP_READS)
            all_orders.extend(g["order"])
            test.assertEqual(len(g["evidence"]), m)
            length = k + 1
            reps = (m + length) // m + 1
            doubled = g["barcode"] * reps
            for pos, ev in enumerate(g["evidence"]):
                test.assertEqual(ev["position"], pos)
                a = reads[ev["prev"] - 1]
                b = reads[ev["next"] - 1]
                test.assertEqual(ev["prev"], g["order"][pos])
                test.assertEqual(ev["next"], g["order"][(pos + 1) % m])
                test.assertEqual(a[1:], b[:-1])
                test.assertEqual(ev["overlap"], a[1:])
                test.assertEqual(ev["overlap_length"], k)
                test.assertEqual(doubled[pos:pos + length],
                                 reads[g["order"][pos] - 1])
            test.assertEqual(canonical_circle(g["barcode"]),
                             g["canonical_barcode"])
            for one in g["order"]:
                test.assertEqual(wit["assignment"][one - 1], rank)
        test.assertEqual(sorted(all_orders), list(range(1, n + 1)))
        test.assertEqual(wit["canonical_barcodes"],
                         [g["canonical_barcode"] for g in groups])
        test.assertEqual(sorted(wit["canonical_barcodes"]),
                         sorted(wit["canonical_barcodes"]))


class DualValidationTests(unittest.TestCase):
    def test_count_bounds(self):
        with self.assertRaises(ValidationError):
            assemble(["AAA"] * 5, 2)  # 少于 6
        with self.assertRaises(ValidationError):
            assemble(["AAA"] * 19, 2)  # 双株上限 18

    def test_bad_barcode_count(self):
        with self.assertRaises(ValidationError):
            assemble(["AAA"] * 6, 0)
        with self.assertRaises(ValidationError):
            assemble(["AAA"] * 6, 3)
        with self.assertRaises(ValidationError):
            assemble(["AAA"] * 6, "2")  # type: ignore[arg-type]

    def test_explicit_one_keeps_single_behavior(self):
        reads = ["AAT", "ATC", "TCG", "CGC", "GCA", "CAA"]
        r = assemble(list(reads), 1)
        self.assertEqual(r["status"], "unique")
        self.assertEqual(r["canonical_barcode"], "AATCGC")


class DualUniqueTests(unittest.TestCase):
    # 两条 4 环，k-mer 字母表互不相交，划分唯一。
    R1 = ["AAT", "ATC", "TCA", "CAA"]   # 环 AATC
    R2 = ["GTT", "TTG", "TGG", "GGT"]   # 环 GT TG.. 规范代表 AACC

    def test_unique_partition(self):
        reads = self.R1 + self.R2
        r = assemble(list(reads), 2)
        self.assertEqual(r["status"], "unique")
        self.assertEqual(r["canonical_barcodes"], ["AACC", "AATC"])
        self.assertEqual(r["read_count"], 8)
        self.assertEqual(r["barcode_count"], 2)
        check_dual_witness(self, r, reads, 2)

    def test_shuffle_and_group_swap_equivalent(self):
        reads = self.R1 + self.R2
        rng = random.Random(20261002)
        pairs = set()
        for _ in range(30):
            sh = reads[:]
            rng.shuffle(sh)
            r = assemble(sh, 2)
            self.assertEqual(r["status"], "unique", sh)
            self.assertEqual(r["canonical_barcodes"], ["AACC", "AATC"])
            pairs.add(tuple(r["canonical_barcodes"]))
            check_dual_witness(self, r, sh, 2)
        self.assertEqual(len(pairs), 1)  # 组交换不制造差异

    def test_assignment_covers_every_read_once(self):
        r = assemble(self.R1 + self.R2, 2)
        self.assertEqual(sorted(r["assignment"]), [1, 1, 1, 1, 2, 2, 2, 2])


class DualIdenticalBarcodeTests(unittest.TestCase):
    """两条条码序列相同、但证据足以分成两份时，必须给唯一双株结论。"""

    READS = ["AAC", "ACA", "CAA"] * 2  # 环 AAC 的两个实例

    def test_unique_even_though_barcodes_identical(self):
        r = assemble(list(self.READS), 2)
        self.assertEqual(r["status"], "unique")
        self.assertEqual(r["canonical_barcodes"], ["AAC", "AAC"])
        self.assertEqual(r["assignment"], [1, 1, 1, 2, 2, 2])
        check_dual_witness(self, r, self.READS, 2)

    def test_shuffled_identical_still_unique(self):
        rng = random.Random(11)
        for _ in range(20):
            sh = self.READS[:]
            rng.shuffle(sh)
            r = assemble(sh, 2)
            self.assertEqual(r["status"], "unique", sh)
            self.assertEqual(r["canonical_barcodes"], ["AAC", "AAC"])
            check_dual_witness(self, r, sh, 2)


class DualAmbiguousTests(unittest.TestCase):
    # 3xAAA 自环 + 3xCCC 自环 + A 环(AATC 结构边: AAC,ACC,CCA,CAA)。
    # 跨环共享顶点 AA/CC 时存在两份不同完整划分：
    #   {'AAA','AACCCCC'} 与 {'AAAAACC','CCC'}
    READS = ["AAA", "AAA", "AAA", "CCC", "CCC", "CCC",
             "AAC", "ACC", "CCA", "CAA"]

    def test_two_distinct_complete_partitions(self):
        r = assemble(list(self.READS), 2)
        self.assertEqual(r["status"], "ambiguous")
        self.assertEqual(len(r["witnesses"]), 2)
        pairs = {tuple(w["canonical_barcodes"]) for w in r["witnesses"]}
        self.assertEqual(len(pairs), 2)
        self.assertEqual(
            pairs,
            {("AAA", "AACCCCC"), ("AAAAACC", "CCC")},
        )
        for w in r["witnesses"]:
            check_dual_witness(self, w, self.READS, 2)

    def test_matches_brute_force(self):
        expected = brute_dual_classes(self.READS, 2)
        self.assertEqual(len(expected), 2)
        r = assemble(list(self.READS), 2)
        got = {tuple(w["canonical_barcodes"]) for w in r["witnesses"]}
        self.assertTrue(got <= expected)


class DualNoSolutionTests(unittest.TestCase):
    def test_globally_imbalanced_unsplittable(self):
        reads = ["AAA", "AAA", "AAA", "AAC", "CCC", "GGG"]
        r = assemble(list(reads), 2)
        self.assertEqual(r["status"], "no_solution")
        codes = {x["code"] for x in r["reasons"]}
        self.assertEqual(codes, {"no_dual_partition"})
        reason = r["reasons"][0]
        self.assertIn("global_degree_imbalance", reason)
        self.assertGreaterEqual(reason["partitions_examined"], 1)

    def test_two_self_loops_plus_isolated_read_splits_fine_else_not(self):
        # 6 条全是互不相同 k-mer 的自环 -> 每组 3 条可拆，应有解
        reads = ["AAA", "AAA", "AAA", "CCC", "CCC", "CCC"]
        r = assemble(list(reads), 2)
        self.assertEqual(r["status"], "unique")
        self.assertEqual(r["canonical_barcodes"], ["AAA", "CCC"])
        # 换成 4 个不同自环，无法凑出每组 >=3 的两个平衡组
        reads2 = ["AAA", "AAA", "CCC", "CCC", "GGG", "GGG"]
        r2 = assemble(list(reads2), 2)
        self.assertEqual(r2["status"], "no_solution")
        self.assertEqual(r2["reasons"][0]["code"], "no_dual_partition")


class DualBruteForceCrossCheck(unittest.TestCase):
    """随机小实例：双株裁决必须与独立全划分暴力参照一致。"""

    def test_random_instances(self):
        rng = random.Random(20261002)
        for _ in range(80):
            n = rng.randint(6, 10)
            ncirc = rng.choice([1, 2, 2, 3])
            base = [n // ncirc] * ncirc
            for i in range(n - sum(base)):
                base[i] += 1
            reads: list[str] = []
            for cnt in base:
                length = rng.randint(3, 5)
                circle = "".join(rng.choice("ACGT") for _ in range(length))
                pool = reads_of(circle, 3)
                if cnt > len(pool):
                    chosen = [rng.choice(pool) for _ in range(cnt)]
                else:
                    chosen = rng.sample(pool, cnt)
                reads.extend(chosen)
            rng.shuffle(reads)
            expected = brute_dual_classes(reads, 2)
            r = assemble(list(reads), 2)
            if not expected:
                self.assertEqual(r["status"], "no_solution", reads)
            elif len(expected) == 1:
                self.assertEqual(r["status"], "unique", reads)
                self.assertEqual(tuple(r["canonical_barcodes"]),
                                 next(iter(expected)), reads)
                check_dual_witness(self, r, reads, 2)
            else:
                self.assertEqual(r["status"], "ambiguous", reads)
                got = {tuple(w["canonical_barcodes"]) for w in r["witnesses"]}
                self.assertEqual(len(got), 2)
                self.assertTrue(got <= expected)
                for w in r["witnesses"]:
                    check_dual_witness(self, w, reads, 2)


if __name__ == "__main__":
    unittest.main(verbosity=2)
