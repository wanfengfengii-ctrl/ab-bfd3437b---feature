"""一次性 verify 容器：

  1. 等待 API /health 就绪（内部轮询，满足"verify 服务等待 API 就绪"）；
  2. 构建自检（compileall 字节码编译）；
  3. 运行单元测试；
  4. 通过 HTTP 对 /assemble 做唯一、歧义、无解三类冒烟；
  5. 汇总结果后自行退出，全部通过退出码 0，否则非零。
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import time
import urllib.error
import urllib.request

API_HOST = os.environ.get("API_HOST", "api")
API_PORT = os.environ.get("API_PORT", "8080")
BASE = f"http://{API_HOST}:{API_PORT}"
ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(ROOT, "app"))

failures: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> bool:
    mark = "PASS" if ok else "FAIL"
    print(f"[{mark}] {name}" + (f" -- {detail}" if detail else ""), flush=True)
    if not ok:
        failures.append(name)
    return ok


def wait_ready(timeout: float = 60.0) -> bool:
    deadline = time.time() + timeout
    last = ""
    while time.time() < deadline:
        try:
            with urllib.request.urlopen(f"{BASE}/health", timeout=3) as resp:
                if resp.status == 200:
                    body = json.loads(resp.read().decode())
                    return body.get("status") == "ok"
        except (urllib.error.URLError, OSError, json.JSONDecodeError) as exc:
            last = str(exc)
            time.sleep(0.5)
    print(f"等待 API 就绪超时: {last}")
    return False


def post_assemble(reads: list[str], barcode_count: int | None = None) -> tuple[int, dict]:
    payload: dict = {"sequences": reads}
    if barcode_count is not None:
        payload["barcode_count"] = barcode_count
    req = urllib.request.Request(
        f"{BASE}/assemble",
        data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    try:
        with urllib.request.urlopen(req, timeout=10) as resp:
            return resp.status, json.loads(resp.read().decode())
    except urllib.error.HTTPError as exc:
        return exc.code, json.loads(exc.read().decode())


def smoke_unique() -> None:
    reads = ["AAT", "ATC", "TCG", "CGC", "GCA", "CAA"]
    code, body = post_assemble(reads)
    ok = (
        code == 200
        and body.get("status") == "unique"
        and body.get("canonical_barcode") == "AATCGC"
        and sorted(body.get("order", [])) == list(range(1, 7))
        and len(body.get("evidence", [])) == 6
        and all(e["overlap_length"] == 2 for e in body["evidence"])
    )
    check("冒烟-唯一拼接", bool(ok),
          f"status={body.get('status')} barcode={body.get('canonical_barcode')}")

    # 反向互补整组提交必须得到同一规范条码（等价归一化）
    from barcode import reverse_complement
    code2, body2 = post_assemble([reverse_complement(s) for s in reads])
    check("冒烟-反向互补等价",
          code2 == 200 and body2.get("canonical_barcode") == "AATCGC",
          f"barcode={body2.get('canonical_barcode')}")


def smoke_ambiguous() -> None:
    reads = ["AAC", "ACC", "AAG", "AGC", "GCC",
             "CCA", "CAA", "CCG", "CGA", "GAA"]
    code, body = post_assemble(reads)
    witnesses = body.get("witnesses", [])
    barcodes = {w.get("canonical_barcode") for w in witnesses}
    ok = (
        code == 200
        and body.get("status") == "ambiguous"
        and len(witnesses) == 2
        and len(barcodes) == 2
        and all(sorted(w.get("order", [])) == list(range(1, 11)) for w in witnesses)
        and all(len(w.get("evidence", [])) == 10 for w in witnesses)
    )
    check("冒烟-歧义双见证", bool(ok),
          f"witnesses={sorted(barcodes)}")


def smoke_no_solution() -> None:
    # 两个原因同时存在：度数失衡 + 多分量
    reads = ["AAA", "AAA", "AAC", "CCC", "GGG", "TTT"]
    code, body = post_assemble(reads)
    reasons = {r.get("code") for r in body.get("reasons", [])}
    ok = (
        code == 200
        and body.get("status") == "no_solution"
        and {"degree_imbalance", "fragmented_graph"} <= reasons
    )
    check("冒烟-无解(度数失衡+非零片段不连通)", bool(ok),
          f"reasons={sorted(reasons)}")

    # 仅多分量但各自平衡
    reads2 = ["AAA", "AAA", "CCC", "CCC", "GGG", "GGG"]
    code2, body2 = post_assemble(reads2)
    reasons2 = {r.get("code") for r in body2.get("reasons", [])}
    check("冒烟-无解(多分量各自平衡)",
          code2 == 200 and body2.get("status") == "no_solution"
          and reasons2 == {"fragmented_graph"},
          f"reasons={sorted(reasons2)}")


def smoke_invalid() -> None:
    code, body = post_assemble(["AAA"] * 5)  # 少于 6 条
    check("冒烟-非法输入 400", code == 400 and body.get("status") == "invalid_input",
          f"code={code}")


def _reads_of(circle: str, length: int = 3) -> list[str]:
    doubled = circle * (length // len(circle) + 2)
    return [doubled[i:i + length] for i in range(len(circle))]


def _check_dual_witness(body: dict, n: int, expected: list[str]) -> bool:
    """核对完整双组见证：两组、每条读数恰好一次、组内闭环证据自洽。"""
    if body.get("canonical_barcodes") != expected:
        return False
    groups = body.get("barcodes", [])
    if len(groups) != 2 or any(len(g.get("order", [])) < 3 for g in groups):
        return False
    all_orders = sorted(i for g in groups for i in g.get("order", []))
    if all_orders != list(range(1, n + 1)):
        return False
    assignment = body.get("assignment", [])
    if [a.get("read_index") for a in assignment] != list(range(1, n + 1)):
        return False
    group_of = {}
    for gi, g in enumerate(groups):
        order = g["order"]
        if len(g.get("evidence", [])) != len(order):
            return False
        for pos, ev in enumerate(g["evidence"]):
            if (ev["prev"], ev["next"]) != (order[pos], order[(pos + 1) % len(order)]):
                return False
            if ev["overlap_length"] != 2:
                return False
        for idx in order:
            group_of[idx] = gi
    if [a.get("group") for a in assignment] != [group_of[i] for i in range(1, n + 1)]:
        return False
    return True


def smoke_dual_unique() -> None:
    # 两环 k-mer 互不相交 -> 唯一完整双组划分
    c1, c2 = "AATCGC", "ACCTTG"
    from barcode import canonical_circle
    expected = sorted([canonical_circle(c1), canonical_circle(c2)])
    reads = _reads_of(c1) + _reads_of(c2)
    code, body = post_assemble(reads, barcode_count=2)
    ok = code == 200 and body.get("status") == "unique" and body.get("barcode_count") == 2
    ok = ok and _check_dual_witness(body, len(reads), expected)
    check("冒烟-双株唯一", bool(ok),
          f"status={body.get('status')} barcodes={body.get('canonical_barcodes')}")

    # 组交换：输入顺序互换两段，结论与归属稳定
    swapped = _reads_of(c2) + _reads_of(c1)
    code2, body2 = post_assemble(swapped, barcode_count=2)
    check("冒烟-双株组交换等价",
          code2 == 200 and body2.get("status") == "unique"
          and body2.get("canonical_barcodes") == expected,
          f"barcodes={body2.get('canonical_barcodes')}")

    # 打乱 + 反向互补同样唯一
    import random
    rng = random.Random(20261001)
    shuffled = list(reads)
    rng.shuffle(shuffled)
    code3, body3 = post_assemble(shuffled, barcode_count=2)
    from barcode import reverse_complement
    rc = [reverse_complement(s) for s in shuffled]
    code4, body4 = post_assemble(rc, barcode_count=2)
    check("冒烟-双株打乱/反向互补等价",
          code3 == 200 and body3.get("status") == "unique"
          and code4 == 200 and body4.get("status") == "unique"
          and body3.get("canonical_barcodes") == expected
          and body4.get("canonical_barcodes") == expected,
          f"{body3.get('canonical_barcodes')} {body4.get('canonical_barcodes')}")

    # 两条条码序列相同、证据足以拆两份：必须唯一，不得因组标签/实例分配误报
    same = _reads_of("AATCGC") * 2
    code5, body5 = post_assemble(same, barcode_count=2)
    ok5 = (
        code5 == 200
        and body5.get("status") == "unique"
        and body5.get("canonical_barcodes") == ["AATCGC", "AATCGC"]
    )
    check("冒烟-双株同码两实例仍唯一", bool(ok5),
          f"status={body5.get('status')} barcodes={body5.get('canonical_barcodes')}")


def smoke_dual_ambiguous() -> None:
    # 共享读数 ATC 的两环：6+6 与 8+4 两种完整划分给出不同规范条码对
    reads = _reads_of("AATCGC") + _reads_of("ACGATC")
    code, body = post_assemble(reads, barcode_count=2)
    witnesses = body.get("witnesses", [])
    pairs = [tuple(w.get("canonical_barcodes")) for w in witnesses]
    ok = (
        code == 200
        and body.get("status") == "ambiguous"
        and body.get("barcode_count") == 2
        and len(witnesses) == 2
        and len(set(pairs)) == 2
        and all(
            sorted(i for g in w.get("barcodes", []) for i in g.get("order", []))
            == list(range(1, len(reads) + 1))
            for w in witnesses
        )
    )
    check("冒烟-双株歧义双见证", bool(ok), f"witnesses={sorted(pairs)}")


def smoke_dual_no_solution() -> None:
    # 单个闭环 6 条：任何 3+3 划分都成不了两个欧拉子图
    reads = _reads_of("AATCGC")
    code, body = post_assemble(reads, barcode_count=2)
    reasons = {r.get("code") for r in body.get("reasons", [])}
    check("冒烟-双株不可拆分无解",
          code == 200 and body.get("status") == "no_solution"
          and reasons == {"no_dual_partition"},
          f"reasons={sorted(reasons)}")

    # 非法条数（双株 6~18）与非法 barcode_count -> 400
    code2, body2 = post_assemble(["AAA"] * 19, barcode_count=2)
    code3, body3 = post_assemble(["AAA"] * 6, barcode_count=3)
    check("冒烟-双株非法输入 400",
          code2 == 400 and body2.get("status") == "invalid_input"
          and code3 == 400 and body3.get("status") == "invalid_input",
          f"codes={code2},{code3}")


def smoke_backward_compat() -> None:
    # 不带 barcode_count 的请求载荷与裁决与原行为一致
    reads = _reads_of("AATCGC")
    code1, body1 = post_assemble(reads)
    code2, body2 = post_assemble(reads, barcode_count=1)
    check("冒烟-单株向后兼容(缺省==1)",
          code1 == 200 and code2 == 200
          and body1 == body2
          and body1.get("status") == "unique"
          and body1.get("canonical_barcode") == "AATCGC"
          and "barcode_count" not in body1,
          f"status={body1.get('status')}")


def main() -> int:
    print(f"verify: target={BASE}", flush=True)
    if not check("等待 API 就绪", wait_ready()):
        return 1

    build = subprocess.run(
        [sys.executable, "-m", "compileall", "-q", "app", "tests", "verify.py"],
        cwd=ROOT,
    )
    check("构建自检 compileall", build.returncode == 0)

    tests = subprocess.run(
        [sys.executable, "-m", "unittest", "discover", "-s", "tests", "-v"],
        cwd=ROOT,
    )
    check("单元测试套件", tests.returncode == 0)

    smoke_unique()
    smoke_ambiguous()
    smoke_no_solution()
    smoke_invalid()
    smoke_dual_unique()
    smoke_dual_ambiguous()
    smoke_dual_no_solution()
    smoke_backward_compat()

    if failures:
        print(f"\nverify 失败 {len(failures)} 项: {failures}", flush=True)
        return 1
    print("\nverify 全部通过", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
