"""一次性 verify 容器：

  1. 等待 API /health 就绪（内部轮询，满足"verify 服务等待 API 就绪"）；
  2. 构建自检（compileall 字节码编译）；
  3. 运行单元测试；
  4. 通过 HTTP 对 /assemble 做唯一、歧义、无解三类冒烟；
  5. 双株模式（barcode_count=2）冒烟：唯一划分、组交换等价、
     同序列双株不误报歧义、歧义双见证、不可拆分、默认兼容；
  6. 汇总结果后自行退出，全部通过退出码 0，否则非零。
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


def post_assemble(reads: list[str], barcode_count=None) -> tuple[int, dict]:
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


def _check_dual_witness(wit: dict, n: int) -> bool:
    """见证自洽：两组、每组 >=3 条、全部读数恰好出现一次、归属一致。"""
    groups = wit.get("groups", [])
    if len(groups) != 2:
        return False
    used: list[int] = []
    assignment = wit.get("assignment", [])
    if len(assignment) != n or set(assignment) != {1, 2}:
        return False
    for rank, g in enumerate(groups, start=1):
        order = g.get("order", [])
        if len(order) < 3:
            return False
        used.extend(order)
        if len(g.get("evidence", [])) != len(order):
            return False
        for one in order:
            if assignment[one - 1] != rank:
                return False
    return sorted(used) == list(range(1, n + 1))


def smoke_dual_unique() -> None:
    # 两条 4 环，k-mer 字母表不相交，划分唯一（GT.. 环规范代表为 AACC）
    r1 = ["AAT", "ATC", "TCA", "CAA"]
    r2 = ["GTT", "TTG", "TGG", "GGT"]
    reads = r1 + r2
    code, body = post_assemble(reads, 2)
    ok = (
        code == 200
        and body.get("status") == "unique"
        and body.get("barcode_count") == 2
        and body.get("canonical_barcodes") == ["AACC", "AATC"]
        and sorted(body.get("assignment", [])) == [1, 1, 1, 1, 2, 2, 2, 2]
        and _check_dual_witness(body, 8)
    )
    check("冒烟-双株唯一划分", bool(ok),
          f"status={body.get('status')} barcodes={body.get('canonical_barcodes')}")

    # 组交换等价：打乱输入序号，规范条码对必须不变、结论仍唯一
    shuffled = ["TGG", "ATC", "GGT", "CAA", "GTT", "AAT", "TTG", "TCA"]
    code2, body2 = post_assemble(shuffled, 2)
    check("冒烟-双株组交换/乱序等价",
          code2 == 200 and body2.get("status") == "unique"
          and body2.get("canonical_barcodes") == ["AACC", "AATC"]
          and _check_dual_witness(body2, 8),
          f"status={body2.get('status')} barcodes={body2.get('canonical_barcodes')}")

    # 两条条码序列相同但证据可分成两份：必须是唯一双株，不得误报歧义
    ident = ["AAC", "ACA", "CAA", "AAC", "ACA", "CAA"]
    code3, body3 = post_assemble(ident, 2)
    check("冒烟-同序列双株仍唯一(重复实例不误报歧义)",
          code3 == 200 and body3.get("status") == "unique"
          and body3.get("canonical_barcodes") == ["AAC", "AAC"]
          and _check_dual_witness(body3, 6),
          f"status={body3.get('status')} barcodes={body3.get('canonical_barcodes')}")


def smoke_dual_ambiguous() -> None:
    reads = ["AAA", "AAA", "AAA", "CCC", "CCC", "CCC",
             "AAC", "ACC", "CCA", "CAA"]
    code, body = post_assemble(reads, 2)
    witnesses = body.get("witnesses", [])
    pairs = {tuple(w.get("canonical_barcodes")) for w in witnesses}
    ok = (
        code == 200
        and body.get("status") == "ambiguous"
        and len(witnesses) == 2
        and pairs == {("AAA", "AACCCCC"), ("AAAAACC", "CCC")}
        and all(_check_dual_witness(w, 10) for w in witnesses)
    )
    check("冒烟-双株歧义两份完整见证", bool(ok), f"pairs={sorted(pairs)}")


def smoke_dual_no_solution() -> None:
    # 全体度数失衡：任何划分都不可能两组同时平衡
    reads = ["AAA", "AAA", "AAA", "AAC", "CCC", "GGG"]
    code, body = post_assemble(reads, 2)
    reasons = body.get("reasons", [])
    ok = (
        code == 200
        and body.get("status") == "no_solution"
        and len(reasons) == 1
        and reasons[0].get("code") == "no_dual_partition"
        and "global_degree_imbalance" in reasons[0]
        and reasons[0].get("partitions_examined", 0) >= 1
    )
    check("冒烟-双株不可拆分(可复核原因)", bool(ok),
          f"reasons={[r.get('code') for r in reasons]}")

    # 各自平衡但无法凑出每组 >=3 的两个连通闭环（4 个不同自环各 2 条）
    reads2 = ["AAA", "AAA", "CCC", "CCC", "GGG", "GGG"]
    code2, body2 = post_assemble(reads2, 2)
    check("冒烟-双株不可拆分(组规模不满足)",
          code2 == 200 and body2.get("status") == "no_solution"
          and body2["reasons"][0]["code"] == "no_dual_partition",
          f"status={body2.get('status')}")


def smoke_dual_compat_and_validation() -> None:
    # 显式 barcode_count=1 与缺省行为一致
    reads = ["AAT", "ATC", "TCG", "CGC", "GCA", "CAA"]
    code, body = post_assemble(reads, 1)
    check("冒烟-显式 barcode_count=1 兼容",
          code == 200 and body.get("status") == "unique"
          and body.get("canonical_barcode") == "AATCGC",
          f"status={body.get('status')}")

    # 双株模式 19 条超出 6~18 上限
    code2, body2 = post_assemble(["AAA"] * 19, 2)
    check("冒烟-双株 19 条非法 400",
          code2 == 400 and body2.get("status") == "invalid_input",
          f"code={code2}")

    # barcode_count 取值非法
    code3, body3 = post_assemble(reads, 3)
    check("冒烟-barcode_count=3 非法 400",
          code3 == 400 and body3.get("status") == "invalid_input",
          f"code={code3}")


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
    smoke_dual_compat_and_validation()

    if failures:
        print(f"\nverify 失败 {len(failures)} 项: {failures}", flush=True)
        return 1
    print("\nverify 全部通过", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
