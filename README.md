# 环状 DNA 条码拼接服务（barcode-assembler）

一批来自同一环状条码的等长 DNA 短序列，判定能否拼成**可信的唯一条码**，
避免把重复片段在不同排列下的差异误报为不同病原株。

## 方法

每条长度 `L = k+1` 的读数视为德布鲁因（de Bruijn）多重图中的一条有向边
`s[:-1] → s[1:]`；**每次出现都是一条独立平行边（一份证据）**。
"每条读数恰好使用一次、相邻重叠长度为 `L-1` 的闭环"等价于该多重图的
欧拉回路。

- **唯一**：所有闭环拼法同属一个规范等价类；
- **歧义**：存在两个以上规范等价类，返回两条不同规范见证；
- **无解**：给出可复核原因——`degree_imbalance`（入度≠出度的 k-mer 清单）
  和/或 `fragmented_graph`（非零度顶点的弱连通分量）。

规范等价归一化：

1. 环的循环移位；
2. 整条环的反向互补；
3. 相同读数（平行边）互换不产生新类别（枚举时同标签只取当前最小输入序号，
   因而**使用次序按输入序号稳定分配**且跨次运行确定）。

枚举从固定起点出发做带桥剪枝的欧拉回路 DFS，按规范代表去重，找到 2 个
不同等价类即停。约束（6~24 条、长度 3~8）下毫秒级完成。

## 双株模式（`barcode_count=2`）

当怀疑某批次混入两株病原时，请求中加 `"barcode_count": 2`，表示**全部读数
来自恰好两条环状条码**（接受 6~18 条读数，每条仍须恰好使用一次）。

求解器**联合枚举完整划分**而不是"先拼一条、再处理余料"：

- 把每条读数实例分入 A/B 两组（读数 1 固定归 A 以消除"两条条码不分先后"的
  组交换对称），每组 3..n−3 条；
- 以 Gray 码走查全部 2^(n−1) 种分配并增量维护度数，只接受两组边子图
  **同时**度数平衡且弱连通的划分，再在组内各自按上述欧拉回路规则枚举闭环；
- 一个完整划分的**规范等价类**＝两组规范条码组成的无序对。循环移位、整条
  反向互补、组交换、相同重复读数（平行边）互换都不产生新类。

裁决：

- **唯一**：只有一个规范条码对。返回 `canonical_barcodes`（按规范条码排序的
  两组）、`groups`（两组各自的条码、稳定输入序号 `order`、组内重叠证据
  `evidence`）以及逐条读数归属 `assignment`（1/2）；
- **歧义**：返回 `witnesses`，含两份规范条码对不同的完整双组见证；
- **无解**：`no_dual_partition`，给出枚举划分数、每组最少读数等可复核信息；
  全体 k-mer 度数失衡时附 `global_degree_imbalance`（两组平衡之和必为 0，
  故全体失衡则不可能拆开）。

两条条码**序列相同但证据足以分成两份**（如同一 3 环各一份重复读数）时，
结论仍是唯一双株——组标签与重复实例的分配差异不会误报歧义。

缺省或显式传 `barcode_count=1` 时，请求、响应与裁决与单株模式完全一致。

## API

零第三方依赖（Python 标准库）。

- `GET /health` → `{"status":"ok",...}`
- `POST /assemble`，请求体：

```json
{ "sequences": ["AAT", "ATC", "TCG", "CGC", "GCA", "CAA"] }
```

双株批次可加选填字段 `barcode_count`（缺省/`1` 为单株，`2` 为双株）：

```json
{ "barcode_count": 2,
  "sequences": ["AAT", "ATC", "TCA", "CAA", "GTT", "TTG", "TGG", "GGT"] }
```

响应 `status` 为 `unique` / `ambiguous` / `no_solution`（非法输入返回 400）。

`unique` 响应包含：

- `canonical_barcode`：规范条码（环旋转+反向互补下的字典序最小代表）；
- `barcode`：与 `order` 同方向的环；
- `order`：使用次序，元素为 **1 起始的输入序号**（每个序号恰好出现一次）；
- `evidence`：相邻重叠证据（`prev`/`next` 序号、`overlap`、`overlap_length = L-1`、
  `appended_base`），最后一条与首条首尾相接。

`ambiguous` 响应在 `witnesses` 中给出两条不同规范见证（结构同上）。

## 运行

```bash
# 启动 API（端口可用 API_PORT 覆盖）
API_PORT=9090 docker compose up --build api

# 一次性 verify：等 API 健康检查通过后，执行构建自检 + 单元测试 +
# 单株唯一/歧义/无解冒烟 + 双株唯一/组交换/同序列双株/歧义/不可拆分冒烟，
# 随后自行退出，退出码 0 表示全部通过
docker compose up --build --exit-code-from verify verify
echo "verify exit code: $?"
```

`docker compose up --abort-on-container-exit verify` 或查看退出码：

```bash
docker compose run verify; echo "exit=$?"
```

本地不使用容器时：

```bash
python3 app/main.py            # 启 API（API_PORT 环境变量可配）
python3 -m unittest discover -s tests -v
python3 verify.py              # 需先启动 API，API_HOST/API_PORT 可配
```

## 目录

```
app/barcode.py    核心算法：校验、建图、欧拉回路枚举、规范归一化
app/main.py       HTTP 服务（/health、/assemble）
tests/            单元测试（含全排列暴力参照交叉验证）
verify.py         一次性 verify 服务（等待就绪→构建→测试→三类冒烟→退出码）
Dockerfile        单一镜像，默认启动 API
docker-compose.yml api + 一次性 verify
```
