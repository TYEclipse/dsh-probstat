# CHANGELOG · dsh-probstat

> 版本口径：patch 修 bug/补测试｜minor 新增用户可见功能｜major 破坏性变更。安装：`dsh plugin --profile web add github:TYEclipse/dsh-probstat`

## [0.3.0] — 2026-10-21
### Minor · R75
- [自主进化] R75 two_sample_test：双样本检验（均值 z / Welch / 合并方差 + 双比例），Welch 分数自由度原值报告、合并方差与合并比例可核对
- `two_sample_test`：四种双样本检验——mean_z（均值差，σ₁/σ₂ 已知）/ mean_welch（σ 未知、方差不齐，Welch t + Welch–Satterthwaite 分数自由度）/ mean_pooled（σ 未知、等方差假设，合并方差 t，df = n₁+n₂−2）/ proportion（双比例差，合并比例正态近似）；可给 `delta0`（零假设差值，默认 0）、`tail`、`alpha`
- 一次给出两组估计、差值、标准误、统计量、自由度、三个 p 值（双侧/左/右）、对应尾部临界值与 α 判定；Welch df **不取整**（取整会静默改变 p 值，`note` 明写）；合并方差口径带 `pooledVariance`、比例口径带 `pooledProportion` 供读者核对假设；比例口径在最小期望单元格 < 5 时提示近似不可靠；标准误为 0（两组全为同一结果）时报错而不是给出无意义的统计量
- 测试 50 → 59；覆盖率 86.97% → 88.12%（基线棘轮上调）；oracle 扩 D 区（`test/oracle/anchors.py`），t 尾概率在本区**改用 Simpson 数值积分**独立实现（与插件的连分式不完全贝塔不同算法），并与公开 t 表值（2.228138852 / 2.042272456 / 1.795884819）以及两种算法互检；`--check` 自检覆盖 p 值域、双侧 = 2×较小尾、Welch df < 合并 df、合并方差 = 自由度加权均值等不变量

## [0.2.0] — 2026-10-02
### Minor · R56
- [自主进化] R56 新增 hypothesis_test（单样本 z/t/比例检验）与 sample_size（目标误差限所需样本量）
- `hypothesis_test`：三种单样本检验——mean_z（σ 已知）/ mean_t（σ 未知，df = n−1）/ proportion（对 p0 的正态近似）；一次给出统计量、标准误、三个 p 值（双侧/左/右）、对应尾部的临界值与 α 判定；比例检验在 n·p0 < 5 时用 `note` 显式提示近似不可靠，不静默给 p 值
- `sample_size`：目标误差限 E 所需样本量——mean_z 闭式 n = (z·σ/E)²；proportion 可用 wald 闭式或 wilson（二分求 Wilson 半宽 ≤ E 的最小 n）；`p` 默认 0.5（最保守），越界与不可达目标报明确错误
- 测试 37 → 50；覆盖率 83.70% → 86.97%（基线棘轮上调）；oracle 迁移到 `test/oracle/anchors.py`（测试文件头 `ORACLE:` 标记，A/B/C 三区覆盖含 legacy 全部数值断言 + `--check` 公开表值自检）

## [0.1.3] — 2026-09-11
### Patch · R31
- [自主进化] 接入版本与覆盖率门禁（工具链）

## [0.1.2] — 2026-09-11
### Patch · R31
- [自主进化] 接入版本与覆盖率门禁（工具链）

## [0.1.1] — 2026-09-11
### Patch · R31
- [自主进化] 接入版本与覆盖率门禁（工具链）

