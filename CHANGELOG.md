# CHANGELOG · dsh-probstat

> 版本口径：patch 修 bug/补测试｜minor 新增用户可见功能｜major 破坏性变更。安装：`dsh plugin --profile web add github:TYEclipse/dsh-probstat`

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

