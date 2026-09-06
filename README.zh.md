# dsh-probstat

DeepSeek Harness（dsh）概率与统计推断数学工具箱。零运行时依赖，四个纯确定性工具，专治 agent 最容易算错的几类题：背错的 z 表值、分布 CDF/分位数算错、置信区间临界值取错、复合事件概率恒等式用错。

## 安装

```bash
dsh plugin --profile web add github:TYEclipse/dsh-probstat
```

## 工具

- `dist_calc`：六种分布（正态/二项/泊松/指数/均匀/几何）的 pdf / cdf / 生存函数 / 分位数 / 统计量（均值·方差·标准差）。
- `z_score`：标准正态表换算（z ↔ 概率/百分位、双侧 p 值、区间概率）。
- `confidence_interval`：均值置信区间（σ 已知走 z、未知走 t，df = n−1）+ 比例 Wilson 区间（小样本与极端比例下依然可靠）。
- `event_probability`：事件概率恒等式——并（独立/互斥）/交（独立）/条件/贝叶斯/补集/至少一次成功。

## 约定

- 几何分布按「首次成功前的失败次数」计数（k ≥ 0，P(X = k) = p·(1−p)^k）。
- 离散分布的分位数 = 最小的 k 使 P(X ≤ k) ≥ q。
- 数值输入非法时返回 `valid: false` + `error` 说明，不抛异常。

## 精度

- 正态 CDF：A&S 7.1.26 近似，绝对误差 ≤ 1.5e-7（远超公开 4 位 z 表的精度）。
- 正态分位数：Acklam 逆正态，相对误差 < 1.15e-9。
- 二项/泊松：递推精确求和；t 分布：连分式不完全贝塔 + Lanczos log-Gamma，已与公开 t 表值（2.2281 / 2.0452 / 1.8125 / 1.984）交叉验证。

## 开发

```bash
pnpm install && pnpm build && pnpm test && pnpm lint
```

测试锚点全部由独立 Python oracle（test/anchors-probstat.py，math.erf / math.comb 精确原语 + 二分）生成并与公开表值交叉验证。

## 许可证

MIT
