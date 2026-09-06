# dsh-probstat

Deterministic probability & statistical-inference math for [DeepSeek Harness](https://github.com/deepseek-ai/deepseek-harness) (`dsh`). Four zero-runtime-dependency tools that replace the math agents get wrong most often: hallucinated z-table values, distribution cdf/quantile errors, wrong critical values, and sloppy compound-event identities.

中文简介：DeepSeek Harness 概率与统计推断数学工具箱——分布计算（正态/二项/泊松/指数/均匀/几何的 pdf/cdf/生存函数/分位数/统计量）、标准正态 z 表换算、置信区间（均值 z/t + 比例 Wilson）、事件概率恒等式（并/交/条件/贝叶斯）。零运行时依赖，纯确定性计算。

## Install

```bash
dsh plugin --profile web add github:TYEclipse/dsh-probstat
```

## Tools

| Tool | What it does |
|---|---|
| `dist_calc` | pdf / cdf / survival / quantile / stats for 6 distributions |
| `z_score` | standard-normal table math (z ↔ probability) |
| `confidence_interval` | mean intervals (z or t) and Wilson proportion intervals |
| `event_probability` | union / intersection / conditional / Bayes / complement / at-least-one |

### dist_calc

Distributions: `normal` (mean, sd — default standard normal), `binomial` (n, p), `poisson` (lambda), `exponential` (rate lambda), `uniform` (a, b), `geometric` (p). Operations: `pdf`, `cdf` (P(X ≤ x)), `survival` (P(X > x)), `quantile` (inverse cdf at probability q), `stats` (mean / variance / sd).

Examples:

- `distribution: "normal", x: 1.96` → cdf `0.9750022` (P(Z ≤ 1.96))
- `distribution: "binomial", operation: "quantile", n: 10, p: 0.5, q: 0.95` → `8`
- `distribution: "poisson", operation: "survival", lambda: 3, x: 5` → `0.08391794`
- `distribution: "geometric", operation: "stats", p: 0.5` → mean `1`, variance `2`

Conventions: the geometric distribution counts **failures before the first success** (support k ≥ 0, P(X = k) = p·(1−p)^k). Quantiles of discrete distributions return the smallest k with P(X ≤ k) ≥ q.

### z_score

Modes: `to_prob` (z → probability + percentile), `from_prob` (probability → z), `two_sided` (z → two-tailed p-value), `between` (z1, z2 → probability between; order is sorted internally).

- `mode: "between", z1: -1.96, z2: 1.96` → `0.9500043`
- `mode: "two_sided", z: 2.5` → p-value `0.01241936`

### confidence_interval

Kinds: `mean_z` (population mean, known sigma), `mean_t` (population mean, unknown sigma; Student-t with df = n − 1), `proportion` (Wilson score interval — stays valid for small samples and extreme proportions, unlike the naive Wald interval). Default confidence level `0.95`.

- `kind: "mean_t", mean: 72.5, sd: 8.2, n: 25` → `[69.115, 75.885]` (t* = 2.0639, df = 24)
- `kind: "proportion", successes: 40, n: 100` → `[0.30940, 0.49800]`
- `kind: "proportion", successes: 0, n: 50` → `[0, 0.07135]`

### event_probability

Operations: `union` (P(A or B); `independent` default true, false treats events as mutually exclusive), `intersection` (independent events), `conditional` (P(A|B) = P(A and B)/P(B)), `bayes` (posterior from prior, true-positive and false-positive rates), `complement`, `at_least_one` (1 − (1−p)^n).

- `operation: "bayes", prior: 0.001, truePositive: 0.99, falsePositive: 0.01` → `0.09016`
- `operation: "at_least_one", p: 0.05, n: 20` → `0.64151`

## Numerics & precision

- Normal cdf: Abramowitz & Stegun 7.1.26 erf fit, absolute error ≤ 1.5e-7 (far beyond published 4-decimal z-tables).
- Normal quantile: Acklam inverse normal, relative error < 1.15e-9.
- Binomial / Poisson: exact pmf via recurrence and cdf by summation (no combinatorics overflow); quantiles by monotone search.
- Student-t: regularized incomplete beta (Lentz continued fraction) + Lanczos log-gamma; cross-checked against published t-table values (2.2281, 2.0452, 1.8125, 1.984).
- Exponential / uniform / geometric: closed forms.
- All validation errors return `valid: false` with an explanatory `error` string; no tool ever throws on bad numeric input.

## Development

```bash
pnpm install
pnpm build
pnpm test   # 37 tests; anchors from test/anchors-probstat.py (independent Python oracle)
pnpm lint
```

## License

MIT
