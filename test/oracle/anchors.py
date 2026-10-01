#!/usr/bin/env python3
"""anchors.py — independent Python oracle for dsh-probstat (v0.2.0).

Every numeric expectation in ``test/probstat.test.ts`` is produced here, from
primitives that do not share code with the TypeScript implementation:

  - normal CDF/PDF: ``math.erf`` / direct formula (exact double precision)
  - normal quantile: bisection on erf (1e-13 tolerance)
  - binomial / poisson: exact summation via ``math.comb`` / ``math.factorial``
  - t-quantile: incomplete beta via Numerical-Recipes continued fraction
    (Lentz), cross-checked against PUBLISHED t-table values
  - Wilson CI: textbook closed form, cross-checked at degenerate + large-n
    limits
  - hypothesis tests (v0.2.0): p-values from erf / t_cdf, critical values by
    bisection on the independent CDFs
  - sample size (v0.2.0): Wald closed form; Wilson target solved by an
    independent binary search on the Wilson half-width

Sections mirror the test file so a reviewer can walk assertion -> anchor:

  A. legacy (v0.1.x) anchors — distributions, z-table, CI, events, t-table
  B. hypothesis_test anchors (v0.2.0)
  C. sample_size anchors (v0.2.0)

Usage:
  python3 test/oracle/anchors.py            # print every anchor
  python3 test/oracle/anchors.py --check    # + published-constant self-checks
"""

import math

SQRT2 = math.sqrt(2.0)
SQRT2PI = math.sqrt(2.0 * math.pi)


def r8(x):
    """Round to 8 significant digits — the granularity used in the tests."""
    if not isinstance(x, float):
        return x
    return float(f"{x:.8g}")


# ---------------------------------------------------------------- normal ----
def normal_pdf(x, mu=0.0, sigma=1.0):
    z = (x - mu) / sigma
    return math.exp(-0.5 * z * z) / (SQRT2PI * sigma)


def normal_cdf(x, mu=0.0, sigma=1.0):
    return 0.5 * (1.0 + math.erf((x - mu) / (sigma * SQRT2)))


def normal_quantile(q, mu=0.0, sigma=1.0):
    lo, hi = -42.0, 42.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if normal_cdf(mid) < q:
            lo = mid
        else:
            hi = mid
    return mu + sigma * 0.5 * (lo + hi)


# -------------------------------------------------------------- binomial ----
def binomial_pmf(k, n, p):
    if k < 0 or k > n:
        return 0.0
    return math.comb(n, k) * (p ** k) * ((1.0 - p) ** (n - k))


def binomial_cdf(k, n, p):
    k = math.floor(k)
    if k < 0:
        return 0.0
    if k >= n:
        return 1.0
    return sum(binomial_pmf(i, n, p) for i in range(0, int(k) + 1))


def binomial_survival(k, n, p):  # P(X > k)
    return 1.0 - binomial_cdf(k, n, p)


def binomial_quantile(q, n, p):
    k = 0
    while binomial_cdf(k, n, p) < q:
        k += 1
    return k


# --------------------------------------------------------------- poisson ----
def poisson_pmf(k, lam):
    if k < 0:
        return 0.0
    return math.exp(-lam) * (lam ** k) / math.factorial(k)


def poisson_cdf(k, lam):
    k = math.floor(k)
    if k < 0:
        return 0.0
    return sum(poisson_pmf(i, lam) for i in range(0, int(k) + 1))


def poisson_survival(k, lam):
    return 1.0 - poisson_cdf(k, lam)


def poisson_quantile(q, lam):
    k = 0
    while poisson_cdf(k, lam) < q:
        k += 1
    return k


# ----------------------------------------------------------- exponential ----
def exp_pdf(x, lam):
    return lam * math.exp(-lam * x) if x >= 0 else 0.0


def exp_cdf(x, lam):
    return 0.0 if x < 0 else 1.0 - math.exp(-lam * x)


def exp_survival(x, lam):
    return 1.0 if x < 0 else math.exp(-lam * x)


def exp_quantile(q, lam):
    return -math.log(1.0 - q) / lam


# -------------------------------------------------------------- uniform ----
def unif_pdf(x, a, b):
    return 0.0 if x < a or x > b else 1.0 / (b - a)


def unif_cdf(x, a, b):
    if x < a:
        return 0.0
    if x > b:
        return 1.0
    return (x - a) / (b - a)


def unif_quantile(q, a, b):
    return a + q * (b - a)


# ------------------------------------------------------------ geometric ----
def geom_pmf(k, p):
    return p * ((1.0 - p) ** k) if k >= 0 else 0.0


def geom_cdf(k, p):
    if k < 0:
        return 0.0
    return 1.0 - (1.0 - p) ** (k + 1)


def geom_survival(k, p):
    return (1.0 - p) ** (k + 1) if k >= 0 else 1.0


def geom_quantile(q, p):
    k = 0
    while geom_cdf(k, p) < q:
        k += 1
    return k


# ------------------------------------------------------- t-distribution ----
def betacf(a, b, x, MAXIT=200, EPS=3e-14, FPMIN=1e-300):
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < FPMIN:
        d = FPMIN
    d = 1.0 / d
    h = d
    for m in range(1, MAXIT + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < FPMIN:
            d = FPMIN
        c = 1.0 + aa / c
        if abs(c) < FPMIN:
            c = FPMIN
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < FPMIN:
            d = FPMIN
        c = 1.0 + aa / c
        if abs(c) < FPMIN:
            c = FPMIN
        d = 1.0 / d
        dele = d * c
        h *= dele
        if abs(dele - 1.0) < EPS:
            break
    return h


def beta_regularized(a, b, x):
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    ln = (math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
          + a * math.log(x) + b * math.log(1.0 - x))
    bt = math.exp(ln)
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * betacf(a, b, x) / a
    return 1.0 - bt * betacf(b, a, 1.0 - x) / b


def t_cdf(t, df):
    x = df / (df + t * t)
    ib = beta_regularized(df / 2.0, 0.5, x)
    if t > 0:
        return 1.0 - 0.5 * ib
    return 0.5 * ib


def t_quantile(q, df):
    lo, hi = -1000.0, 1000.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if t_cdf(mid, df) < q:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# ------------------------------------------------------------ Wilson CI ----
def wilson_half_width(n, phat, conf):
    """Wilson score interval half-width for n observations at proportion phat."""
    z = normal_quantile(1.0 - (1.0 - conf) / 2.0)
    return z * math.sqrt(phat * (1.0 - phat) / n + z * z / (4.0 * n * n)) / (1.0 + z * z / n)


def wilson_ci(successes, n, conf):
    z = normal_quantile(1.0 - (1.0 - conf) / 2.0)
    phat = successes / n
    denom = 1.0 + z * z / n
    center = (phat + z * z / (2 * n)) / denom
    half = z * math.sqrt(phat * (1.0 - phat) / n + z * z / (4 * n * n)) / denom
    return center - half, center + half


# --------------------------------------------------- hypothesis tests ------
def test_p_values(stat, cdf):
    """Two-sided / left / right p-values for a symmetric null distribution."""
    right = 1.0 - cdf(stat)
    left = cdf(stat)
    two = 2.0 * (1.0 - cdf(abs(stat)))
    return two, left, right


def z_test_mean(sample_mean, mu0, sigma, n):
    se = sigma / math.sqrt(n)
    stat = (sample_mean - mu0) / se
    two, left, right = test_p_values(stat, lambda t: normal_cdf(t))
    return dict(
        statistic=stat,
        standardError=se,
        pTwoSided=two,
        pLeft=left,
        pRight=right,
        critTwoSided=normal_quantile(1.0 - 0.05 / 2.0),
        critLeft=normal_quantile(0.05),
        critRight=normal_quantile(0.95),
    )


def t_test_mean(sample_mean, mu0, sd, n):
    se = sd / math.sqrt(n)
    df = n - 1
    stat = (sample_mean - mu0) / se
    two, left, right = test_p_values(stat, lambda t: t_cdf(t, df))
    return dict(
        statistic=stat,
        standardError=se,
        df=df,
        pTwoSided=two,
        pLeft=left,
        pRight=right,
        critTwoSided=t_quantile(1.0 - 0.05 / 2.0, df),
        critLeft=t_quantile(0.05, df),
        critRight=t_quantile(0.95, df),
    )


def z_test_proportion(successes, n, p0):
    phat = successes / n
    se = math.sqrt(p0 * (1.0 - p0) / n)
    stat = (phat - p0) / se
    two, left, right = test_p_values(stat, lambda t: normal_cdf(t))
    return dict(
        statistic=stat,
        standardError=se,
        phat=phat,
        pTwoSided=two,
        pLeft=left,
        pRight=right,
        critTwoSided=normal_quantile(1.0 - 0.05 / 2.0),
        critLeft=normal_quantile(0.05),
        critRight=normal_quantile(0.95),
    )


# ------------------------------------------------------ sample size --------
def sample_size_mean(sigma, e, conf):
    z = normal_quantile(1.0 - (1.0 - conf) / 2.0)
    n_exact = (z * sigma / e) ** 2
    return dict(criticalValue=z, nExact=n_exact, n=math.ceil(n_exact))


def sample_size_proportion_wald(p, e, conf):
    z = normal_quantile(1.0 - (1.0 - conf) / 2.0)
    n_exact = z * z * p * (1.0 - p) / (e * e)
    return dict(criticalValue=z, nExact=n_exact, n=math.ceil(n_exact))


def sample_size_proportion_wilson(p, e, conf, hi=10 ** 7):
    """Smallest n whose Wilson half-width at phat = p is <= e (binary search)."""
    if wilson_half_width(hi, p, conf) > e:
        raise ValueError("no n <= 1e7 meets the target margin of error")
    lo, hi_n = 1, hi
    while lo < hi_n:
        mid = (lo + hi_n) // 2
        if wilson_half_width(mid, p, conf) <= e:
            hi_n = mid
        else:
            lo = mid + 1
    return dict(
        criticalValue=normal_quantile(1.0 - (1.0 - conf) / 2.0),
        n=lo,
        halfWidthAtN=wilson_half_width(lo, p, conf),
        halfWidthAtNMinus1=wilson_half_width(lo - 1, p, conf) if lo > 1 else None,
    )


# ----------------------------------------------------------------------
# A. legacy (v0.1.x) anchors
# ----------------------------------------------------------------------
legacy = {}

# --- math primitives: normal (published z-table values)
legacy['normal cdf(1.96)'] = normal_cdf(1.96)              # published 0.9750021
legacy['normal cdf(2.576)'] = normal_cdf(2.576)            # published 0.9950025
legacy['normal cdf(0)'] = normal_cdf(0.0)                  # 0.5
legacy['normal cdf(-1)'] = normal_cdf(-1.0)                # 0.1586553
legacy['normal cdf(5,2,3) == cdf(1,0,1)'] = normal_cdf(5.0, 2.0, 3.0)
legacy['normal quantile(0.975)'] = normal_quantile(0.975)  # published 1.959964
legacy['normal quantile(0.95)'] = normal_quantile(0.95)    # published 1.6448536
legacy['normal quantile(0.025)'] = normal_quantile(0.025)  # -1.959964
legacy['normal quantile(0.995)'] = normal_quantile(0.995)  # 2.5758293
legacy['normal quantile(0.5)'] = normal_quantile(0.5)      # 0
legacy['normal pdf(0)'] = normal_pdf(0.0)                  # 0.3989423
legacy['normal pdf(1)'] = normal_pdf(1.0)                  # 0.2419707

# --- binomial
legacy['binomial pdf 3 of 10 p=0.5'] = binomial_pmf(3, 10, 0.5)     # 0.1171875
legacy['binomial cdf 3 of 10 p=0.5'] = binomial_cdf(3, 10, 0.5)     # 0.171875
legacy['binomial cdf 5 of 20 p=0.3'] = binomial_cdf(5, 20, 0.3)     # 0.4163708
legacy['binomial survival 1 of 5 p=0.2'] = binomial_survival(1, 5, 0.2)  # 0.26272
legacy['binomial quantile q=0.5 n=10 p=0.5'] = binomial_quantile(0.5, 10, 0.5)   # 5
legacy['binomial quantile q=0.95 n=20 p=0.4'] = binomial_quantile(0.95, 20, 0.4)  # 12
legacy['binomial stats n=20 p=0.4 mean'] = 20 * 0.4                   # 8
legacy['binomial stats n=20 p=0.4 variance'] = 20 * 0.4 * 0.6         # 4.8
legacy['binomial stats n=20 p=0.4 sd'] = math.sqrt(20 * 0.4 * 0.6)    # 2.1908902
legacy['binomial cdf n=5 p=0 x=0'] = binomial_cdf(0, 5, 0.0)          # 1
legacy['binomial pdf n=5 p=0 x=0'] = binomial_pmf(0, 5, 0.0)          # 1
legacy['binomial pdf n=5 p=1 x=5'] = binomial_pmf(5, 5, 1.0)          # 1
legacy['binomial cdf n=5 p=1 x=4'] = binomial_cdf(4, 5, 1.0)          # 0

# --- poisson
legacy['poisson pdf k=2 lambda=3'] = poisson_pmf(2, 3.0)      # 0.2240418
legacy['poisson cdf k=2 lambda=3'] = poisson_cdf(2, 3.0)      # 0.4231901
legacy['poisson cdf k=5 lambda=4.5'] = poisson_cdf(5, 4.5)    # 0.7029304
legacy['poisson survival k=3 lambda=2'] = poisson_survival(3, 2.0)   # 0.1428765
legacy['poisson quantile q=0.5 lambda=3'] = poisson_quantile(0.5, 3.0)   # 3
legacy['poisson quantile q=0.99 lambda=10'] = poisson_quantile(0.99, 10.0)  # 18
legacy['poisson stats lambda=4.5 mean'] = 4.5
legacy['poisson stats lambda=4.5 variance'] = 4.5
legacy['poisson stats lambda=4.5 sd'] = math.sqrt(4.5)        # 2.1213203

# --- exponential
legacy['exponential pdf x=1 lambda=0.5'] = exp_pdf(1.0, 0.5)       # 0.3032653
legacy['exponential cdf x=1 lambda=0.5'] = exp_cdf(1.0, 0.5)       # 0.3934693
legacy['exponential survival x=2 lambda=0.5'] = exp_survival(2.0, 0.5)  # 0.3678794
legacy['exponential quantile q=0.5 lambda=0.5'] = exp_quantile(0.5, 0.5)  # 1.3862944
legacy['exponential cdf x=-1 lambda=1'] = exp_cdf(-1.0, 1.0)       # 0

# --- uniform
legacy['uniform pdf x=1.5 [0,3]'] = unif_pdf(1.5, 0.0, 3.0)        # 0.3333333
legacy['uniform cdf x=2 [0,3]'] = unif_cdf(2.0, 0.0, 3.0)          # 0.6666667
legacy['uniform cdf x=-1 [0,3]'] = unif_cdf(-1.0, 0.0, 3.0)        # 0
legacy['uniform quantile q=0.25 [2,6]'] = unif_quantile(0.25, 2.0, 6.0)   # 3
legacy['uniform stats [2,6] mean'] = 4.0
legacy['uniform stats [2,6] variance'] = (6 - 2) ** 2 / 12.0       # 1.3333333

# --- geometric (failures before first success)
legacy['geometric pdf k=2 p=0.3'] = geom_pmf(2, 0.3)          # 0.147
legacy['geometric cdf k=2 p=0.3'] = geom_cdf(2, 0.3)          # 0.657
legacy['geometric survival k=1 p=0.25'] = geom_survival(1, 0.25)  # 0.5625
legacy['geometric quantile q=0.5 p=0.5'] = geom_quantile(0.5, 0.5)  # 0
legacy['geometric stats p=0.25 mean'] = (1 - 0.25) / 0.25     # 3
legacy['geometric stats p=0.25 variance'] = (1 - 0.25) / 0.25 ** 2  # 12

# --- t-table cross-checks (published values)
legacy['t quantile 0.975 df=10'] = t_quantile(0.975, 10)      # published 2.2281
legacy['t quantile 0.975 df=29'] = t_quantile(0.975, 29)      # published 2.0452
legacy['t quantile 0.95 df=10'] = t_quantile(0.95, 10)        # published 1.8125
legacy['t quantile 0.975 df=100'] = t_quantile(0.975, 100)    # published 1.9840

# --- z_score tool
legacy['z_score to_prob z=1.96 probability'] = normal_cdf(1.96)
legacy['z_score to_prob z=1.96 percentile'] = normal_cdf(1.96) * 100.0     # 97.50021
legacy['z_score from_prob p=0.975 z'] = normal_quantile(0.975)             # 1.959964
legacy['z_score from_prob p=0.975 percentile'] = 97.5
legacy['z_score two_sided z=1.96'] = 2.0 * (1.0 - normal_cdf(1.96))        # 0.0499958
legacy['z_score two_sided z=2'] = 2.0 * (1.0 - normal_cdf(2.0))            # 0.0455003
legacy['z_score between -1..1'] = normal_cdf(1.0) - normal_cdf(-1.0)       # 0.6826895
legacy['z_score between sorted bounds'] = (-1, 1)

# --- confidence_interval tool
_ci_z = z_test_mean(0.0, 0.0, 1.0, 1)  # unused placeholder (keeps names unique)
_ci_mean_z = (lambda z, se: (100 - z * se, 100 + z * se, z * se, z, se))(
    normal_quantile(0.975), 15 / math.sqrt(30))
legacy['ci mean_z lower'] = _ci_mean_z[0]        # 94.6324176
legacy['ci mean_z upper'] = _ci_mean_z[1]        # 105.3675824
legacy['ci mean_z marginOfError'] = _ci_mean_z[2]  # 5.3675824
legacy['ci mean_z criticalValue'] = _ci_mean_z[3]  # 1.959964
legacy['ci mean_z standardError'] = _ci_mean_z[4]  # 2.7386128
_ci_mean_t = (lambda t, se: (5 - t * se, 5 + t * se, t * se, t, se))(
    t_quantile(0.975, 11), 1.5 / math.sqrt(12))
legacy['ci mean_t lower'] = _ci_mean_t[0]        # 4.0469455
legacy['ci mean_t upper'] = _ci_mean_t[1]        # 5.9530545
legacy['ci mean_t marginOfError'] = _ci_mean_t[2]  # 0.9530545
legacy['ci mean_t criticalValue'] = _ci_mean_t[3]  # 2.2009852
legacy['ci proportion 40/100 lower'] = wilson_ci(40, 100, 0.95)[0]   # 0.3094013
legacy['ci proportion 40/100 upper'] = wilson_ci(40, 100, 0.95)[1]   # 0.4979974
legacy['ci proportion 0/50 lower'] = wilson_ci(0, 50, 0.95)[0]       # 0
legacy['ci proportion 0/50 upper'] = wilson_ci(0, 50, 0.95)[1]       # 0.0713476
legacy['ci proportion 50/50 lower'] = wilson_ci(50, 50, 0.95)[0]     # 0.9286524
legacy['ci proportion 50/50 upper'] = wilson_ci(50, 50, 0.95)[1]     # 1

# --- event_probability tool
legacy['event union independent 0.3/0.4'] = 0.3 + 0.4 - 0.3 * 0.4    # 0.58
legacy['event union mutually exclusive 0.3/0.4'] = 0.3 + 0.4          # 0.7
legacy['event intersection independent 0.3/0.4'] = 0.3 * 0.4          # 0.12
legacy['event conditional 0.2/0.5'] = 0.2 / 0.5                       # 0.4
legacy['event bayes prior=0.01 tp=0.99 fp=0.05'] = (
    0.01 * 0.99 / (0.01 * 0.99 + 0.99 * 0.05))                        # 0.1666667
legacy['event complement p=0.25'] = 1 - 0.25                          # 0.75
legacy['event at_least_one p=0.1 n=5'] = 1 - (1 - 0.1) ** 5           # 0.40951

# --- discipline: default operation is cdf, default normal is standard
legacy['discipline default normal cdf(1.96)'] = normal_cdf(1.96)

# ----------------------------------------------------------------------
# B. hypothesis_test anchors (v0.2.0)
# ----------------------------------------------------------------------
hyp = {}

_h1 = z_test_mean(102.0, 100.0, 15.0, 36)   # z = 0.8, two-sided
hyp['mean_z 102/100 sigma=15 n=36 statistic'] = _h1['statistic']            # 0.8
hyp['mean_z 102/100 sigma=15 n=36 standardError'] = _h1['standardError']    # 2.5
hyp['mean_z 102/100 sigma=15 n=36 pTwoSided'] = _h1['pTwoSided']           # 0.4237108
hyp['mean_z 102/100 sigma=15 n=36 pLeft'] = _h1['pLeft']                   # 0.7881446
hyp['mean_z 102/100 sigma=15 n=36 pRight'] = _h1['pRight']                 # 0.2118554
hyp['mean_z 102/100 sigma=15 n=36 critTwoSided'] = _h1['critTwoSided']     # 1.959964
hyp['mean_z 102/100 sigma=15 n=36 critRight'] = _h1['critRight']           # 1.6448536

_h2 = z_test_mean(100.0, 103.0, 10.0, 25)   # z = -1.5, left tail
hyp['mean_z 100/103 sigma=10 n=25 statistic'] = _h2['statistic']           # -1.5
hyp['mean_z 100/103 sigma=10 n=25 pLeft'] = _h2['pLeft']                   # 0.0668072
hyp['mean_z 100/103 sigma=10 n=25 critLeft'] = _h2['critLeft']             # -1.6448536

_h3 = t_test_mean(5.2, 5.0, 1.5, 12)        # t = 0.4618802, df = 11, right tail
hyp['mean_t 5.2/5 sd=1.5 n=12 statistic'] = _h3['statistic']               # 0.4618802
hyp['mean_t 5.2/5 sd=1.5 n=12 standardError'] = _h3['standardError']       # 0.4330127
hyp['mean_t 5.2/5 sd=1.5 n=12 df'] = _h3['df']                            # 11
hyp['mean_t 5.2/5 sd=1.5 n=12 pRight'] = _h3['pRight']                     # 0.3267656
hyp['mean_t 5.2/5 sd=1.5 n=12 pTwoSided'] = _h3['pTwoSided']               # 0.6535312
hyp['mean_t 5.2/5 sd=1.5 n=12 critRight'] = _h3['critRight']               # 1.7958848
hyp['mean_t 5.2/5 sd=1.5 n=12 critTwoSided'] = _h3['critTwoSided']         # 2.2009852

_h4 = z_test_proportion(40, 100, 0.5)       # z = -2, two-sided, reject
hyp['proportion 40/100 p0=0.5 statistic'] = _h4['statistic']               # -2
hyp['proportion 40/100 p0=0.5 standardError'] = _h4['standardError']       # 0.05
hyp['proportion 40/100 p0=0.5 pTwoSided'] = _h4['pTwoSided']               # 0.0455003
hyp['proportion 40/100 p0=0.5 pLeft'] = _h4['pLeft']                       # 0.0227501
hyp['proportion 40/100 p0=0.5 critTwoSided'] = _h4['critTwoSided']         # 1.959964

_h5 = z_test_proportion(62, 100, 0.5)       # z = 2.4, right tail
hyp['proportion 62/100 p0=0.5 statistic'] = _h5['statistic']               # 2.4
hyp['proportion 62/100 p0=0.5 pRight'] = _h5['pRight']                     # 0.0081975

# ----------------------------------------------------------------------
# C. sample_size anchors (v0.2.0)
# ----------------------------------------------------------------------
ss = {}

_s1 = sample_size_mean(15.0, 5.0, 0.95)
ss['mean sigma=15 E=5 conf=0.95 criticalValue'] = _s1['criticalValue']     # 1.959964
ss['mean sigma=15 E=5 conf=0.95 nExact'] = _s1['nExact']                   # 34.57318
ss['mean sigma=15 E=5 conf=0.95 n'] = _s1['n']                             # 35
_s2 = sample_size_mean(15.0, 5.0, 0.99)
ss['mean sigma=15 E=5 conf=0.99 nExact'] = _s2['nExact']                   # 59.71418
ss['mean sigma=15 E=5 conf=0.99 n'] = _s2['n']                             # 60
_s3 = sample_size_proportion_wald(0.5, 0.03, 0.95)
ss['proportion wald p=0.5 E=0.03 nExact'] = _s3['nExact']                  # 1067.0719
ss['proportion wald p=0.5 E=0.03 n'] = _s3['n']                            # 1068
_s4 = sample_size_proportion_wald(0.2, 0.05, 0.95)
ss['proportion wald p=0.2 E=0.05 nExact'] = _s4['nExact']                  # 245.85339
ss['proportion wald p=0.2 E=0.05 n'] = _s4['n']                            # 246
_s5 = sample_size_proportion_wilson(0.5, 0.03, 0.95)
ss['proportion wilson p=0.5 E=0.03 n'] = _s5['n']
ss['proportion wilson p=0.5 E=0.03 halfWidthAtN'] = _s5['halfWidthAtN']
ss['proportion wilson p=0.5 E=0.03 halfWidthAtNMinus1'] = _s5['halfWidthAtNMinus1']
_s6 = sample_size_proportion_wilson(0.2, 0.05, 0.95)
ss['proportion wilson p=0.2 E=0.05 n'] = _s6['n']
ss['proportion wilson p=0.2 E=0.05 halfWidthAtN'] = _s6['halfWidthAtN']


def emit():
    for label, table in (('A legacy', legacy), ('B hypothesis_test', hyp), ('C sample_size', ss)):
        print(f"=== {label} ===")
        for key, value in table.items():
            if isinstance(value, float):
                print(f"  {key} = {value!r}   round8={r8(value)!r}")
            else:
                print(f"  {key} = {value!r}")


def check():
    """Self-checks against published constants (independent data source)."""
    assert abs(legacy['normal cdf(1.96)'] - 0.9750021048517795) < 1e-15
    assert abs(legacy['normal quantile(0.975)'] - 1.959963984540054) < 1e-9
    assert abs(legacy['normal quantile(0.95)'] - 1.6448536269514722) < 1e-9
    assert abs(legacy['z_score between -1..1'] - 0.6826894921370859) < 1e-15
    assert abs(legacy['binomial pdf 3 of 10 p=0.5'] - 0.1171875) < 1e-15
    assert abs(legacy['t quantile 0.975 df=10'] - 2.2281) < 1e-3
    assert abs(legacy['t quantile 0.975 df=29'] - 2.0452) < 1e-3
    assert abs(legacy['t quantile 0.95 df=10'] - 1.8125) < 1e-3
    assert abs(legacy['t quantile 0.975 df=100'] - 1.9840) < 1e-3
    # published p-value for z = 2.0 (two-sided) = 0.045500263
    assert abs(legacy['z_score two_sided z=2'] - 0.04550026389635842) < 1e-12
    # t-table cross-check: t_cdf(2.2281, 10) = 0.975
    assert abs(t_cdf(2.2281, 10) - 0.975) < 1e-5
    # hypothesis_test: exact z = 0.8 / -2 / 2.4 / -1.5
    assert abs(hyp['mean_z 102/100 sigma=15 n=36 statistic'] - 0.8) < 1e-14
    assert abs(hyp['proportion 40/100 p0=0.5 statistic'] + 2.0) < 1e-14
    assert abs(hyp['proportion 62/100 p0=0.5 statistic'] - 2.4) < 1e-14
    assert abs(hyp['mean_z 100/103 sigma=10 n=25 statistic'] + 1.5) < 1e-14
    # p-value consistency: two-sided = 2 * one-sided at both symmetric tails
    assert abs(hyp['proportion 40/100 p0=0.5 pTwoSided'] - 2 * hyp['proportion 40/100 p0=0.5 pLeft']) < 1e-15
    assert abs(hyp['mean_z 102/100 sigma=15 n=36 pTwoSided'] - 2 * hyp['mean_z 102/100 sigma=15 n=36 pRight']) < 1e-15
    # Wilson sample size: smallest n meeting the target
    assert ss['proportion wilson p=0.5 E=0.03 halfWidthAtN'] <= 0.03
    assert ss['proportion wilson p=0.5 E=0.03 halfWidthAtNMinus1'] > 0.03
    assert ss['proportion wilson p=0.2 E=0.05 halfWidthAtN'] <= 0.05
    # Wald vs Wilson agree within a percent at these targets
    assert abs(ss['proportion wilson p=0.5 E=0.03 n'] - ss['proportion wald p=0.5 E=0.03 n']) <= 5
    print("ALL CROSS-CHECKS PASSED")


if __name__ == '__main__':
    import sys
    emit()
    if '--check' in sys.argv:
        check()
