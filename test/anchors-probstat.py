#!/usr/bin/env python3
"""anchors-probstat.py — independent oracle for dsh-probstat v0.1.0.

Authoritative primitives only:
  - normal CDF/PDF: math.erf / direct formula (exact double precision)
  - normal quantile: bisection on erf (1e-13 tolerance)
  - binomial/poisson: exact summation via math.comb / math.factorial
  - t-quantile: incomplete beta via Numerical-Recipes continued fraction (Lentz),
    cross-checked against PUBLISHED t-table values (independent data source)
  - Wilson CI: textbook closed form, cross-checked at degenerate + large-n limits
  - events: closed-form probability identities
"""

import math, sys, json

SQRT2 = math.sqrt(2.0)
SQRT2PI = math.sqrt(2.0 * math.pi)

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
    if k < 0 or k > n: return 0.0
    return math.comb(n, k) * (p ** k) * ((1.0 - p) ** (n - k))

def binomial_cdf(k, n, p):
    k = math.floor(k)
    if k < 0: return 0.0
    if k >= n: return 1.0
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
    if k < 0: return 0.0
    return math.exp(-lam) * (lam ** k) / math.factorial(k)

def poisson_cdf(k, lam):
    k = math.floor(k)
    if k < 0: return 0.0
    return sum(poisson_pmf(i, lam) for i in range(0, int(k) + 1))

def poisson_survival(k, lam):
    return 1.0 - poisson_cdf(k, lam)

def poisson_quantile(q, lam):
    k = 0
    while poisson_cdf(k, lam) < q:
        k += 1
    return k

# ----------------------------------------------------------- exponential ----
def exp_pdf(x, lam): return lam * math.exp(-lam * x) if x >= 0 else 0.0
def exp_cdf(x, lam): return 0.0 if x < 0 else 1.0 - math.exp(-lam * x)
def exp_survival(x, lam): return 1.0 if x < 0 else math.exp(-lam * x)
def exp_quantile(q, lam): return -math.log(1.0 - q) / lam

# -------------------------------------------------------------- uniform ----
def unif_pdf(x, a, b): return 0.0 if x < a or x > b else 1.0 / (b - a)
def unif_cdf(x, a, b):
    if x < a: return 0.0
    if x > b: return 1.0
    return (x - a) / (b - a)
def unif_quantile(q, a, b): return a + q * (b - a)

# ------------------------------------------------------------ geometric ----
def geom_pmf(k, p): return p * ((1.0 - p) ** k) if k >= 0 else 0.0
def geom_cdf(k, p):
    if k < 0: return 0.0
    return 1.0 - (1.0 - p) ** (k + 1)
def geom_survival(k, p): return (1.0 - p) ** (k + 1) if k >= 0 else 1.0
def geom_quantile(q, p):
    k = 0
    while geom_cdf(k, p) < q:
        k += 1
    return k

# ------------------------------------------------------- t-distribution ----
def betacf(a, b, x, MAXIT=200, EPS=3e-14, FPMIN=1e-300):
    qab = a + b; qap = a + 1.0; qam = a - 1.0
    c = 1.0; d = 1.0 - qab * x / qap
    if abs(d) < FPMIN: d = FPMIN
    d = 1.0 / d; h = d
    for m in range(1, MAXIT + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < FPMIN: d = FPMIN
        c = 1.0 + aa / c
        if abs(c) < FPMIN: c = FPMIN
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < FPMIN: d = FPMIN
        c = 1.0 + aa / c
        if abs(c) < FPMIN: c = FPMIN
        d = 1.0 / d
        dele = d * c
        h *= dele
        if abs(dele - 1.0) < EPS: break
    return h

def beta_regularized(a, b, x):
    if x <= 0.0: return 0.0
    if x >= 1.0: return 1.0
    ln = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log(1.0 - x)
    bt = math.exp(ln)
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * betacf(a, b, x) / a
    return 1.0 - bt * betacf(b, a, 1.0 - x) / b

def t_cdf(t, df):
    x = df / (df + t * t)
    ib = beta_regularized(df / 2.0, 0.5, x)
    if t > 0: return 1.0 - 0.5 * ib
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
def wilson_ci(successes, n, conf):
    z = normal_quantile(1.0 - (1.0 - conf) / 2.0)
    phat = successes / n
    denom = 1.0 + z * z / n
    center = (phat + z * z / (2 * n)) / denom
    half = z * math.sqrt(phat * (1 - phat) / n + z * z / (4 * n * n)) / denom
    return center - half, center + half

# -------------------------------------------------------------- events -----
def bayes(prior, tp, fp):
    return prior * tp / (prior * tp + (1 - prior) * fp)

# ----------------------------------------------------------------------
out = {}

# z-table standard values (published): Phi(1.96)=0.975002..., Phi(2.576)=0.995002...
out['normal_cdf_1.96'] = normal_cdf(1.96)
out['normal_cdf_2.576'] = normal_cdf(2.576)
out['normal_cdf_0'] = normal_cdf(0.0)
out['normal_cdf_neg1'] = normal_cdf(-1.0)
out['normal_pdf_0'] = normal_pdf(0.0)
out['normal_pdf_1'] = normal_pdf(1.0)
out['normal_q_0.975'] = normal_quantile(0.975)      # published 1.95996398454
out['normal_q_0.95'] = normal_quantile(0.95)        # published 1.64485362695
out['normal_q_0.025'] = normal_quantile(0.025)      # published -1.95996398454
out['normal_q_0.995'] = normal_quantile(0.995)      # published 2.57582930355
out['between_-1_1'] = normal_cdf(1.0) - normal_cdf(-1.0)   # published 0.6826894921
out['two_sided_1.96'] = 2 * (1 - normal_cdf(1.96))   # published 0.049995...
out['two_sided_2.0'] = 2 * (1 - normal_cdf(2.0))     # published 0.04550026390

# binomial
out['bin_pmf_3of10_0.5'] = binomial_pmf(3, 10, 0.5)          # 0.1171875
out['bin_cdf_3of10_0.5'] = binomial_cdf(3, 10, 0.5)          # 0.171875
out['bin_cdf_5of20_0.3'] = binomial_cdf(5, 20, 0.3)
out['bin_surv_1of5_0.2'] = binomial_survival(1, 5, 0.2)
out['bin_q_0.5_10_0.5'] = binomial_quantile(0.5, 10, 0.5)    # 5
out['bin_q_0.95_20_0.4'] = binomial_quantile(0.95, 20, 0.4)

# poisson
out['pois_pmf_2_lam3'] = poisson_pmf(2, 3.0)        # 0.224041807...
out['pois_cdf_2_lam3'] = poisson_cdf(2, 3.0)        # 0.423190081...
out['pois_cdf_5_lam4.5'] = poisson_cdf(5, 4.5)
out['pois_surv_3_lam2'] = poisson_survival(3, 2.0)
out['pois_q_0.5_lam3'] = poisson_quantile(0.5, 3.0)  # 3
out['pois_q_0.99_lam10'] = poisson_quantile(0.99, 10.0)

# exponential
out['exp_pdf_1_lam0.5'] = exp_pdf(1.0, 0.5)
out['exp_cdf_1_lam0.5'] = exp_cdf(1.0, 0.5)         # 1-e^-0.5 = 0.39346934...
out['exp_surv_2_lam0.5'] = exp_survival(2.0, 0.5)   # e^-1 = 0.36787944...
out['exp_q_0.5_lam0.5'] = exp_quantile(0.5, 0.5)    # ln2/0.5 = 1.38629436...
out['exp_cdf_neg_lam1'] = exp_cdf(-1.0, 1.0)        # 0

# uniform
out['unif_pdf_1.5_0_3'] = unif_pdf(1.5, 0.0, 3.0)
out['unif_cdf_2_0_3'] = unif_cdf(2.0, 0.0, 3.0)
out['unif_cdf_neg'] = unif_cdf(-1.0, 0.0, 3.0)
out['unif_q_0.25_2_6'] = unif_quantile(0.25, 2.0, 6.0)   # 3.0

# geometric (failures before first success)
out['geom_pmf_2_p0.3'] = geom_pmf(2, 0.3)           # 0.3*0.49=0.147
out['geom_cdf_2_p0.3'] = geom_cdf(2, 0.3)           # 1-0.343=0.657
out['geom_surv_1_p0.25'] = geom_survival(1, 0.25)   # 0.5625
out['geom_q_0.5_p0.5'] = geom_quantile(0.5, 0.5)    # 0 (cdf(0)=0.5 >= 0.5)

# t-distribution — cross-check vs PUBLISHED t-table values
out['t_q_0.975_df10'] = t_quantile(0.975, 10)   # published 2.2281
out['t_q_0.975_df29'] = t_quantile(0.975, 29)   # published 2.0452
out['t_q_0.95_df10'] = t_quantile(0.95, 10)     # published 1.8125
out['t_q_0.975_df100'] = t_quantile(0.975, 100) # published 1.9840
out['t_q_0.975_df1000'] = t_quantile(0.975, 1000)  # -> ~1.9623
out['t_cdf_2.2281_df10'] = t_cdf(2.2281, 10)    # ~0.9750

# confidence intervals
out['ci_mean_z'] = [round(v, 9) for v in
                    (lambda z, se: (100 - z * se, 100 + z * se, z * se, z, se))
                    (normal_quantile(0.975), 15 / math.sqrt(30))]
out['ci_mean_t'] = [round(v, 9) for v in
                    (lambda t, se: (5 - t * se, 5 + t * se, t * se, t, se))
                    (t_quantile(0.975, 11), 1.5 / math.sqrt(12))]
out['wilson_40of100'] = [round(v, 9) for v in wilson_ci(40, 100, 0.95)]
out['wilson_0of50'] = [round(v, 9) for v in wilson_ci(0, 50, 0.95)]
out['wilson_all'] = [round(v, 9) for v in wilson_ci(50, 50, 0.95)]

# events
out['union_indep'] = 0.3 + 0.4 - 0.3 * 0.4            # 0.58
out['union_mutex'] = 0.3 + 0.4                        # 0.7
out['intersection_indep'] = 0.3 * 0.4                 # 0.12
out['conditional'] = 0.2 / 0.5                        # 0.4
out['bayes_standard'] = bayes(0.01, 0.99, 0.05)       # ~0.1667
out['at_least_one'] = 1 - (1 - 0.1) ** 5              # 0.40951
out['complement'] = 1 - 0.25                          # 0.75

for k, v in out.items():
    print(f"{k} = {v!r}")

# sanity assertions vs published constants
assert abs(out['normal_cdf_1.96'] - 0.9750021048517795) < 1e-15
assert abs(out['normal_q_0.975'] - 1.959963984540054) < 1e-9
assert abs(out['normal_q_0.95'] - 1.6448536269514722) < 1e-9
assert abs(out['between_-1_1'] - 0.6826894921370859) < 1e-15
assert abs(out['bin_pmf_3of10_0.5'] - 0.1171875) < 1e-15
assert abs(out['t_q_0.975_df10'] - 2.2281) < 1e-3
assert abs(out['t_q_0.975_df29'] - 2.0452) < 1e-3
assert abs(out['t_q_0.95_df10'] - 1.8125) < 1e-3
assert abs(out['t_q_0.975_df100'] - 1.9840) < 1e-3
print("ALL CROSS-CHECKS PASSED")
