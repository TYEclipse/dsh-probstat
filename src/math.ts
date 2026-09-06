/**
 * dsh-probstat/math — deterministic probability math primitives.
 *
 * Zero runtime dependencies, pure arithmetic:
 *   - normal:  pdf / cdf (A&S 7.1.26 erf fit, |err| ≤ 1.5e-7) / quantile
 *              (Acklam inverse normal, relative error < 1.15e-9)
 *   - binomial / poisson: exact pmf via recurrence, cdf by summation,
 *     quantile by monotone search (no combinatorics overflow)
 *   - exponential / uniform / geometric: closed forms
 *   - Student-t: cdf via regularized incomplete beta (Lentz continued
 *     fraction) + Lanczos log-gamma; quantile by bisection
 *
 * @module dsh-probstat/math
 */

const SQRT2PI = Math.sqrt(2 * Math.PI)

/* ------------------------------------------------------------------ */
/* normal                                                              */
/* ------------------------------------------------------------------ */

/** Error function approximation (Abramowitz & Stegun 7.1.26, |err| ≤ 1.5e-7). */
function erf(x: number): number {
  if (x === 0) return 0
  const sign = x < 0 ? -1 : 1
  const ax = Math.abs(x)
  const t = 1 / (1 + 0.3275911 * ax)
  const y = 1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-ax * ax)
  return sign * y
}

/** Normal probability density function. */
export function normalPdf(x: number, mean: number, sd: number): number {
  const z = (x - mean) / sd
  return Math.exp(-0.5 * z * z) / (SQRT2PI * sd)
}

/** Normal cumulative distribution function P(X <= x). */
export function normalCdf(x: number, mean: number, sd: number): number {
  return 0.5 * (1 + erf((x - mean) / (sd * Math.SQRT2)))
}

/* Acklam inverse-normal coefficients (Peter J. Acklam). */
const A0 = -3.969683028665376e1
const A1 = 2.209460984245205e2
const A2 = -2.759285104469687e2
const A3 = 1.38357751867269e2
const A4 = -3.066479806614716e1
const A5 = 2.506628277459239e0
const B0 = -5.447609879822406e1
const B1 = 1.615858368580409e2
const B2 = -1.556989798598866e2
const B3 = 6.680131188771972e1
const B4 = -1.328068155288572e1
const C0 = -7.784894002430293e-3
const C1 = -3.223964580411365e-1
const C2 = -2.400758277161838e0
const C3 = -2.549732539343734e0
const C4 = 4.374664141464968e0
const C5 = 2.938163982698783e0
const D0 = 7.784695709041462e-3
const D1 = 3.224671290700398e-1
const D2 = 2.445134137142996e0
const D3 = 3.754408661907416e0
const LOW = 0.02425

/** Normal quantile function (inverse CDF); p must lie strictly in (0, 1). */
export function normalQuantile(p: number, mean: number, sd: number): number {
  let z: number
  if (p < LOW) {
    const q = Math.sqrt(-2 * Math.log(p))
    z = (((((C0 * q + C1) * q + C2) * q + C3) * q + C4) * q + C5)
      / ((((D0 * q + D1) * q + D2) * q + D3) * q + 1)
  } else if (p <= 1 - LOW) {
    const q = p - 0.5
    const r = q * q
    z = (((((A0 * r + A1) * r + A2) * r + A3) * r + A4) * r + A5) * q
      / (((((B0 * r + B1) * r + B2) * r + B3) * r + B4) * r + 1)
  } else {
    const q = Math.sqrt(-2 * Math.log(1 - p))
    z = -(((((C0 * q + C1) * q + C2) * q + C3) * q + C4) * q + C5)
      / ((((D0 * q + D1) * q + D2) * q + D3) * q + 1)
  }
  return mean + sd * z
}

/* ------------------------------------------------------------------ */
/* binomial                                                            */
/* ------------------------------------------------------------------ */

/** Binomial probability mass function P(X = k), exact via recurrence. */
export function binomialPmf(k: number, n: number, p: number): number {
  if (!Number.isInteger(k) || k < 0 || k > n) return 0
  if (p === 0) return k === 0 ? 1 : 0
  if (p === 1) return k === n ? 1 : 0
  let b = Math.pow(1 - p, n)
  for (let i = 0; i < k; i++) {
    b = (b * (n - i) * p) / ((i + 1) * (1 - p))
  }
  return b
}

/** Binomial cumulative distribution function P(X <= k). */
export function binomialCdf(k: number, n: number, p: number): number {
  if (k < 0) return 0
  if (k >= n) return 1
  const top = Math.floor(k)
  let sum = 0
  let b = Math.pow(1 - p, n)
  if (p === 0) return 1
  if (p === 1) return top >= n ? 1 : 0
  for (let i = 0; i <= top; i++) {
    sum += b
    b = (b * (n - i) * p) / ((i + 1) * (1 - p))
  }
  return sum
}

/** Binomial quantile: smallest k with P(X <= k) >= q. */
export function binomialQuantile(q: number, n: number, p: number): number {
  if (p === 1) return n
  let k = 0
  let cdf = binomialCdf(0, n, p)
  while (cdf < q) {
    k++
    cdf = binomialCdf(k, n, p)
  }
  return k
}

/* ------------------------------------------------------------------ */
/* poisson                                                             */
/* ------------------------------------------------------------------ */

/** Poisson probability mass function P(X = k). */
export function poissonPmf(k: number, lambda: number): number {
  if (!Number.isInteger(k) || k < 0) return 0
  if (lambda === 0) return k === 0 ? 1 : 0
  let term = Math.exp(-lambda)
  for (let i = 0; i < k; i++) {
    term = (term * lambda) / (i + 1)
  }
  return term
}

/** Poisson cumulative distribution function P(X <= k). */
export function poissonCdf(k: number, lambda: number): number {
  if (k < 0) return 0
  const top = Math.floor(k)
  let sum = 0
  let term = Math.exp(-lambda)
  for (let i = 0; i <= top; i++) {
    sum += term
    term = (term * lambda) / (i + 1)
  }
  return sum
}

/** Poisson quantile: smallest k with P(X <= k) >= q. */
export function poissonQuantile(q: number, lambda: number): number {
  let k = 0
  let cdf = poissonCdf(0, lambda)
  while (cdf < q) {
    k++
    cdf = poissonCdf(k, lambda)
  }
  return k
}

/* ------------------------------------------------------------------ */
/* exponential / uniform / geometric                                   */
/* ------------------------------------------------------------------ */

/** Exponential pdf with rate lambda. */
export function exponentialPdf(x: number, lambda: number): number {
  return x < 0 ? 0 : lambda * Math.exp(-lambda * x)
}

/** Exponential cdf P(X <= x) with rate lambda. */
export function exponentialCdf(x: number, lambda: number): number {
  return x < 0 ? 0 : 1 - Math.exp(-lambda * x)
}

/** Exponential survival P(X > x) with rate lambda. */
export function exponentialSurvival(x: number, lambda: number): number {
  return x < 0 ? 1 : Math.exp(-lambda * x)
}

/** Exponential quantile. */
export function exponentialQuantile(q: number, lambda: number): number {
  return -Math.log(1 - q) / lambda
}

/** Uniform pdf on [a, b]. */
export function uniformPdf(x: number, a: number, b: number): number {
  return x < a || x > b ? 0 : 1 / (b - a)
}

/** Uniform cdf on [a, b]. */
export function uniformCdf(x: number, a: number, b: number): number {
  if (x < a) return 0
  if (x > b) return 1
  return (x - a) / (b - a)
}

/** Uniform quantile on [a, b]. */
export function uniformQuantile(q: number, a: number, b: number): number {
  return a + q * (b - a)
}

/**
 * Geometric pmf — support convention: k = number of FAILURES before the
 * first success, k >= 0. P(X = k) = p * (1 - p)^k.
 */
export function geometricPmf(k: number, p: number): number {
  if (!Number.isInteger(k) || k < 0) return 0
  return p * Math.pow(1 - p, k)
}

/** Geometric cdf P(X <= k), failures-before-success convention. */
export function geometricCdf(k: number, p: number): number {
  if (k < 0) return 0
  return 1 - Math.pow(1 - p, Math.floor(k) + 1)
}

/** Geometric survival P(X > k), failures-before-success convention. */
export function geometricSurvival(k: number, p: number): number {
  if (k < 0) return 1
  return Math.pow(1 - p, Math.floor(k) + 1)
}

/** Geometric quantile: smallest k with P(X <= k) >= q. */
export function geometricQuantile(q: number, p: number): number {
  let k = 0
  let cdf = geometricCdf(0, p)
  while (cdf < q) {
    k++
    cdf = geometricCdf(k, p)
  }
  return k
}

/* ------------------------------------------------------------------ */
/* Student-t distribution                                              */
/* ------------------------------------------------------------------ */

/* Lanczos coefficients (g = 7, n = 9), rounded to 15 significant digits. */
const LANCZOS = [
  0.999999999999810,
  676.520368121885,
  -1259.139216722403,
  771.323428777653,
  -176.615029162141,
  12.507343278687,
  -0.138571095266,
  9.98436957801957e-6,
  1.50563273514931e-7,
]

/** Natural log of the gamma function (Lanczos approximation). */
export function logGamma(z: number): number {
  if (z < 0.5) {
    return Math.log(Math.PI) - Math.log(Math.sin(Math.PI * z)) - logGamma(1 - z)
  }
  const x = z - 1
  let sum = 0
  for (const [i, coef] of LANCZOS.entries()) {
    sum += i === 0 ? coef : coef / (x + i)
  }
  const t = x + 7.5
  return 0.5 * Math.log(2 * Math.PI) + (x + 0.5) * Math.log(t) - t + Math.log(sum)
}

/**
 * Continued fraction for the incomplete beta function
 * (Numerical Recipes betacf, modified Lentz's method).
 */
function betacf(a: number, b: number, x: number): number {
  const MAXIT = 200
  const EPS = 3e-14
  const FPMIN = 1e-300
  const qab = a + b
  const qap = a + 1
  const qam = a - 1
  let c = 1
  let d = 1 - (qab * x) / qap
  if (Math.abs(d) < FPMIN) d = FPMIN
  d = 1 / d
  let h = d
  for (let m = 1; m <= MAXIT; m++) {
    const m2 = 2 * m
    let aa = (m * (b - m) * x) / ((qam + m2) * (a + m2))
    d = 1 + aa * d
    if (Math.abs(d) < FPMIN) d = FPMIN
    c = 1 + aa / c
    if (Math.abs(c) < FPMIN) c = FPMIN
    d = 1 / d
    h *= d * c
    aa = (-(a + m) * (qab + m) * x) / ((a + m2) * (qap + m2))
    d = 1 + aa * d
    if (Math.abs(d) < FPMIN) d = FPMIN
    c = 1 + aa / c
    if (Math.abs(c) < FPMIN) c = FPMIN
    d = 1 / d
    const dele = d * c
    h *= dele
    if (Math.abs(dele - 1) < EPS) break
  }
  return h
}

/** Regularized incomplete beta I_x(a, b). */
function betaRegularized(a: number, b: number, x: number): number {
  if (x <= 0) return 0
  if (x >= 1) return 1
  const ln = logGamma(a + b) - logGamma(a) - logGamma(b) + a * Math.log(x) + b * Math.log(1 - x)
  const bt = Math.exp(ln)
  if (x < (a + 1) / (a + b + 2)) {
    return (bt * betacf(a, b, x)) / a
  }
  return 1 - (bt * betacf(b, a, 1 - x)) / b
}

/** Student-t cumulative distribution function with df degrees of freedom. */
export function tCdf(t: number, df: number): number {
  const x = df / (df + t * t)
  const ib = betaRegularized(df / 2, 0.5, x)
  return t > 0 ? 1 - 0.5 * ib : 0.5 * ib
}

/** Student-t quantile function (bisection on tCdf); q in (0, 1), df > 0. */
export function tQuantile(q: number, df: number): number {
  let lo = -1000
  let hi = 1000
  for (let i = 0; i < 200; i++) {
    const mid = 0.5 * (lo + hi)
    if (tCdf(mid, df) < q) {
      lo = mid
    } else {
      hi = mid
    }
  }
  return 0.5 * (lo + hi)
}
