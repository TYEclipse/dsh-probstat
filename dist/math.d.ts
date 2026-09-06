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
/** Normal probability density function. */
export declare function normalPdf(x: number, mean: number, sd: number): number;
/** Normal cumulative distribution function P(X <= x). */
export declare function normalCdf(x: number, mean: number, sd: number): number;
/** Normal quantile function (inverse CDF); p must lie strictly in (0, 1). */
export declare function normalQuantile(p: number, mean: number, sd: number): number;
/** Binomial probability mass function P(X = k), exact via recurrence. */
export declare function binomialPmf(k: number, n: number, p: number): number;
/** Binomial cumulative distribution function P(X <= k). */
export declare function binomialCdf(k: number, n: number, p: number): number;
/** Binomial quantile: smallest k with P(X <= k) >= q. */
export declare function binomialQuantile(q: number, n: number, p: number): number;
/** Poisson probability mass function P(X = k). */
export declare function poissonPmf(k: number, lambda: number): number;
/** Poisson cumulative distribution function P(X <= k). */
export declare function poissonCdf(k: number, lambda: number): number;
/** Poisson quantile: smallest k with P(X <= k) >= q. */
export declare function poissonQuantile(q: number, lambda: number): number;
/** Exponential pdf with rate lambda. */
export declare function exponentialPdf(x: number, lambda: number): number;
/** Exponential cdf P(X <= x) with rate lambda. */
export declare function exponentialCdf(x: number, lambda: number): number;
/** Exponential survival P(X > x) with rate lambda. */
export declare function exponentialSurvival(x: number, lambda: number): number;
/** Exponential quantile. */
export declare function exponentialQuantile(q: number, lambda: number): number;
/** Uniform pdf on [a, b]. */
export declare function uniformPdf(x: number, a: number, b: number): number;
/** Uniform cdf on [a, b]. */
export declare function uniformCdf(x: number, a: number, b: number): number;
/** Uniform quantile on [a, b]. */
export declare function uniformQuantile(q: number, a: number, b: number): number;
/**
 * Geometric pmf — support convention: k = number of FAILURES before the
 * first success, k >= 0. P(X = k) = p * (1 - p)^k.
 */
export declare function geometricPmf(k: number, p: number): number;
/** Geometric cdf P(X <= k), failures-before-success convention. */
export declare function geometricCdf(k: number, p: number): number;
/** Geometric survival P(X > k), failures-before-success convention. */
export declare function geometricSurvival(k: number, p: number): number;
/** Geometric quantile: smallest k with P(X <= k) >= q. */
export declare function geometricQuantile(q: number, p: number): number;
/** Natural log of the gamma function (Lanczos approximation). */
export declare function logGamma(z: number): number;
/** Student-t cumulative distribution function with df degrees of freedom. */
export declare function tCdf(t: number, df: number): number;
/** Student-t quantile function (bisection on tCdf); q in (0, 1), df > 0. */
export declare function tQuantile(q: number, df: number): number;
//# sourceMappingURL=math.d.ts.map