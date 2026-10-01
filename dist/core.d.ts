/**
 * dsh-probstat/core — engine functions behind the four tools.
 *
 * Each engine validates its inputs, computes with the math primitives and
 * assembles a result object. Discipline: only assign keys that hold values
 * (lossless-JSON gate), never `{ key: undefined }`.
 *
 * @module dsh-probstat/core
 */
export type DistKind = 'normal' | 'binomial' | 'poisson' | 'exponential' | 'uniform' | 'geometric';
export type DistOperation = 'pdf' | 'cdf' | 'survival' | 'quantile' | 'stats';
export interface DistCalcResult {
    valid: boolean;
    distribution?: string;
    operation?: string;
    result?: number;
    mean?: number;
    variance?: number;
    sd?: number;
    note?: string;
    error?: string;
}
export interface DistCalcArgs {
    distribution: string;
    operation?: string;
    x?: number;
    q?: number;
    n?: number;
    p?: number;
    lambda?: number;
    mean?: number;
    sd?: number;
    a?: number;
    b?: number;
}
/** Distribution calculator: pdf / cdf / survival / quantile / stats. */
export declare function distCalc(args: DistCalcArgs): DistCalcResult;
export type ZScoreMode = 'to_prob' | 'from_prob' | 'two_sided' | 'between';
export interface ZScoreResult {
    valid: boolean;
    mode?: string;
    z?: number;
    z1?: number;
    z2?: number;
    probability?: number;
    percentile?: number;
    pValue?: number;
    error?: string;
}
export interface ZScoreArgs {
    mode: string;
    z?: number;
    p?: number;
    z1?: number;
    z2?: number;
}
/** Standard-normal table math: z <-> probability conversions. */
export declare function zScore(args: ZScoreArgs): ZScoreResult;
export type CiKind = 'mean_z' | 'mean_t' | 'proportion';
export interface ConfidenceResult {
    valid: boolean;
    kind?: string;
    lower?: number;
    upper?: number;
    marginOfError?: number;
    criticalValue?: number;
    standardError?: number;
    center?: number;
    n?: number;
    conf?: number;
    error?: string;
}
export interface ConfidenceArgs {
    kind: string;
    mean?: number;
    sigma?: number;
    sd?: number;
    n?: number;
    successes?: number;
    conf?: number;
}
/** Confidence interval construction: mean (z / t) or proportion (Wilson). */
export declare function confidenceInterval(args: ConfidenceArgs): ConfidenceResult;
export type EventOperation = 'union' | 'intersection' | 'conditional' | 'bayes' | 'complement' | 'at_least_one';
export interface EventResult {
    valid: boolean;
    operation?: string;
    result?: number;
    formula?: string;
    error?: string;
}
export interface EventArgs {
    operation: string;
    p?: number;
    pA?: number;
    pB?: number;
    independent?: boolean;
    pIntersection?: number;
    pCondition?: number;
    prior?: number;
    truePositive?: number;
    falsePositive?: number;
    n?: number;
}
/** Event probability identities: union / intersection / conditional / Bayes. */
export declare function eventProbability(args: EventArgs): EventResult;
export type TestKind = 'mean_z' | 'mean_t' | 'proportion';
export type TestTail = 'two-sided' | 'left' | 'right';
export interface HypothesisResult {
    valid: boolean;
    kind?: string;
    tail?: string;
    statistic?: number;
    standardError?: number;
    df?: number;
    pValue?: number;
    pValueTwoSided?: number;
    pValueLeft?: number;
    pValueRight?: number;
    criticalValue?: number;
    alpha?: number;
    reject?: boolean;
    note?: string;
    error?: string;
}
export interface HypothesisArgs {
    kind: string;
    tail?: string;
    alpha?: number;
    sampleMean?: number;
    mu0?: number;
    sigma?: number;
    sd?: number;
    n?: number;
    successes?: number;
    p0?: number;
}
/**
 * One-sample hypothesis tests: mean with known sigma (z), mean with unknown
 * sigma (Student-t, df = n - 1) and a single proportion (normal approximation).
 * Statistic, all three p-values, the tail-specific critical value and the
 * alpha decision are reported together so no step is left to recall.
 */
export declare function hypothesisTest(args: HypothesisArgs): HypothesisResult;
export type SampleKind = 'mean_z' | 'proportion';
export type ProportionMethod = 'wald' | 'wilson';
export interface SampleSizeResult {
    valid: boolean;
    kind?: string;
    method?: string;
    n?: number;
    nExact?: number;
    criticalValue?: number;
    marginOfError?: number;
    conf?: number;
    assumedP?: number;
    note?: string;
    error?: string;
}
export interface SampleSizeArgs {
    kind: string;
    sigma?: number;
    p?: number;
    method?: string;
    marginOfError?: number;
    conf?: number;
}
/**
 * Required sample size for a target margin of error: mean with known sigma
 * (closed form) and proportion (Wald closed form or the Wilson-consistent
 * smallest n, solved by binary search on the Wilson half-width).
 */
export declare function sampleSize(args: SampleSizeArgs): SampleSizeResult;
//# sourceMappingURL=core.d.ts.map