/**
 * dsh-probstat/core — engine functions behind the four tools.
 *
 * Each engine validates its inputs, computes with the math primitives and
 * assembles a result object. Discipline: only assign keys that hold values
 * (lossless-JSON gate), never `{ key: undefined }`.
 *
 * @module dsh-probstat/core
 */

import {
  binomialCdf,
  binomialPmf,
  binomialQuantile,
  exponentialCdf,
  exponentialPdf,
  exponentialQuantile,
  exponentialSurvival,
  geometricCdf,
  geometricPmf,
  geometricQuantile,
  geometricSurvival,
  normalCdf,
  normalPdf,
  normalQuantile,
  poissonCdf,
  poissonPmf,
  poissonQuantile,
  tCdf,
  tQuantile,
  uniformCdf,
  uniformPdf,
  uniformQuantile,
} from './math.ts'

export type DistKind = 'normal' | 'binomial' | 'poisson' | 'exponential' | 'uniform' | 'geometric'
export type DistOperation = 'pdf' | 'cdf' | 'survival' | 'quantile' | 'stats'

export interface DistCalcResult {
  valid: boolean
  distribution?: string
  operation?: string
  result?: number
  mean?: number
  variance?: number
  sd?: number
  note?: string
  error?: string
}

export interface DistCalcArgs {
  distribution: string
  operation?: string
  x?: number
  q?: number
  n?: number
  p?: number
  lambda?: number
  mean?: number
  sd?: number
  a?: number
  b?: number
}

function isIntFinite(value: number | undefined, min: number): value is number {
  return value !== undefined && Number.isFinite(value) && Number.isInteger(value) && value >= min
}

function inRange(value: number | undefined, lo: number, hi: number): value is number {
  return value !== undefined && Number.isFinite(value) && value >= lo && value <= hi
}

/** Distribution calculator: pdf / cdf / survival / quantile / stats. */
export function distCalc(args: DistCalcArgs): DistCalcResult {
  const dist = args.distribution as DistKind
  const op = (args.operation ?? 'cdf') as DistOperation
  const KNOWN: readonly DistKind[] = ['normal', 'binomial', 'poisson', 'exponential', 'uniform', 'geometric']
  if (!KNOWN.includes(dist)) {
    return { valid: false, error: `unknown distribution "${args.distribution}" (expected one of normal, binomial, poisson, exponential, uniform, geometric)` }
  }
  if (!['pdf', 'cdf', 'survival', 'quantile', 'stats'].includes(op)) {
    return { valid: false, error: `unknown operation "${args.operation}" (expected pdf, cdf, survival, quantile or stats)` }
  }

  const needsX = op === 'pdf' || op === 'cdf' || op === 'survival'
  if (needsX && args.x === undefined) {
    return { valid: false, error: `operation "${op}" requires the x parameter` }
  }
  if (op === 'quantile' && args.q === undefined) {
    return { valid: false, error: 'operation "quantile" requires the q parameter (probability in (0, 1))' }
  }
  if (op === 'quantile' && !inRange(args.q, 0, 1)) {
    return { valid: false, error: `q must lie strictly in (0, 1) (got ${args.q})` }
  }

  const result: DistCalcResult = { valid: true, distribution: dist, operation: op }

  switch (dist) {
    case 'normal': {
      const mean = args.mean ?? 0
      const sd = args.sd ?? 1
      if (!Number.isFinite(mean)) return { valid: false, error: `mean must be a finite number (got ${args.mean})` }
      if (sd === undefined || !Number.isFinite(sd) || sd <= 0) {
        return { valid: false, error: `sd must be a positive number (got ${args.sd})` }
      }
      if (op === 'stats') {
        result.mean = mean
        result.variance = sd * sd
        result.sd = sd
        return result
      }
      if (op === 'pdf') result.result = normalPdf(args.x!, mean, sd)
      else if (op === 'cdf') result.result = normalCdf(args.x!, mean, sd)
      else if (op === 'survival') result.result = 1 - normalCdf(args.x!, mean, sd)
      else result.result = normalQuantile(args.q!, mean, sd)
      return result
    }
    case 'binomial': {
      if (!isIntFinite(args.n, 1)) return { valid: false, error: `n must be an integer >= 1 (got ${args.n})` }
      if (!inRange(args.p, 0, 1)) return { valid: false, error: `p must lie in [0, 1] (got ${args.p})` }
      const n = args.n
      const p = args.p
      if (op === 'stats') {
        result.mean = n * p
        result.variance = n * p * (1 - p)
        result.sd = Math.sqrt(result.variance)
        return result
      }
      if (op === 'quantile') {
        result.result = binomialQuantile(args.q!, n, p)
        result.note = 'quantile = smallest k with P(X <= k) >= q'
        return result
      }
      if (!Number.isInteger(args.x)) {
        return { valid: false, error: `x must be an integer for a binomial distribution (got ${args.x})` }
      }
      if (args.x! < 0 || args.x! > n) {
        return { valid: false, error: `x must lie in [0, n] = [0, ${n}] (got ${args.x})` }
      }
      if (op === 'pdf') result.result = binomialPmf(args.x!, n, p)
      else if (op === 'cdf') result.result = binomialCdf(args.x!, n, p)
      else result.result = 1 - binomialCdf(args.x!, n, p)
      return result
    }
    case 'poisson': {
      if (args.lambda === undefined || !Number.isFinite(args.lambda) || args.lambda < 0) {
        return { valid: false, error: `lambda must be a non-negative number (got ${args.lambda})` }
      }
      const lambda = args.lambda
      if (op === 'stats') {
        result.mean = lambda
        result.variance = lambda
        result.sd = Math.sqrt(lambda)
        return result
      }
      if (op === 'quantile') {
        result.result = poissonQuantile(args.q!, lambda)
        result.note = 'quantile = smallest k with P(X <= k) >= q'
        return result
      }
      if (!Number.isInteger(args.x)) {
        return { valid: false, error: `x must be an integer for a poisson distribution (got ${args.x})` }
      }
      if (args.x! < 0) return { valid: false, error: `x must be >= 0 (got ${args.x})` }
      if (op === 'pdf') result.result = poissonPmf(args.x!, lambda)
      else if (op === 'cdf') result.result = poissonCdf(args.x!, lambda)
      else result.result = 1 - poissonCdf(args.x!, lambda)
      return result
    }
    case 'exponential': {
      if (args.lambda === undefined || !Number.isFinite(args.lambda) || args.lambda <= 0) {
        return { valid: false, error: `lambda (rate) must be a positive number (got ${args.lambda})` }
      }
      const lambda = args.lambda
      if (op === 'stats') {
        result.mean = 1 / lambda
        result.variance = 1 / (lambda * lambda)
        result.sd = 1 / lambda
        return result
      }
      if (op === 'quantile') {
        result.result = exponentialQuantile(args.q!, lambda)
        return result
      }
      if (op === 'pdf') result.result = exponentialPdf(args.x!, lambda)
      else if (op === 'cdf') result.result = exponentialCdf(args.x!, lambda)
      else result.result = exponentialSurvival(args.x!, lambda)
      if (args.x! < 0) result.note = 'x < 0: density 0, P(X <= x) = 0, P(X > x) = 1'
      return result
    }
    case 'uniform': {
      if (args.a === undefined || !Number.isFinite(args.a)) return { valid: false, error: `a must be a finite number (got ${args.a})` }
      if (args.b === undefined || !Number.isFinite(args.b)) return { valid: false, error: `b must be a finite number (got ${args.b})` }
      if (args.a >= args.b) return { valid: false, error: `a must be less than b (got a=${args.a}, b=${args.b})` }
      const a = args.a
      const b = args.b
      if (op === 'stats') {
        result.mean = (a + b) / 2
        result.variance = ((b - a) * (b - a)) / 12
        result.sd = (b - a) / Math.sqrt(12)
        return result
      }
      if (op === 'quantile') {
        result.result = uniformQuantile(args.q!, a, b)
        return result
      }
      if (op === 'pdf') result.result = uniformPdf(args.x!, a, b)
      else if (op === 'cdf') result.result = uniformCdf(args.x!, a, b)
      else result.result = 1 - uniformCdf(args.x!, a, b)
      return result
    }
    case 'geometric': {
      if (!inRange(args.p, 0, 1) || args.p === 0) {
        return { valid: false, error: `p must lie in (0, 1] (got ${args.p})` }
      }
      const p = args.p
      result.note = 'support convention: k = number of failures before the first success (k >= 0)'
      if (op === 'stats') {
        result.mean = (1 - p) / p
        result.variance = (1 - p) / (p * p)
        result.sd = Math.sqrt(1 - p) / p
        return result
      }
      if (op === 'quantile') {
        result.result = geometricQuantile(args.q!, p)
        result.note = 'quantile = smallest k with P(X <= k) >= q; failures-before-success convention'
        return result
      }
      if (!Number.isInteger(args.x)) {
        return { valid: false, error: `x must be an integer for a geometric distribution (got ${args.x})` }
      }
      if (args.x! < 0) return { valid: false, error: `x must be >= 0 (got ${args.x})` }
      if (op === 'pdf') result.result = geometricPmf(args.x!, p)
      else if (op === 'cdf') result.result = geometricCdf(args.x!, p)
      else result.result = geometricSurvival(args.x!, p)
      return result
    }
  }
}

/* ------------------------------------------------------------------ */
/* z_score                                                             */
/* ------------------------------------------------------------------ */

export type ZScoreMode = 'to_prob' | 'from_prob' | 'two_sided' | 'between'

export interface ZScoreResult {
  valid: boolean
  mode?: string
  z?: number
  z1?: number
  z2?: number
  probability?: number
  percentile?: number
  pValue?: number
  error?: string
}

export interface ZScoreArgs {
  mode: string
  z?: number
  p?: number
  z1?: number
  z2?: number
}

/** Standard-normal table math: z <-> probability conversions. */
export function zScore(args: ZScoreArgs): ZScoreResult {
  const mode = args.mode
  switch (mode) {
    case 'to_prob': {
      if (args.z === undefined || !Number.isFinite(args.z)) return { valid: false, error: `to_prob requires a finite z (got ${args.z})` }
      const prob = normalCdf(args.z, 0, 1)
      return { valid: true, mode, z: args.z, probability: prob, percentile: prob * 100 }
    }
    case 'from_prob': {
      if (!inRange(args.p, 0, 1) || args.p === 0 || args.p === 1) {
        return { valid: false, error: `from_prob requires p strictly in (0, 1) (got ${args.p})` }
      }
      return { valid: true, mode, z: normalQuantile(args.p, 0, 1), probability: args.p, percentile: args.p * 100 }
    }
    case 'two_sided': {
      if (args.z === undefined || !Number.isFinite(args.z)) return { valid: false, error: `two_sided requires a finite z (got ${args.z})` }
      const az = Math.abs(args.z)
      return { valid: true, mode, z: args.z, pValue: 2 * (1 - normalCdf(az, 0, 1)) }
    }
    case 'between': {
      if (args.z1 === undefined || !Number.isFinite(args.z1)) return { valid: false, error: `between requires a finite z1 (got ${args.z1})` }
      if (args.z2 === undefined || !Number.isFinite(args.z2)) return { valid: false, error: `between requires a finite z2 (got ${args.z2})` }
      const lo = Math.min(args.z1, args.z2)
      const hi = Math.max(args.z1, args.z2)
      return { valid: true, mode, z1: lo, z2: hi, probability: normalCdf(hi, 0, 1) - normalCdf(lo, 0, 1) }
    }
    default:
      return { valid: false, error: `unknown mode "${args.mode}" (expected to_prob, from_prob, two_sided or between)` }
  }
}

/* ------------------------------------------------------------------ */
/* confidence_interval                                                 */
/* ------------------------------------------------------------------ */

export type CiKind = 'mean_z' | 'mean_t' | 'proportion'

export interface ConfidenceResult {
  valid: boolean
  kind?: string
  lower?: number
  upper?: number
  marginOfError?: number
  criticalValue?: number
  standardError?: number
  center?: number
  n?: number
  conf?: number
  error?: string
}

export interface ConfidenceArgs {
  kind: string
  mean?: number
  sigma?: number
  sd?: number
  n?: number
  successes?: number
  conf?: number
}

/** Confidence interval construction: mean (z / t) or proportion (Wilson). */
export function confidenceInterval(args: ConfidenceArgs): ConfidenceResult {
  const kind = args.kind
  const conf = args.conf ?? 0.95
  if (!inRange(conf, 0, 1) || conf === 0 || conf === 1) {
    return { valid: false, error: `conf must lie strictly in (0, 1) (got ${args.conf})` }
  }
  const tail = (1 - conf) / 2

  if (kind === 'mean_z') {
    if (args.mean === undefined || !Number.isFinite(args.mean)) return { valid: false, error: `mean is required (got ${args.mean})` }
    if (args.sigma === undefined || !Number.isFinite(args.sigma) || args.sigma <= 0) {
      return { valid: false, error: `sigma must be a positive number (got ${args.sigma})` }
    }
    if (!isIntFinite(args.n, 1)) return { valid: false, error: `n must be an integer >= 1 (got ${args.n})` }
    const se = args.sigma / Math.sqrt(args.n)
    const z = normalQuantile(1 - tail, 0, 1)
    return {
      valid: true,
      kind,
      center: args.mean,
      lower: args.mean - z * se,
      upper: args.mean + z * se,
      marginOfError: z * se,
      criticalValue: z,
      standardError: se,
      n: args.n,
      conf,
    }
  }

  if (kind === 'mean_t') {
    if (args.mean === undefined || !Number.isFinite(args.mean)) return { valid: false, error: `mean is required (got ${args.mean})` }
    if (args.sd === undefined || !Number.isFinite(args.sd) || args.sd <= 0) {
      return { valid: false, error: `sd must be a positive number (got ${args.sd})` }
    }
    if (!isIntFinite(args.n, 2)) return { valid: false, error: `n must be an integer >= 2 (got ${args.n})` }
    const se = args.sd / Math.sqrt(args.n)
    const t = tQuantile(1 - tail, args.n - 1)
    return {
      valid: true,
      kind,
      center: args.mean,
      lower: args.mean - t * se,
      upper: args.mean + t * se,
      marginOfError: t * se,
      criticalValue: t,
      standardError: se,
      n: args.n,
      conf,
    }
  }

  if (kind === 'proportion') {
    if (!isIntFinite(args.n, 1)) return { valid: false, error: `n must be an integer >= 1 (got ${args.n})` }
    if (!isIntFinite(args.successes, 0) || args.successes > args.n) {
      return { valid: false, error: `successes must be an integer in [0, n] = [0, ${args.n}] (got ${args.successes})` }
    }
    const n = args.n
    const k = args.successes
    const z = normalQuantile(1 - tail, 0, 1)
    const phat = k / n
    const denom = 1 + (z * z) / n
    const center = (phat + (z * z) / (2 * n)) / denom
    const half = (z * Math.sqrt(phat * (1 - phat) / n + (z * z) / (4 * n * n))) / denom
    return {
      valid: true,
      kind,
      center,
      lower: center - half,
      upper: center + half,
      marginOfError: half,
      criticalValue: z,
      standardError: Math.sqrt(phat * (1 - phat) / n),
      n,
      conf,
    }
  }

  return { valid: false, error: `unknown kind "${args.kind}" (expected mean_z, mean_t or proportion)` }
}

/* ------------------------------------------------------------------ */
/* event_probability                                                   */
/* ------------------------------------------------------------------ */

export type EventOperation = 'union' | 'intersection' | 'conditional' | 'bayes' | 'complement' | 'at_least_one'

export interface EventResult {
  valid: boolean
  operation?: string
  result?: number
  formula?: string
  error?: string
}

export interface EventArgs {
  operation: string
  p?: number
  pA?: number
  pB?: number
  independent?: boolean
  pIntersection?: number
  pCondition?: number
  prior?: number
  truePositive?: number
  falsePositive?: number
  n?: number
}

/** Event probability identities: union / intersection / conditional / Bayes. */
export function eventProbability(args: EventArgs): EventResult {
  switch (args.operation) {
    case 'union': {
      if (!inRange(args.pA, 0, 1)) return { valid: false, error: `pA must lie in [0, 1] (got ${args.pA})` }
      if (!inRange(args.pB, 0, 1)) return { valid: false, error: `pB must lie in [0, 1] (got ${args.pB})` }
      const independent = args.independent ?? true
      if (independent) {
        return {
          valid: true,
          operation: args.operation,
          result: args.pA + args.pB - args.pA * args.pB,
          formula: 'P(A or B) = P(A) + P(B) - P(A)P(B) (independent events)',
        }
      }
      return {
        valid: true,
        operation: args.operation,
        result: args.pA + args.pB,
        formula: 'P(A or B) = P(A) + P(B) (treated as mutually exclusive when independent=false)',
      }
    }
    case 'intersection': {
      if (!inRange(args.pA, 0, 1)) return { valid: false, error: `pA must lie in [0, 1] (got ${args.pA})` }
      if (!inRange(args.pB, 0, 1)) return { valid: false, error: `pB must lie in [0, 1] (got ${args.pB})` }
      const independent = args.independent ?? true
      if (!independent) {
        return { valid: false, error: 'P(A and B) cannot be computed from marginal probabilities alone; use the conditional operation with P(A and B) given directly' }
      }
      return {
        valid: true,
        operation: args.operation,
        result: args.pA * args.pB,
        formula: 'P(A and B) = P(A)P(B) (independent events)',
      }
    }
    case 'conditional': {
      if (!inRange(args.pIntersection, 0, 1)) return { valid: false, error: `pIntersection must lie in [0, 1] (got ${args.pIntersection})` }
      if (!inRange(args.pCondition, 0, 1) || args.pCondition === 0) {
        return { valid: false, error: `pCondition must lie in (0, 1] (got ${args.pCondition})` }
      }
      if (args.pIntersection > args.pCondition) {
        return { valid: false, error: `pIntersection cannot exceed pCondition (got ${args.pIntersection} > ${args.pCondition})` }
      }
      return {
        valid: true,
        operation: args.operation,
        result: args.pIntersection / args.pCondition,
        formula: 'P(A|B) = P(A and B) / P(B)',
      }
    }
    case 'bayes': {
      if (!inRange(args.prior, 0, 1)) return { valid: false, error: `prior must lie in [0, 1] (got ${args.prior})` }
      if (!inRange(args.truePositive, 0, 1)) return { valid: false, error: `truePositive (P(E|H)) must lie in [0, 1] (got ${args.truePositive})` }
      if (!inRange(args.falsePositive, 0, 1)) return { valid: false, error: `falsePositive (P(E|not H)) must lie in [0, 1] (got ${args.falsePositive})` }
      const denom = args.prior * args.truePositive + (1 - args.prior) * args.falsePositive
      if (denom === 0) return { valid: false, error: 'P(E) evaluates to 0; posterior is undefined' }
      return {
        valid: true,
        operation: args.operation,
        result: (args.prior * args.truePositive) / denom,
        formula: 'P(H|E) = P(H)P(E|H) / [P(H)P(E|H) + P(not H)P(E|not H)]',
      }
    }
    case 'complement': {
      if (!inRange(args.p, 0, 1)) return { valid: false, error: `p must lie in [0, 1] (got ${args.p})` }
      return { valid: true, operation: args.operation, result: 1 - args.p, formula: 'P(not A) = 1 - P(A)' }
    }
    case 'at_least_one': {
      if (!inRange(args.p, 0, 1)) return { valid: false, error: `p must lie in [0, 1] (got ${args.p})` }
      if (!isIntFinite(args.n, 1)) return { valid: false, error: `n must be an integer >= 1 (got ${args.n})` }
      return {
        valid: true,
        operation: args.operation,
        result: 1 - Math.pow(1 - args.p, args.n),
        formula: 'P(at least one success in n trials) = 1 - (1 - p)^n',
      }
    }
    default:
      return { valid: false, error: `unknown operation "${args.operation}" (expected union, intersection, conditional, bayes, complement or at_least_one)` }
  }
}

/* ------------------------------------------------------------------ */
/* hypothesis_test                                                     */
/* ------------------------------------------------------------------ */

export type TestKind = 'mean_z' | 'mean_t' | 'proportion'
export type TestTail = 'two-sided' | 'left' | 'right'

const TAILS: readonly TestTail[] = ['two-sided', 'left', 'right']

export interface HypothesisResult {
  valid: boolean
  kind?: string
  tail?: string
  statistic?: number
  standardError?: number
  df?: number
  pValue?: number
  pValueTwoSided?: number
  pValueLeft?: number
  pValueRight?: number
  criticalValue?: number
  alpha?: number
  reject?: boolean
  note?: string
  error?: string
}

export interface HypothesisArgs {
  kind: string
  tail?: string
  alpha?: number
  sampleMean?: number
  mu0?: number
  sigma?: number
  sd?: number
  n?: number
  successes?: number
  p0?: number
}

function pValueForTail(tail: TestTail, two: number, left: number, right: number): number {
  if (tail === 'left') return left
  if (tail === 'right') return right
  return two
}

function criticalForTail(tail: TestTail, alpha: number, quantile: (p: number) => number): number {
  if (tail === 'left') return quantile(alpha)
  if (tail === 'right') return quantile(1 - alpha)
  return quantile(1 - alpha / 2)
}

/**
 * One-sample hypothesis tests: mean with known sigma (z), mean with unknown
 * sigma (Student-t, df = n - 1) and a single proportion (normal approximation).
 * Statistic, all three p-values, the tail-specific critical value and the
 * alpha decision are reported together so no step is left to recall.
 */
export function hypothesisTest(args: HypothesisArgs): HypothesisResult {
  const kind = args.kind
  const alpha = args.alpha ?? 0.05
  if (!inRange(alpha, 0, 1) || alpha === 0 || alpha === 1) {
    return { valid: false, error: `alpha must lie strictly in (0, 1) (got ${args.alpha})` }
  }
  const tail = (args.tail ?? 'two-sided') as TestTail
  if (!TAILS.includes(tail)) {
    return { valid: false, error: `unknown tail "${args.tail}" (expected two-sided, left or right)` }
  }

  /** Assemble the shared result: p-values from the null cdf, critical value per tail. */
  const finish = (
    statistic: number,
    standardError: number,
    df: number | undefined,
    quantile: (p: number) => number,
    note?: string,
  ): HypothesisResult => {
    const cdf = df === undefined ? (t: number) => normalCdf(t, 0, 1) : (t: number) => tCdf(t, df)
    const left = cdf(statistic)
    const right = 1 - left
    const two = 2 * (1 - cdf(Math.abs(statistic)))
    const pValue = pValueForTail(tail, two, left, right)
    const result: HypothesisResult = {
      valid: true,
      kind,
      tail,
      statistic,
      standardError,
      pValue,
      pValueTwoSided: two,
      pValueLeft: left,
      pValueRight: right,
      criticalValue: criticalForTail(tail, alpha, quantile),
      alpha,
      reject: pValue < alpha,
    }
    if (df !== undefined) result.df = df
    if (note !== undefined) result.note = note
    return result
  }

  if (kind === 'mean_z') {
    if (args.sampleMean === undefined || !Number.isFinite(args.sampleMean)) {
      return { valid: false, error: `sampleMean is required (got ${args.sampleMean})` }
    }
    if (args.mu0 === undefined || !Number.isFinite(args.mu0)) {
      return { valid: false, error: `mu0 (null-hypothesis mean) is required (got ${args.mu0})` }
    }
    if (args.sigma === undefined || !Number.isFinite(args.sigma) || args.sigma <= 0) {
      return { valid: false, error: `sigma (known population sd) must be positive (got ${args.sigma})` }
    }
    if (!isIntFinite(args.n, 1)) return { valid: false, error: `n must be an integer >= 1 (got ${args.n})` }
    const se = args.sigma / Math.sqrt(args.n)
    return finish((args.sampleMean - args.mu0) / se, se, undefined, (p) => normalQuantile(p, 0, 1))
  }

  if (kind === 'mean_t') {
    if (args.sampleMean === undefined || !Number.isFinite(args.sampleMean)) {
      return { valid: false, error: `sampleMean is required (got ${args.sampleMean})` }
    }
    if (args.mu0 === undefined || !Number.isFinite(args.mu0)) {
      return { valid: false, error: `mu0 (null-hypothesis mean) is required (got ${args.mu0})` }
    }
    if (args.sd === undefined || !Number.isFinite(args.sd) || args.sd <= 0) {
      return { valid: false, error: `sd (sample sd) must be positive (got ${args.sd})` }
    }
    if (!isIntFinite(args.n, 2)) return { valid: false, error: `n must be an integer >= 2 (got ${args.n})` }
    const df = args.n - 1
    const se = args.sd / Math.sqrt(args.n)
    return finish((args.sampleMean - args.mu0) / se, se, df, (p) => tQuantile(p, df))
  }

  if (kind === 'proportion') {
    if (!isIntFinite(args.n, 1)) return { valid: false, error: `n must be an integer >= 1 (got ${args.n})` }
    if (!isIntFinite(args.successes, 0) || args.successes > args.n) {
      return { valid: false, error: `successes must be an integer in [0, n] = [0, ${args.n}] (got ${args.successes})` }
    }
    if (!inRange(args.p0, 0, 1) || args.p0 === 0 || args.p0 === 1) {
      return { valid: false, error: `p0 (null-hypothesis proportion) must lie strictly in (0, 1) (got ${args.p0})` }
    }
    const n = args.n
    const p0 = args.p0
    const se = Math.sqrt((p0 * (1 - p0)) / n)
    const statistic = (args.successes / n - p0) / se
    const expected = Math.min(n * p0, n * (1 - p0))
    const note = expected < 5
      ? `normal approximation is shaky here: n*p0 = ${Number(expected.toPrecision(6))} (< 5); treat the p-value as indicative only`
      : undefined
    return finish(statistic, se, undefined, (p) => normalQuantile(p, 0, 1), note)
  }

  return { valid: false, error: `unknown kind "${args.kind}" (expected mean_z, mean_t or proportion)` }
}

/* ------------------------------------------------------------------ */
/* sample_size                                                         */
/* ------------------------------------------------------------------ */

export type SampleKind = 'mean_z' | 'proportion'
export type ProportionMethod = 'wald' | 'wilson'

/** Upper bound of the Wilson sample-size search (guards against runaway loops). */
const SAMPLE_SIZE_MAX_N = 1e7

export interface SampleSizeResult {
  valid: boolean
  kind?: string
  method?: string
  n?: number
  nExact?: number
  criticalValue?: number
  marginOfError?: number
  conf?: number
  assumedP?: number
  note?: string
  error?: string
}

export interface SampleSizeArgs {
  kind: string
  sigma?: number
  p?: number
  method?: string
  marginOfError?: number
  conf?: number
}

/**
 * Required sample size for a target margin of error: mean with known sigma
 * (closed form) and proportion (Wald closed form or the Wilson-consistent
 * smallest n, solved by binary search on the Wilson half-width).
 */
export function sampleSize(args: SampleSizeArgs): SampleSizeResult {
  const kind = args.kind
  const conf = args.conf ?? 0.95
  if (!inRange(conf, 0, 1) || conf === 0 || conf === 1) {
    return { valid: false, error: `conf must lie strictly in (0, 1) (got ${args.conf})` }
  }
  if (args.marginOfError === undefined || !Number.isFinite(args.marginOfError) || args.marginOfError <= 0) {
    return { valid: false, error: `marginOfError must be a positive number (got ${args.marginOfError})` }
  }
  const e = args.marginOfError
  const z = normalQuantile(1 - (1 - conf) / 2, 0, 1)

  if (kind === 'mean_z') {
    if (args.sigma === undefined || !Number.isFinite(args.sigma) || args.sigma <= 0) {
      return { valid: false, error: `sigma (known population sd) must be positive (got ${args.sigma})` }
    }
    const nExact = ((z * args.sigma) / e) ** 2
    return { valid: true, kind, criticalValue: z, marginOfError: e, conf, nExact, n: Math.ceil(nExact) }
  }

  if (kind === 'proportion') {
    const method = (args.method ?? 'wald') as ProportionMethod
    if (method !== 'wald' && method !== 'wilson') {
      return { valid: false, error: `unknown method "${args.method}" (expected wald or wilson)` }
    }
    const p = args.p ?? 0.5
    if (!Number.isFinite(p) || p <= 0 || p >= 1) {
      return { valid: false, error: `p (expected proportion) must lie strictly in (0, 1) (got ${args.p})` }
    }
    const assumed = args.p === undefined
      ? 'p defaults to 0.5 (the conservative worst case: p(1-p) is maximal)'
      : undefined

    if (method === 'wald') {
      const nExact = (z * z * p * (1 - p)) / (e * e)
      const result: SampleSizeResult = {
        valid: true,
        kind,
        method,
        criticalValue: z,
        marginOfError: e,
        conf,
        assumedP: p,
        nExact,
        n: Math.ceil(nExact),
      }
      result.note = assumed ?? 'Wald closed form n = z^2 p(1-p) / E^2'
      return result
    }

    // Wilson: half-width is strictly decreasing in n for fixed p and conf,
    // so the smallest n meeting the target is found by binary search.
    const halfWidth = (n: number): number => {
      const w = (z * Math.sqrt((p * (1 - p)) / n + (z * z) / (4 * n * n))) / (1 + (z * z) / n)
      return w
    }
    if (halfWidth(SAMPLE_SIZE_MAX_N) > e) {
      return { valid: false, error: `no sample size up to ${SAMPLE_SIZE_MAX_N} reaches margin of error ${e}; loosen marginOfError` }
    }
    let lo = 1
    let hi = SAMPLE_SIZE_MAX_N
    while (lo < hi) {
      const mid = Math.floor((lo + hi) / 2)
      if (halfWidth(mid) <= e) hi = mid
      else lo = mid + 1
    }
    const result: SampleSizeResult = {
      valid: true,
      kind,
      method,
      criticalValue: z,
      marginOfError: e,
      conf,
      assumedP: p,
      n: lo,
    }
    result.note = `Wilson target: smallest n whose Wilson half-width (${fmtNumber(halfWidth(lo))}) is <= ${fmtNumber(e)}${assumed === undefined ? '' : `; ${assumed}`}`
    return result
  }

  return { valid: false, error: `unknown kind "${args.kind}" (expected mean_z or proportion)` }
}

/** Compact numeric rendering for note strings (never exponent-soup). */
function fmtNumber(value: number): string {
  return String(Number(value.toPrecision(8)))
}
