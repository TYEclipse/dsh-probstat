/**
 * Tool definitions for dsh-probstat: four deterministic probability and
 * statistical-inference tools exposed to every agent via defineTool.
 *
 * @module dsh-probstat/tools
 */

import { defineTool, type ToolDefinition } from '@deepseek-ai/dsh-tools'
import {
  confidenceInterval,
  distCalc,
  eventProbability,
  hypothesisTest,
  sampleSize,
  zScore,
  type ConfidenceResult,
  type DistCalcResult,
  type EventResult,
  type HypothesisResult,
  type SampleSizeResult,
  type ZScoreResult,
} from './core.ts'

export interface ToolSet {
  dist_calc: ToolDefinition
  z_score: ToolDefinition
  confidence_interval: ToolDefinition
  event_probability: ToolDefinition
  hypothesis_test: ToolDefinition
  sample_size: ToolDefinition
}

function fmt(value: number): string {
  if (!Number.isFinite(value)) return String(value)
  if (value !== 0 && (Math.abs(value) >= 1e7 || Math.abs(value) < 1e-4)) {
    return value.toExponential(6)
  }
  return String(Number(value.toPrecision(8)))
}

/* ------------------------------------------------------------------ */
/* render helpers                                                      */
/* ------------------------------------------------------------------ */

function renderDistCalc(value: unknown): string {
  const v = value as DistCalcResult
  if (!v.valid) return `dist_calc failed: ${v.error}`
  const lines: string[] = [`${v.distribution} distribution (operation: ${v.operation})`]
  if (v.result !== undefined) lines.push(`result: ${fmt(v.result)}`)
  if (v.mean !== undefined) lines.push(`mean: ${fmt(v.mean)}`)
  if (v.variance !== undefined) lines.push(`variance: ${fmt(v.variance)}`)
  if (v.sd !== undefined) lines.push(`sd: ${fmt(v.sd)}`)
  if (v.note !== undefined) lines.push(`note: ${v.note}`)
  return lines.join('\n')
}

function renderZScore(value: unknown): string {
  const v = value as ZScoreResult
  if (!v.valid) return `z_score failed: ${v.error}`
  const lines: string[] = [`mode: ${v.mode}`]
  if (v.z !== undefined) lines.push(`z: ${fmt(v.z)}`)
  if (v.z1 !== undefined && v.z2 !== undefined) lines.push(`z range: [${fmt(v.z1)}, ${fmt(v.z2)}]`)
  if (v.probability !== undefined) lines.push(`probability: ${fmt(v.probability)}`)
  if (v.percentile !== undefined) lines.push(`percentile: ${fmt(v.percentile)} %`)
  if (v.pValue !== undefined) lines.push(`two-sided p-value: ${fmt(v.pValue)}`)
  return lines.join('\n')
}

function renderConfidence(value: unknown): string {
  const v = value as ConfidenceResult
  if (!v.valid) return `confidence_interval failed: ${v.error}`
  const kind = v.kind === 'proportion' ? 'proportion (Wilson)' : v.kind
  const lines: string[] = [`${kind} ${fmt((v.conf ?? 1) * 100)}% confidence interval`]
  lines.push(`[${fmt(v.lower!)} , ${fmt(v.upper!)}]`)
  if (v.center !== undefined) lines.push(`center: ${fmt(v.center)}`)
  lines.push(`margin of error: ${fmt(v.marginOfError!)}`)
  lines.push(`critical value: ${fmt(v.criticalValue!)}`)
  lines.push(`standard error: ${fmt(v.standardError!)}`)
  return lines.join('\n')
}

function renderEvent(value: unknown): string {
  const v = value as EventResult
  if (!v.valid) return `event_probability failed: ${v.error}`
  const lines: string[] = [`${v.operation}: ${fmt(v.result!)}`]
  if (v.formula !== undefined) lines.push(`formula: ${v.formula}`)
  return lines.join('\n')
}

function renderHypothesis(value: unknown): string {
  const v = value as HypothesisResult
  if (!v.valid) return `hypothesis_test failed: ${v.error}`
  const lines: string[] = [`${v.kind} test — ${v.tail} tail, alpha = ${fmt(v.alpha!)}`]
  lines.push(`statistic: ${fmt(v.statistic!)}`)
  if (v.df !== undefined) lines.push(`df: ${fmt(v.df)}`)
  lines.push(`standard error: ${fmt(v.standardError!)}`)
  lines.push(`p-value (${v.tail}): ${fmt(v.pValue!)}`)
  lines.push(`p-values: two-sided ${fmt(v.pValueTwoSided!)} | left ${fmt(v.pValueLeft!)} | right ${fmt(v.pValueRight!)}`)
  lines.push(`critical value (${v.tail}): ${fmt(v.criticalValue!)}`)
  lines.push(v.reject
    ? `decision: reject the null hypothesis at alpha = ${fmt(v.alpha!)}`
    : `decision: do not reject the null hypothesis at alpha = ${fmt(v.alpha!)}`)
  if (v.note !== undefined) lines.push(`note: ${v.note}`)
  return lines.join('\n')
}

function renderSampleSize(value: unknown): string {
  const v = value as SampleSizeResult
  if (!v.valid) return `sample_size failed: ${v.error}`
  const lines: string[] = [
    `${v.kind} sample size — ${fmt((v.conf ?? 1) * 100)}% confidence, margin of error ${fmt(v.marginOfError!)}`,
  ]
  if (v.method !== undefined) lines.push(`method: ${v.method}`)
  lines.push(`required n: ${v.n}`)
  if (v.nExact !== undefined) lines.push(`exact (uncorrected) n: ${fmt(v.nExact)}`)
  if (v.assumedP !== undefined) lines.push(`assumed proportion p: ${fmt(v.assumedP)}`)
  lines.push(`critical value: ${fmt(v.criticalValue!)}`)
  if (v.note !== undefined) lines.push(`note: ${v.note}`)
  return lines.join('\n')
}

/* ------------------------------------------------------------------ */
/* tool registration                                                   */
/* ------------------------------------------------------------------ */

export function buildProbstatTools(): ToolSet {
  const dist_calc = defineTool({
    name: 'dist_calc',
    description: 'Deterministic probability distribution math. Six distributions: normal (mean/sd, default '
      + 'standard normal), binomial (n, p), poisson (lambda), exponential (rate lambda), uniform (a, b) and '
      + 'geometric (p; support counts failures before the first success). Five operations: pdf (density/mass '
      + 'at x), cdf (P(X <= x)), survival (P(X > x)), quantile (inverse CDF at probability q) and stats '
      + '(mean/variance/sd). Use this instead of mental arithmetic for distribution questions.',
    parameters: {
      distribution: {
        type: 'string',
        required: true,
        enum: ['normal', 'binomial', 'poisson', 'exponential', 'uniform', 'geometric'],
        description: 'Which distribution to evaluate.',
      },
      operation: {
        type: 'string',
        enum: ['pdf', 'cdf', 'survival', 'quantile', 'stats'],
        description: 'Which quantity to compute (default: cdf).',
      },
      x: { type: 'number', description: 'Point to evaluate (required for pdf / cdf / survival).' },
      q: { type: 'number', description: 'Probability in (0, 1) for quantile lookups.' },
      n: { type: 'number', description: 'Number of trials (binomial; integer >= 1).' },
      p: { type: 'number', description: 'Success probability per trial (binomial: [0, 1]; geometric: (0, 1]).' },
      lambda: { type: 'number', description: 'Rate parameter (poisson: >= 0; exponential: > 0).' },
      mean: { type: 'number', description: 'Normal mean (default 0).' },
      sd: { type: 'number', description: 'Normal standard deviation, > 0 (default 1).' },
      a: { type: 'number', description: 'Uniform lower bound (must be less than b).' },
      b: { type: 'number', description: 'Uniform upper bound.' },
    },
    output: {
      schema: {
        type: 'object',
        additionalProperties: false,
        properties: {
          valid: { type: 'boolean', required: true },
          distribution: { type: 'string' },
          operation: { type: 'string' },
          result: { type: 'number' },
          mean: { type: 'number' },
          variance: { type: 'number' },
          sd: { type: 'number' },
          note: { type: 'string' },
          error: { type: 'string' },
        },
      },
      render: (_args: unknown, value: unknown) => [{ type: 'text', text: renderDistCalc(value) }],
    },
    async execute(args: {
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
    }): Promise<DistCalcResult> {
      return distCalc(args)
    },
  })

  const z_score = defineTool({
    name: 'z_score',
    description: 'Standard-normal (z-table) math. Four modes: "to_prob" (z -> P(Z <= z) and percentile), '
      + '"from_prob" (probability p in (0, 1) -> z), "two_sided" (z -> two-tailed p-value 2*(1-Phi(|z|))) '
      + 'and "between" (z1, z2 -> P(z1 < Z < z2); order does not matter). Use this instead of recalling '
      + 'z-table values from memory.',
    parameters: {
      mode: { type: 'string', required: true, enum: ['to_prob', 'from_prob', 'two_sided', 'between'], description: 'Which z-table conversion to perform.' },
      z: { type: 'number', description: 'Standard score (to_prob / two_sided modes).' },
      p: { type: 'number', description: 'Probability in (0, 1) (from_prob mode).' },
      z1: { type: 'number', description: 'Lower bound (between mode; order is sorted internally).' },
      z2: { type: 'number', description: 'Upper bound (between mode; order is sorted internally).' },
    },
    output: {
      schema: {
        type: 'object',
        additionalProperties: false,
        properties: {
          valid: { type: 'boolean', required: true },
          mode: { type: 'string' },
          z: { type: 'number' },
          z1: { type: 'number' },
          z2: { type: 'number' },
          probability: { type: 'number' },
          percentile: { type: 'number' },
          pValue: { type: 'number' },
          error: { type: 'string' },
        },
      },
      render: (_args: unknown, value: unknown) => [{ type: 'text', text: renderZScore(value) }],
    },
    async execute(args: { mode: string; z?: number; p?: number; z1?: number; z2?: number }): Promise<ZScoreResult> {
      return zScore(args)
    },
  })

  const confidence_interval = defineTool({
    name: 'confidence_interval',
    description: 'Confidence interval construction. Three kinds: "mean_z" (population mean, known sigma, '
      + 'z critical values), "mean_t" (population mean, unknown sigma, Student-t critical values with df = n-1) '
      + 'and "proportion" (Wilson score interval, reliable even for small samples and extreme proportions). '
      + 'Give mean/sigma (or sd), sample size n, and optional confidence level (default 0.95). '
      + 'Use this instead of mental arithmetic for standard errors and critical values.',
    parameters: {
      kind: { type: 'string', required: true, enum: ['mean_z', 'mean_t', 'proportion'], description: 'Which interval to build.' },
      mean: { type: 'number', description: 'Sample mean (mean_z / mean_t kinds).' },
      sigma: { type: 'number', description: 'Known population standard deviation, > 0 (mean_z kind).' },
      sd: { type: 'number', description: 'Sample standard deviation, > 0 (mean_t kind).' },
      n: { type: 'number', description: 'Sample size: integer >= 1 (mean_z / proportion) or >= 2 (mean_t).' },
      successes: { type: 'number', description: 'Observed successes, integer in [0, n] (proportion kind).' },
      conf: { type: 'number', description: 'Confidence level in (0, 1) (default 0.95).' },
    },
    output: {
      schema: {
        type: 'object',
        additionalProperties: false,
        properties: {
          valid: { type: 'boolean', required: true },
          kind: { type: 'string' },
          lower: { type: 'number' },
          upper: { type: 'number' },
          marginOfError: { type: 'number' },
          criticalValue: { type: 'number' },
          standardError: { type: 'number' },
          center: { type: 'number' },
          n: { type: 'number' },
          conf: { type: 'number' },
          error: { type: 'string' },
        },
      },
      render: (_args: unknown, value: unknown) => [{ type: 'text', text: renderConfidence(value) }],
    },
    async execute(args: {
      kind: string
      mean?: number
      sigma?: number
      sd?: number
      n?: number
      successes?: number
      conf?: number
    }): Promise<ConfidenceResult> {
      return confidenceInterval(args)
    },
  })

  const event_probability = defineTool({
    name: 'event_probability',
    description: 'Event probability identities. Six operations: "union" (P(A or B); independent default true, '
      + 'false treats events as mutually exclusive), "intersection" (P(A and B), independent events), '
      + '"conditional" (P(A|B) = P(A and B)/P(B)), "bayes" (posterior P(H|E) from prior, true-positive and '
      + 'false-positive rates), "complement" (1 - p) and "at_least_one" (1 - (1-p)^n over n trials). '
      + 'All probabilities in [0, 1]. Use this instead of mental arithmetic for compound event questions.',
    parameters: {
      operation: {
        type: 'string',
        required: true,
        enum: ['union', 'intersection', 'conditional', 'bayes', 'complement', 'at_least_one'],
        description: 'Which identity to apply.',
      },
      p: { type: 'number', description: 'Event probability in [0, 1] (complement / at_least_one operations).' },
      pA: { type: 'number', description: 'P(A) in [0, 1] (union / intersection operations).' },
      pB: { type: 'number', description: 'P(B) in [0, 1] (union / intersection operations).' },
      independent: { type: 'boolean', description: 'Whether A and B are independent (union / intersection; default true).' },
      pIntersection: { type: 'number', description: 'P(A and B) in [0, 1] (conditional operation).' },
      pCondition: { type: 'number', description: 'P(B) in (0, 1] (conditional operation).' },
      prior: { type: 'number', description: 'P(H) in [0, 1] (bayes operation).' },
      truePositive: { type: 'number', description: 'P(E|H) in [0, 1] (bayes operation).' },
      falsePositive: { type: 'number', description: 'P(E|not H) in [0, 1] (bayes operation).' },
      n: { type: 'number', description: 'Number of trials, integer >= 1 (at_least_one operation).' },
    },
    output: {
      schema: {
        type: 'object',
        additionalProperties: false,
        properties: {
          valid: { type: 'boolean', required: true },
          operation: { type: 'string' },
          result: { type: 'number' },
          formula: { type: 'string' },
          error: { type: 'string' },
        },
      },
      render: (_args: unknown, value: unknown) => [{ type: 'text', text: renderEvent(value) }],
    },
    async execute(args: {
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
    }): Promise<EventResult> {
      return eventProbability(args)
    },
  })

  const hypothesis_test = defineTool({
    name: 'hypothesis_test',
    description: 'One-sample hypothesis tests with the p-value and the decision spelled out. Three kinds: '
      + '"mean_z" (mean, known population sigma), "mean_t" (mean, unknown sigma — Student-t with df = n - 1) '
      + 'and "proportion" (single proportion against p0, normal approximation). Give the sample statistic, the '
      + 'null value and n; optional tail ("two-sided" default, "left", "right") and alpha (default 0.05). '
      + 'Returns the test statistic, the standard error, all three p-values, the tail-specific critical value '
      + 'and whether the null is rejected. Use this instead of recalling critical values or doing p-value '
      + 'arithmetic by hand.',
    parameters: {
      kind: {
        type: 'string',
        required: true,
        enum: ['mean_z', 'mean_t', 'proportion'],
        description: 'Which one-sample test to run.',
      },
      sampleMean: { type: 'number', description: 'Observed sample mean (mean_z / mean_t kinds).' },
      mu0: { type: 'number', description: 'Null-hypothesis mean mu0 (mean_z / mean_t kinds).' },
      sigma: { type: 'number', description: 'Known population standard deviation, > 0 (mean_z kind).' },
      sd: { type: 'number', description: 'Sample standard deviation, > 0 (mean_t kind).' },
      n: { type: 'number', description: 'Sample size: integer >= 1 (mean_z / proportion) or >= 2 (mean_t).' },
      successes: { type: 'number', description: 'Observed successes, integer in [0, n] (proportion kind).' },
      p0: { type: 'number', description: 'Null-hypothesis proportion, strictly in (0, 1) (proportion kind).' },
      tail: { type: 'string', enum: ['two-sided', 'left', 'right'], description: 'Which p-value to report (default two-sided).' },
      alpha: { type: 'number', description: 'Significance level in (0, 1) (default 0.05).' },
    },
    output: {
      schema: {
        type: 'object',
        additionalProperties: false,
        properties: {
          valid: { type: 'boolean', required: true },
          kind: { type: 'string' },
          tail: { type: 'string' },
          statistic: { type: 'number' },
          standardError: { type: 'number' },
          df: { type: 'number' },
          pValue: { type: 'number' },
          pValueTwoSided: { type: 'number' },
          pValueLeft: { type: 'number' },
          pValueRight: { type: 'number' },
          criticalValue: { type: 'number' },
          alpha: { type: 'number' },
          reject: { type: 'boolean' },
          note: { type: 'string' },
          error: { type: 'string' },
        },
      },
      render: (_args: unknown, value: unknown) => [{ type: 'text', text: renderHypothesis(value) }],
    },
    async execute(args: {
      kind: string
      sampleMean?: number
      mu0?: number
      sigma?: number
      sd?: number
      n?: number
      successes?: number
      p0?: number
      tail?: string
      alpha?: number
    }): Promise<HypothesisResult> {
      return hypothesisTest(args)
    },
  })

  const sample_size = defineTool({
    name: 'sample_size',
    description: 'Sample size needed to hit a target margin of error. Two kinds: "mean_z" (mean estimate with '
      + 'known sigma: n = (z*sigma/E)^2) and "proportion" (proportion estimate; method "wald" default uses '
      + 'n = z^2*p*(1-p)/E^2, method "wilson" returns the smallest n whose Wilson half-width meets E). Give '
      + 'marginOfError (E > 0) and optionally sigma, p (default 0.5, the conservative worst case) and conf '
      + '(default 0.95). Returns the required integer n plus the exact uncorrected value. Use this instead of '
      + 'guessing study sizes.',
    parameters: {
      kind: {
        type: 'string',
        required: true,
        enum: ['mean_z', 'proportion'],
        description: 'Which quantity is being estimated.',
      },
      sigma: { type: 'number', description: 'Known population standard deviation, > 0 (mean_z kind).' },
      p: { type: 'number', description: 'Expected proportion, strictly in (0, 1) (proportion kind; default 0.5).' },
      method: { type: 'string', enum: ['wald', 'wilson'], description: 'Sample-size formula for the proportion kind (default wald).' },
      marginOfError: { type: 'number', required: true, description: 'Half-width of the interval you are willing to accept, > 0.' },
      conf: { type: 'number', description: 'Confidence level in (0, 1) (default 0.95).' },
    },
    output: {
      schema: {
        type: 'object',
        additionalProperties: false,
        properties: {
          valid: { type: 'boolean', required: true },
          kind: { type: 'string' },
          method: { type: 'string' },
          n: { type: 'number' },
          nExact: { type: 'number' },
          criticalValue: { type: 'number' },
          marginOfError: { type: 'number' },
          conf: { type: 'number' },
          assumedP: { type: 'number' },
          note: { type: 'string' },
          error: { type: 'string' },
        },
      },
      render: (_args: unknown, value: unknown) => [{ type: 'text', text: renderSampleSize(value) }],
    },
    async execute(args: {
      kind: string
      sigma?: number
      p?: number
      method?: string
      marginOfError?: number
      conf?: number
    }): Promise<SampleSizeResult> {
      return sampleSize(args)
    },
  })

  return { dist_calc, z_score, confidence_interval, event_probability, hypothesis_test, sample_size }
}
