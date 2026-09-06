/**
 * Tool definitions for dsh-probstat: four deterministic probability and
 * statistical-inference tools exposed to every agent via defineTool.
 *
 * @module dsh-probstat/tools
 */
import { defineTool } from '@deepseek-ai/dsh-tools';
import { confidenceInterval, distCalc, eventProbability, zScore, } from "./core.js";
function fmt(value) {
    if (!Number.isFinite(value))
        return String(value);
    if (value !== 0 && (Math.abs(value) >= 1e7 || Math.abs(value) < 1e-4)) {
        return value.toExponential(6);
    }
    return String(Number(value.toPrecision(8)));
}
/* ------------------------------------------------------------------ */
/* render helpers                                                      */
/* ------------------------------------------------------------------ */
function renderDistCalc(value) {
    const v = value;
    if (!v.valid)
        return `dist_calc failed: ${v.error}`;
    const lines = [`${v.distribution} distribution (operation: ${v.operation})`];
    if (v.result !== undefined)
        lines.push(`result: ${fmt(v.result)}`);
    if (v.mean !== undefined)
        lines.push(`mean: ${fmt(v.mean)}`);
    if (v.variance !== undefined)
        lines.push(`variance: ${fmt(v.variance)}`);
    if (v.sd !== undefined)
        lines.push(`sd: ${fmt(v.sd)}`);
    if (v.note !== undefined)
        lines.push(`note: ${v.note}`);
    return lines.join('\n');
}
function renderZScore(value) {
    const v = value;
    if (!v.valid)
        return `z_score failed: ${v.error}`;
    const lines = [`mode: ${v.mode}`];
    if (v.z !== undefined)
        lines.push(`z: ${fmt(v.z)}`);
    if (v.z1 !== undefined && v.z2 !== undefined)
        lines.push(`z range: [${fmt(v.z1)}, ${fmt(v.z2)}]`);
    if (v.probability !== undefined)
        lines.push(`probability: ${fmt(v.probability)}`);
    if (v.percentile !== undefined)
        lines.push(`percentile: ${fmt(v.percentile)} %`);
    if (v.pValue !== undefined)
        lines.push(`two-sided p-value: ${fmt(v.pValue)}`);
    return lines.join('\n');
}
function renderConfidence(value) {
    const v = value;
    if (!v.valid)
        return `confidence_interval failed: ${v.error}`;
    const kind = v.kind === 'proportion' ? 'proportion (Wilson)' : v.kind;
    const lines = [`${kind} ${fmt((v.conf ?? 1) * 100)}% confidence interval`];
    lines.push(`[${fmt(v.lower)} , ${fmt(v.upper)}]`);
    if (v.center !== undefined)
        lines.push(`center: ${fmt(v.center)}`);
    lines.push(`margin of error: ${fmt(v.marginOfError)}`);
    lines.push(`critical value: ${fmt(v.criticalValue)}`);
    lines.push(`standard error: ${fmt(v.standardError)}`);
    return lines.join('\n');
}
function renderEvent(value) {
    const v = value;
    if (!v.valid)
        return `event_probability failed: ${v.error}`;
    const lines = [`${v.operation}: ${fmt(v.result)}`];
    if (v.formula !== undefined)
        lines.push(`formula: ${v.formula}`);
    return lines.join('\n');
}
/* ------------------------------------------------------------------ */
/* tool registration                                                   */
/* ------------------------------------------------------------------ */
export function buildProbstatTools() {
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
            render: (_args, value) => [{ type: 'text', text: renderDistCalc(value) }],
        },
        async execute(args) {
            return distCalc(args);
        },
    });
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
            render: (_args, value) => [{ type: 'text', text: renderZScore(value) }],
        },
        async execute(args) {
            return zScore(args);
        },
    });
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
            render: (_args, value) => [{ type: 'text', text: renderConfidence(value) }],
        },
        async execute(args) {
            return confidenceInterval(args);
        },
    });
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
            render: (_args, value) => [{ type: 'text', text: renderEvent(value) }],
        },
        async execute(args) {
            return eventProbability(args);
        },
    });
    return { dist_calc, z_score, confidence_interval, event_probability };
}
//# sourceMappingURL=tools.js.map