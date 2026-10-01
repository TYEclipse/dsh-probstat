/**
 * dsh-probstat — deterministic probability & statistical inference math for
 * DeepSeek Harness.
 *
 * Six zero-dependency tools (pure arithmetic):
 *   dist_calc           — pdf / cdf / survival / quantile / stats for normal,
 *                         binomial, poisson, exponential, uniform, geometric
 *   z_score             — standard-normal table math (z <-> probability)
 *   confidence_interval — mean (z / t) and proportion (Wilson) intervals
 *   event_probability   — union / intersection / conditional / Bayes /
 *                         complement / at-least-one identities
 *   hypothesis_test     — one-sample z / t / proportion tests with p-values,
 *                         critical values and the alpha decision
 *   sample_size         — required n for a target margin of error (mean with
 *                         known sigma; proportion via Wald or Wilson)
 *
 * Agents get distribution math wrong in predictable ways: hallucinated
 * z-table values, normal cdf/quantile off by decimals, wrong critical values
 * and sloppy compound-event identities. These tools replace all of it with
 * deterministic formulas.
 *
 * @module dsh-probstat
 */
import type { Context } from '@deepseek-ai/cordis';
import z from '@deepseek-ai/schemastery';
/** Stable Cordis plugin name (also the config key under `plugins:`). */
export declare const name = "dsh-probstat";
/** Services required before tool registration can start. */
export declare const inject: string[];
/** Plugin configuration (no tunables today — kept for schema completeness). */
export interface Config {
    /** Reserved for future options. */
    _?: never;
}
export declare const Config: z<Config>;
/** Mount the probstat tools on every live agent and every future one. */
export declare function apply(ctx: Context, _config: Config): void;
//# sourceMappingURL=index.d.ts.map