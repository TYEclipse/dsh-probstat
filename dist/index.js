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
import z from '@deepseek-ai/schemastery';
import { buildProbstatTools } from "./tools.js";
/** Stable Cordis plugin name (also the config key under `plugins:`). */
export const name = 'dsh-probstat';
/** Services required before tool registration can start. */
export const inject = ['agents', 'tools'];
export const Config = z.object({});
/** Register every probstat tool on one agent; returns the disposer. */
function decorate(agent, tools) {
    const disposers = Object.values(tools).map((definition) => agent.ctx.tools.register(definition));
    return () => {
        for (const dispose of disposers) {
            try {
                dispose();
            }
            catch {
                // already disposed
            }
        }
    };
}
/** Mount the probstat tools on every live agent and every future one. */
export function apply(ctx, _config) {
    const tools = buildProbstatTools();
    const disposers = new Set();
    const decorateAgent = (agent) => {
        try {
            disposers.add(decorate(agent, tools));
        }
        catch (error) {
            ctx.logger('probstat').warn(`tool registration for agent ${agent.id} failed: ${error instanceof Error ? error.message : String(error)}`);
        }
    };
    for (const agent of ctx.agents.list())
        decorateAgent(agent);
    const off = ctx.on('agent/created', ({ agent }) => decorateAgent(agent));
    ctx.effect(() => () => {
        off();
        for (const dispose of disposers) {
            try {
                dispose();
            }
            catch {
                // already disposed
            }
        }
        disposers.clear();
    });
}
//# sourceMappingURL=index.js.map