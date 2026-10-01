/**
 * Tool definitions for dsh-probstat: four deterministic probability and
 * statistical-inference tools exposed to every agent via defineTool.
 *
 * @module dsh-probstat/tools
 */
import { type ToolDefinition } from '@deepseek-ai/dsh-tools';
export interface ToolSet {
    dist_calc: ToolDefinition;
    z_score: ToolDefinition;
    confidence_interval: ToolDefinition;
    event_probability: ToolDefinition;
    hypothesis_test: ToolDefinition;
    sample_size: ToolDefinition;
}
export declare function buildProbstatTools(): ToolSet;
//# sourceMappingURL=tools.d.ts.map