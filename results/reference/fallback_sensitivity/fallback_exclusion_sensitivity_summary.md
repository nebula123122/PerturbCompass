# V39：排除回退靶点后的敏感性分析

- 状态：**V39_V36_CONCLUSION_ROBUST_TO_FALLBACK_TARGET_EXCLUSION**
- 三种状态仍全部通过 V36 的 2% regret 门槛：**是**。
- 这是看到 V36 结果后做的稳健性检查，不是新的盲测。

## Rest

- 保留 449 个，删除 36 个回退靶点。
- selector regret：-0.887%（95% CI -0.932% 至 -0.280%）。
- selector 相对 scalar 改善：2.055%。
- residual 相对 low-rank 改善：1.191%。

## Stim8hr

- 保留 449 个，删除 36 个回退靶点。
- selector regret：-1.037%（95% CI -1.195% 至 -0.492%）。
- selector 相对 scalar 改善：2.559%。
- residual 相对 low-rank 改善：1.494%。

## Stim48hr

- 保留 449 个，删除 36 个回退靶点。
- selector regret：0.253%（95% CI -0.000% 至 1.033%）。
- selector 相对 scalar 改善：1.054%。
- residual 相对 low-rank 改善：1.006%。

