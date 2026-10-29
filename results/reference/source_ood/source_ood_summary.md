# V38 source-defined OOD 压力测试

- 状态：**V38_SOURCE_OOD_SENSITIVITY_SUPPORTS_SELECTOR**
- 三状态 5% regret 联合门槛：**通过**。

## Rest

- regret：-0.117%（95% CI -0.213% 至 0.001%）
- residual 胜过 low-rank 的 source clusters：7/7

## Stim8hr

- regret：-0.446%（95% CI -0.572% 至 -0.083%）
- residual 胜过 low-rank 的 source clusters：7/7

## Stim48hr

- regret：-0.255%（95% CI -0.304% 至 0.356%）
- residual 胜过 low-rank 的 source clusters：7/7

该分析只支持 source geometry 外推敏感性，不等同于按目标真值聚类的 EL-PERTURB OOD，也不构成新的独立研究。
