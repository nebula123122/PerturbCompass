# V36 原代 CD4 T 细胞三状态前瞻确认

- 状态：**V36_PROSPECTIVE_SELECTOR_CONFIRMED_ALL_THREE_CONTEXTS**
- 联合门槛：**通过**（三种状态 regret 95% CI 上界均 ≤2%）。

## Rest

- anchor 决策：ambiguous_probability_model_average；权重：direct_program=0.749, scalar=0.000, residual_program=0.159, low_rank=0.000, gene_ridge=0.092
- anchor 赢家：direct_program；稳定率：74.87%
- held-out 最佳固定家族：residual_program
- selector regret：-0.819%（95% CI -0.949% 至 -0.223%）
- selector 相对 scalar：2.057%
- residual 相对 low-rank：1.146%

## Stim8hr

- anchor 决策：ambiguous_probability_model_average；权重：direct_program=0.651, scalar=0.000, residual_program=0.312, low_rank=0.000, gene_ridge=0.037
- anchor 赢家：direct_program；稳定率：65.10%
- held-out 最佳固定家族：residual_program
- selector regret：-0.942%（95% CI -1.203% 至 -0.401%）
- selector 相对 scalar：2.589%
- residual 相对 low-rank：1.441%

## Stim48hr

- anchor 决策：dominant_single_family；权重：direct_program=1.000, scalar=0.000, residual_program=0.000, low_rank=0.000, gene_ridge=0.000
- anchor 赢家：direct_program；稳定率：82.45%
- held-out 最佳固定家族：residual_program
- selector regret：0.240%（95% CI -0.000% 至 0.992%）
- selector 相对 scalar：1.004%
- residual 相对 low-rank：0.971%

三种状态和四位供者来自同一项研究，因此联合门槛检验跨状态稳定性，不等于三项独立外部复制。
