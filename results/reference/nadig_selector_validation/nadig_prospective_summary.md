# V35 Nadig 双屏前瞻确认

- 状态：**V35_PROSPECTIVE_SELECTOR_NOT_CONFIRMED_BOTH_SCREENS**
- 设计：256 个 anchor 只用于选法；1004 个 held-out 扰动只打开一次。
- 联合主门槛：**未通过**（两屏 regret 95% CI 上界均 ≤2%）。

## hepg2

- anchor 决策：uncertain_equivalence_ensemble；成员：scalar, residual_program
- anchor 赢家：residual_program；bootstrap 稳定率：68.54%
- test 最佳固定家族：residual_program
- selector regret：3.767%（95% CI 1.046% 至 6.569%）
- selector 相对 scalar：28.553%（95% CI 26.742% 至 30.240%；Holm p=1.99998e-05）

## jurkat

- anchor 决策：confident_single_family；成员：residual_program
- anchor 赢家：residual_program；bootstrap 稳定率：92.00%
- test 最佳固定家族：residual_program
- selector regret：0.000%（95% CI -0.000% 至 0.000%）
- selector 相对 scalar：25.061%（95% CI 20.952% 至 29.194%；Holm p=1.99998e-05）

HepG2 与 Jurkat 是同一外部研究中的两个不同细胞屏，不等于两项独立研究；目标级置信区间只描述屏内扰动稳定性。
