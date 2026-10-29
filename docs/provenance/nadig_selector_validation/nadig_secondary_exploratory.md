# V35 结果后次要分析（探索性）

**这些分析是在 V35 主结果已知后进行，不能冒充前瞻确认。**

## hepg2

- residual_program 相对 low_rank：2.699%（95% CI 1.988% 至 3.381%；Holm p=7.99992e-05）。
- 若事后删除强制 scalar：regret 0.000%（95% CI -0.000% 至 0.000%）。
- 若事后平均 residual+low_rank：regret 1.387%（95% CI 1.019% 至 1.760%）。

## jurkat

- residual_program 相对 low_rank：3.322%（95% CI 2.613% 至 4.000%；Holm p=7.99992e-05）。
- 若事后删除强制 scalar：regret 0.000%（95% CI -0.000% 至 0.000%）。
- 若事后平均 residual+low_rank：regret 1.718%（95% CI 1.344% 至 2.083%）。

可用于下一轮预注册的候选修订：不确定时只在接近最优的强家族之间集成，不再无条件加入 scalar。
