# Results draft

Values transcribed from byte-preserved result records. This is a manuscript working draft, not a new analysis or publication approval.

**Primary CD4 T confirmation**

| Condition | APMA decision | Regret | 95% interval |
| --- | --- | ---: | --- |
| Rest | ambiguous_probability_model_average | -0.819% | -0.949% to -0.223% |
| Stim8hr | ambiguous_probability_model_average | -0.942% | -1.203% to -0.401% |
| Stim48hr | dominant_single_family | 0.240% | -0.000% to 0.992% |

All three condition-specific upper bounds meet the 2% rule. Each condition has 485 held-out perturbation targets. They come from one external study.

**Anchor budget**

| Condition | 32 versus 256 anchors: relative MSE increase | 95% interval |
| --- | ---: | --- |
| Rest | 1.359% | 1.159% to 1.573% |
| Stim8hr | 1.563% | 1.257% to 1.942% |
| Stim48hr | 0.356% | 0.099% to 0.614% |

Thirty-two anchors passes both declared 5% rules in all three conditions, as do the larger tested budgets. It is the smallest tested passing budget, not an identified universal minimum.

**Supporting and negative evidence**

The source-defined seven-cluster OOD sensitivity and fallback exclusion remain supportive. The replicate analysis evaluates residual-program versus low-rank consistency; it is not a donor-disjoint refit of the complete APMA pipeline. Large-scale K562 validation and Feng external validation support program-residual signal, while the Adamson primary gate and initial Nadig selector joint gate did not pass. Their complete records remain in the research history and provenance.

The current evidence does not estimate between-study generalization variance, establish unique biology/noise separation, or replace direct comparisons with modern models and broader biological metrics.
