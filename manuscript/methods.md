# Methods draft

Draft based on the frozen protocols, their amendments and the retained implementation. Accountable authors must review the final manuscript and declarations.

**Task and data**

The prediction task maps source-context perturbation response profiles to target-context summary effects for perturbations already measured in the source. The primary confirmation uses Replogle K562 CRISPRi source profiles and author-provided primary human CD4 T differential-expression estimates in Rest, Stim8hr and Stim48hr. Target `log_fc` is DESeq2 log2 fold change; no additional effect-layer transform is applied. The author data schema is documented [here](https://github.com/emdann/GWT_perturbseq_analysis_2025/blob/master/metadata/data_sharing_readme.md).

The final confirmation population comprises 741 eligible perturbation targets and 6,851 common response genes. A metadata-stratified partition assigns 256 targets to calibration anchors and 485 to held-out evaluation, shared across all three conditions. The original 1,100-target availability threshold failed; the final analysis continued after an outcome-blind metadata amendment lowered the feasibility threshold to 700. Missing allowed metadata rows were excluded without imputation. This is blinded metadata-adaptive confirmation, not an unchanged preregistration.

**Calibration families and anchor decisions**

The five candidates are source-feature scalar calibration, direct-program calibration, low-rank Ridge transfer, low-rank plus program-residual transfer, and gene-wise Ridge calibration. Source-derived gene programs use a 32-component source embedding and 16-cluster K-means. The low-rank family uses up to 16 components. Model settings, coefficient bounds and numerical seeds are specified in `configs/cd4t_confirmation.yaml`.

Target transformations and fitting use training anchors only. Residual-program coefficients are estimated against inner cross-fitted low-rank residuals. Source representations may use all eligible source profiles because those are available at deployment. Targets whose gene symbol is absent from the common response axis receive a deterministic source-only program assignment; 52 eligible targets use that rule.

Anchors undergo five-fold out-of-fold evaluation. Per-target loss is MSE on the 200 response genes with the largest absolute source effect. Twenty thousand target bootstraps estimate each family's frequency of minimum mean OOF loss. If the original mean-loss winner has frequency at least 0.80, APMA uses that family alone. Otherwise, it averages the prediction vectors using the bootstrap frequencies. No scalar baseline is forced into an ambiguous-case average. These frequencies describe resampling stability conditional on the anchor data and candidate set.

**Primary inference**

Within each condition, regret is mean APMA held-out loss divided by the minimum mean loss among fixed candidate families, minus one. The fixed-family oracle is recomputed within each of 20,000 held-out-target bootstraps. Confirmation requires the upper endpoint of the percentile 95% interval to be at most 0.02 in all three conditions. Negative regret is possible because an average of prediction vectors may outperform each fixed family. Secondary comparisons and Holm corrections are specified in the active protocol and byte-preserved original protocol.

**Replicate consistency**

Predictions fitted to pooled anchors are evaluated against held-out guide-level and donor-pair effects without fitting on those replicate outcomes. The retained replicate implementation refits the pooled-anchor families using its own fixed seed offsets (2000 plus the condition index); it does not load the exact fitted prediction matrices from the primary confirmation, which uses offsets 999 plus the condition index. Family definitions, hyperparameters and the anchor/held-out partition are inherited unchanged. This is a consistency check of the refitted family contrast, rather than an identity check of stored primary predictions. The primary comparison is residual-program versus low-rank MSE. Guide consensus averages the two guide losses within target; donor consensus averages the six donor-pair losses within eligible target. Bootstrap resamples preserve complete target replicate vectors. Guide-consensus target counts are 255, 285 and 258; donor-consensus counts are 183, 217 and 198. The original feasibility thresholds failed and were amended without reading replicate effect values. Donor-pair modalities share donors and are not six independent replications.

**Sensitivity analyses**

Source-defined OOD evaluates fixed source functional clusters with all anchors in the tested cluster removed. The original eight-cluster partition failed the minimum held-out-size rule; a source-only feasibility adaptation selected seven clusters while retaining the original size thresholds. This tests source-geometry extrapolation within the same study.

Fallback exclusion is a post hoc sensitivity removing 36 held-out fallback targets. The anchor-budget analysis evaluates nested source-only farthest-first subsets of 32, 64, 128 and 256 anchors on the same 485 held-out targets. Both regret and relative MSE increase versus 256 anchors use a 5% upper-confidence-bound rule. The first anchor-budget run was invalid because fold ordering and seeds prevented the 256-anchor control from reproducing confirmation. The correction and invalid outputs are retained. Thirty-two is the smallest tested passing budget; lower budgets were not evaluated.

**Evidence scope**

These sensitivity analyses reuse the primary external study. Target-level uncertainty is not between-study uncertainty. Replicate-consistent improvement supports a reproducible component but does not uniquely separate biological signal and measurement noise. See the original protocols and amendments under `docs/provenance/` and the complete development history under `docs/research_history.html`.
