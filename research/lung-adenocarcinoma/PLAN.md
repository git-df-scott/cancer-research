# EGFR-mutant lung adenocarcinoma: residual disease and relapse

Status: initial evidence audit completed on Scott's subsequent authorization, 2026-09-07. See `AUDIT.md`, `SOURCES.json`, and `SEARCH_LOG.md`. No model implementation, biological-data analysis, simulations, or intervention experiments have started. Scott explicitly requested stopping after the audit; resume only when asked. No overnight or scheduled work is authorized.

Branch: `codex/lung-adenocarcinoma-persisters`

## Objective

Work toward durable cancer eradication by understanding and experimentally targeting the malignant cells that survive treatment. The initial scope is EGFR-mutant lung adenocarcinoma under osimertinib exposure. Select the exact genotype and experimental system after inspecting available data.

Research question: Can an experimentally supported intervention strategy eliminate residual malignant populations before resistant regrowth, under realistic exposure and toxicity constraints?

Near-term success is one independently testable prediction about residual-cell survival or regrowth. A simulated zero-cell count is not evidence of patient cure. A well-supported negative finding or demonstration that the available data cannot distinguish mechanisms is also useful.

## Working definitions and competing explanations

- Treatment-sensitive cells: cells susceptible under the specified assay and exposure conditions.
- Drug-tolerant persisters: survivors whose tolerance may be reversible; distinguish dividing and non-dividing populations where measurements support it.
- Stable resistance: resistance persisting under a defined withdrawal/rechallenge assessment; do not automatically equate it with a known genetic mutation.
- Residual disease: viable malignant cells remaining after initial treatment response.
- Durable control: a prespecified period without regrowth following a defined observation or withdrawal protocol. It does not imply lifelong eradication.

Test selection of pre-existing states against treatment-induced adaptation and mixtures of both. Do not assume that all relapse passes through a single persister state, or that a small state-based model adequately represents every experimental system.

## Planned sequence

1. **Inspect evidence and data access.** Read original methods, supplements, source data and repository instructions. Inventory genotype, model system, treatment exposure, time points, biological replicates, lineage measurements, cell counts, withdrawal/rechallenge results and missing observations. Confirm download permissions and usable formats before choosing the benchmark.
2. **Audit prior art.** Separate established persister mechanisms from unresolved questions. Treat FAK–YAP/TEAD targeting as an existing experimental benchmark, not our discovery. Select one specific gap only after this audit.
3. **Choose an identifiable baseline.** Start with the smallest model supported by the measurements. Compare alternative explanations; add spatial structure only if data and the question require it. Distinguish stochastic population extinction from numerical truncation.
4. **Write an analysis plan before fitting.** Fix endpoints, calibration/validation separation, uncertainty treatment, controls and failure criteria. Split at the level of independent experiments or biological models where possible, rather than randomly splitting correlated cells.
5. **Reproduce an external experiment.** Require agreement with initial killing, residual populations and regrowth where available. Track failure as well as success. Do not tune toward an attractive intervention ranking.
6. **Validate and challenge.** Test independent observations, parameter identifiability, alternate state-transition assumptions, pre-existing resistance, observation limits and model mismatch. Stop intervention ranking if plausible models make incompatible predictions unsupported by discriminating data.
7. **Evaluate a bounded intervention question.** Use measured perturbation effects and pharmacologically plausible exposure. Assess normal-cell effects where evidence exists; explicitly mark unavailable safety information. Compare against treatment alone and relevant established combinations. A lower residual burden alone does not prove prevention of relapse.
8. **Prepare one laboratory validation proposal.** State the discriminating prediction, controls, viable-cell and regrowth endpoints, and results that would falsify it. Specialist collaboration is needed to assess feasibility. Do not contact anyone without Scott's authorization.

## Initial sources to inspect

These sources motivated the project; full methods, data availability and novelty checks remain to be done.

- Oren et al. *Cycling cancer persister cells arise from lineages with distinct programs*. Nature (2021). https://doi.org/10.1038/s41586-021-03796-6 . Starting evidence for lineage-resolved dividing and non-dividing persister behaviour.
- *Focal adhesion kinase-YAP signaling axis drives drug-tolerant persister cells and residual disease in lung cancer*. Nature Communications (2024). https://doi.org/10.1038/s41467-024-47423-0 . Starting preclinical benchmark for perturbing persistence and assessing regrowth; not clinical cure evidence.

## Reproducibility and relationship to the existing project

Keep lung-project artifacts under `research/lung-adenocarcinoma/` initially. Preserve the lymphoma work and its results. Do not repurpose its model without an explicit justification based on lung-cancer biology.

Record source versions and access dates, data provenance, environment, seeds, parameter choices, failed checks and deviations. Never overwrite original data. Bound compute and preserve results incrementally. Separate observations, calculations, hypotheses and interpretations in all reports.

## Next session, only if authorized

Read `AUDIT.md` and applicable repository instructions. The audit supersedes the original generic simulator starting point: ResMap already provides a cross-context benchmark. Consider whether the BRD4 discrepancy between ResMap and sequential-treatment work supports a matched comparison, with Oren's lineage findings as a constraint on interpreting survivors. Further protocol appraisal, data analysis and modeling require Scott's next instruction. Do not begin a large simulation sweep merely because compute is available.
