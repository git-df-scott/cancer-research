# EGFR-mutant lung adenocarcinoma: what is established, what remains unresolved

Audit date: 2026-09-07. Scope: osimertinib-associated residual disease and relapse, principally common EGFR exon 19 deletion and L858R contexts. This is a focused, reproducible evidence audit, not a systematic review or a patent clearance opinion. Research decisions, not patient-specific treatment advice.

**Status: audit complete; stopped at Scott's request. No model code, simulations, or biological-data analysis started.**

## Main finding

There are credible routes toward eliminating residual disease, with eradication reported in selected mouse models. The unresolved problem is obtaining durable, safe eradication across heterogeneous human disease. Neither targeting persisters, adding a second drug after initial shrinkage, comparing targets across contexts, nor fitting sensitive–persister–resistant populations is intrinsically new.

The most consequential new source for this project is **ResMap (Sun et al., Science Advances, June 2026)**: 94 compounds, 57 previously reported targets, four lung-cancer models and two oxygen conditions. It already supplies the standardized comparison we might otherwise have proposed. Twelve broadly active targets emerged initially; nine reproduced in follow-up validation. Its experiments are in vitro, and its tolerability proxies are not human safety measurements. [S01]

**Recommendation for a future session:** investigate a bounded discrepancy between existing studies before building a simulator. BRD4 inhibition is beneficial in selected sequential-treatment experiments, whereas ResMap reports pro-persistence effects in some PC9 conditions. Establish whether assay timing, cell state, concentration, oxygen, or measurement differences explain the discrepancy. This is a candidate research question, not an established new mechanism or a proven novel project. [S01, S02]

## Exact question, assumptions and definitions

Can a treatment strategy eliminate viable residual malignant populations before resistant regrowth at tolerable exposure?

Start with a defined experimental system; do not combine every EGFR genotype, treatment line, disease stage or anatomical site into one calibration. A reversible tolerant phenotype, stable resistance, and histological transformation are distinct possibilities. A high expression score does not establish a functional dependency. A surface marker can be useful for recognition without being necessary for survival.

Residual disease means viable cancer remaining after treatment. Molecular residual disease detected in blood is an assay-defined observation, not a direct count of all residual cells. No detected signal does not establish absence of cancer. For this project, a useful result is an independent, falsifiable prediction about survival or regrowth. A cure requires clinical evidence of durable freedom from disease; mouse non-regrowth and simulated extinction do not establish it.

## Evidence ledger: what has already been done

| Approach | Primary evidence and result | Limits and next falsification check |
|---|---|---|
| Stronger first-line systemic treatment | FLAURA2 randomized 557 patients: median overall survival 47.5 versus 37.6 months for osimertinib plus chemotherapy versus osimertinib; HR 0.77, 95% CI 0.61–0.96. Grade ≥3 adverse events: 70% versus 34%. MARIPOSA compared 429 patients per principal arm: amivantamab–lazertinib improved OS, HR 0.75, 95% CI 0.61–0.92; grade ≥3 adverse events 80% versus 52%. [S03, S04] | Strong clinical evidence for benefit with added toxicity. These trials do not demonstrate universal eradication or identify which persister mechanism caused benefit. Do not compare survival numbers between trials as if randomized head-to-head. |
| Identify different surviving lineages | Oren et al. distinguished cycling and non-cycling persisters in lineage-tracked experiments. In their PC9 experiments, GPX4 inhibition reduced total persisters but enriched the cycling fraction among survivors. [S05] | Strong evidence within tested models. A larger fraction is not necessarily a larger absolute number. Test absolute survivor counts and subsequent regrowth, not just proportions. |
| Sequential and upfront drug combinations | Criscione et al. screened PC9 and six additional EGFR-mutant lines. BRD4 inhibition was more effective sequentially in selected lines; other targets favored upfront treatment. AURKB combinations delayed regrowth in LU5221/PC9 models but did not show that benefit in H1975 xenografts. [S02] | Direct counterexample to assuming one useful target works in all models. Assay timing and exposure must be matched before pooling effect sizes. |
| FAK–YAP/TEAD | Haderk et al. and an independent functional-genomics study supported this pathway using perturbations, organoids and/or animal systems. Haderk included immune-context experiments and withdrawal/regrowth experiments. [S06, S07] | Already established preclinical direction, not a new target. Small model-specific experiments and shared human datasets limit generalization. Independently test failures and normal-tissue effects. |
| Early, brief IGF-1R inhibition | Wang et al. reported eradication and no regrowth after cessation in selected AXL-low cell-derived and patient-derived xenografts. [S08] | Already tests a finite intervention aimed at durable elimination. Biomarker-defined, small mouse models; not all EGFR cancers and not human cure. Test outside the selected AXL-low context. |
| Restore apoptosis | Navitoclax plus osimertinib reached a phase Ib study, n=27. The 100% response rate applies to nine expansion patients, not the whole cohort. Prior-osimertinib patients did not respond in this study. Platelet effects demonstrate a real normal-cell constraint. [S09] | Single arm, different treatment context from first-line residual disease; efficacy cannot be assigned to navitoclax alone. Do not convert a nine-person response rate into an eradication claim. |
| Ferroptosis | Hangauer et al. established GPX4 vulnerability in persister models across cancers. Higuchi et al. later identified survivors relying on residual FSP1 and altered oxidative metabolism. [S10, S11] | An existing strategy with escape routes. Much evidence is outside EGFR lung cancer. Genetic loss, chemical inhibition, achievable exposure and safe systemic treatment are not interchangeable. |
| Immune elimination of TROP2-positive residual cells | Baldacci et al. compared a TROP2 ADC with CAR-T treatment. In the DFCI-243 residual-disease experiment, 4/10 mice had complete responses and 1/10 a durable partial response to 200 days; some animals relapsed. [S12] | Encouraging finite-follow-up result, not uniform cure. Immunodeficient xenografts do not establish human immunity, organ coverage or safety. Surviving antigen-positive disease makes access/function a competing explanation to antigen loss. |
| Improve TROP2 recognition and safety | Brea et al. engineered TROP2 CAR approaches; the authors explicitly state their binders do not react with murine TROP2, preventing evaluation of corresponding normal-tissue toxicity in those mice. Logic-gated approaches are already proposed. [S13] | More efficient killing is not enough. Establish recognition of residual human cancer versus vulnerable normal human tissues. Do not propose generic dual-antigen gating as novel. |
| TROP2 ADC treatment at the persister stage | Liao et al., Cancer Cell 2026, report DTP suppression and delayed relapse with sacituzumab tirumotecan plus osimertinib, with preliminary phase 2 activity. Abstract verified; full clinical tables not obtained. [S14] | This prevents generalizing the weaker sacituzumab-govitecan result to all TROP2 ADCs. Agents, payloads, settings and endpoints differ. No verified human cure claim here. |
| Block the supporting microenvironment | Liu et al. report damaged mitochondrial transfer from tumour cells to a fibroblast population and preclinical reversal of tolerance when transfer was disrupted. [S15] | Stromal support is already being targeted. Association in human samples and mechanistic model evidence have different strength. Check prevalence and whether benefit persists across stromal contexts. |
| Suppress adaptive evolution | Isozaki et al. link therapy-induced APOBEC3A to persistent-cell evolution; deletion delayed resistance. [S16] | Preventing one route to new mutations is not killing pre-existing resistant cells or eliminating the reservoir. |
| Local consolidation and molecular monitoring | A single-arm radiation-plus-osimertinib study enrolled 42, with 32 receiving radiation; median PFS 32.3 months, 95% CI 21.9–51.7. ADAURA's exploratory MRD subset contained 220 patients and showed molecular detection preceding imaging events by a median 4.7 months. [S17, S18] | Consolidation is already tested. Single-arm selection limits causality. ADAURA is resected disease, not metastatic persister-model validation; its molecular subset is not the entire randomized cohort. |

Other located primary directions include PRMT1, BET/BRD2, SOS1, and immune-suppressive CCL20 perturbation. These add prior art, not a reason to assemble an untested many-drug combination. Exact citations and inspection depth are in `SOURCES.json`.

## Trial records checked directly

ClinicalTrials.gov records retrieved through its API. Status is what the saved record states, not independent confirmation of activity at every site. Registration is not evidence of successful dosing, target engagement or efficacy.

| Record | What is already registered | Status in retrieved record |
|---|---|---|
| [NCT04665206](https://clinicaltrials.gov/study/NCT04665206) | VT3989 study includes an osimertinib combination | Recruiting; updated 2026-04-02 |
| [NCT05228015](https://clinicaltrials.gov/study/NCT05228015) | IK-930 includes an osimertinib combination arm | Terminated; sponsor strategic reasons; updated 2024-11-18. No posted results. This does not demonstrate failure of TEAD biology. |
| [NCT02520778](https://clinicaltrials.gov/study/NCT02520778) | Navitoclax–osimertinib | Completed; actual n=27; updated 2024-08-21 |
| [NCT04780568](https://clinicaltrials.gov/study/NCT04780568) | Tegavivint–osimertinib | Active, not recruiting; actual n=24; updated 2026-08-13 |
| [NCT05401110](https://clinicaltrials.gov/study/NCT05401110) | Carotuximab–osimertinib | Recruiting; estimated n=60; updated 2026-08-06 |

## Where a contribution could be made

### 1. Explain when an apparent persister treatment protects survivors

**Priority: highest for an initial computational investigation.** Compare the BRD4 discrepancy described above. First establish protocol comparability, then ask whether initial growth suppression, persistence maintenance and later regrowth need distinct response models. Oxygen, concentration, cell-line provenance and counting methods are alternative explanations.

**Candidate falsification:** if matched-condition measurements agree, there is no unexplained sign reversal to claim. If a state-based explanation fails to predict held-out conditions better than a simple exposure-and-growth baseline, reject it. The intended contribution is an experimentally validated boundary of benefit and harm, not another list of targets. An early protective/late vulnerable state is a hypothesis only.

### 2. Identify which survivors remain after a promising intervention

Use lineage information to distinguish fewer survivors from elimination of the lineages that regrow. The most consequential question is whether an intervention reliably removes all observed regrowth-capable populations while leaving a feasible normal-tissue margin. No assay can certify coverage of every cell in a patient.

**Candidate falsification:** apparent complementary killing disappears after measuring absolute cell counts or under a different patient-derived model. Separate-drug bulk screens cannot identify whether the same individual cells survive both drugs; missing joint measurements cannot be manufactured by simulation.

### 3. Define a selective immune-recognition window in residual disease

The TROP2 studies justify asking whether treatment creates a measurable period of useful tumour-versus-normal recognition and whether that period covers residual disease across sites. This is a translation problem already recognized by the original authors, not an untouched idea. A specific validated marker relationship or failure mechanism could still be valuable.

**Candidate falsification:** normal human tissue shares the proposed recognition pattern, or an important residual malignant population lacks it. RNA coexpression is insufficient evidence of surface-protein density or safe recognition. This track requires suitable human-tissue evidence and specialist experimental collaboration.

### What I did not establish

I did not establish that any target or combination has never been tried. The audit found substantial prior art for every broad idea considered. No complete match was established for the narrowly framed BRD4 reconciliation project, but unpublished work, patents, conference material and incompletely indexed literature remain. A generic target atlas is directly anticipated by ResMap. A generic three-population simulation also has prior art; two additional modeling papers were located but not fully verified, so no modeling novelty is claimed.

## Data feasibility and independence

| Resource | Access verified | Useful for | Critical limit |
|---|---|---|---|
| ResMap | Full text and data-availability statement inspected; supplementary benchmark tables identified | Existing reference for cross-context perturbation effects | Six-day screen, in vitro. Underlying tables not yet downloaded or analyzed. No longitudinal lineage coverage or human therapeutic index. |
| [GSE193259](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193259), Criscione | Public series manifest retrieved: 112 sample records, RNA/ATAC/ChIP data; paper supplies separate drug-screen supplements | Expression and chromatin context for upstream/sequential response experiments | 112 records are not 112 patients or independent models. Raw archive and response tables not yet analyzed. |
| [GSE150949](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150949), Oren | Public manifest: 19 sample records; processed PC9 matrix and lineage metadata listed | Lineage-resolved persistence and cycling | Barcode lineages, cells and biological replicates are different statistical units. |
| [GSE302284](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE302284), Baldacci | Public manifest: six sample records; raw archive listed | TROP2-associated residual states | Other data require author request. No permission to contact authors has been given. |
| Haderk: PRJNA766057; Maynard: PRJNA591860 | Deposits and code locations checked in full-text availability statements | Mechanistic and human residual-disease context | Repository/data contents not fully inspected. Haderk reuses Maynard patient data: these are not independent patient validations. |

Maynard's study sampled 49 biopsies from 30 patients across targeted-therapy contexts. It is not an osimertinib-only randomized intervention dataset. A model must not use individual cells as independent patient replicates or infer cell-transition rates directly from cross-sectional expression. [S19]

## Recommended decision, and stopping point

Keep the chosen cancer. Replace the initial generic simulator proposal with a future **benchmark-and-discrepancy investigation** anchored in ResMap, Criscione and Oren. The first future decision is whether their underlying observations permit a matched comparison; if not, the useful output is a precisely specified missing experiment.

The most plausible path toward cure suggested by this audit is durable elimination of residual disease with an acceptable normal-tissue margin, before additional resistance emerges. This is a research interpretation, not evidence that a particular combination will accomplish it. Existing human trials establish survival gains; selected animal studies establish finite-follow-up eradication; neither closes the translation gap.

**Stop here.** No fitting, candidate ranking by new analysis, simulations, code implementation, laboratory work or outreach follows without Scott's next instruction.

## Source key

Full author lists, indexed publication dates, verified identifiers, inspection depth and supported claims are preserved in `SOURCES.json`. Online and issue dates sometimes differ.

- **S01:** Sun et al. *ResMap: A community resource for systematic mapping of therapy-persistent residual cancer cell dependencies across contexts*. Science Advances, 2026. https://doi.org/10.1126/sciadv.aed7476
- **S02:** Criscione et al. *The landscape of therapeutic vulnerabilities in EGFR inhibitor osimertinib drug tolerant persister cells*. npj Precision Oncology, 2022. https://doi.org/10.1038/s41698-022-00337-w
- **S03:** Jänne et al. *Survival with Osimertinib plus Chemotherapy in EGFR-Mutated Advanced NSCLC*. NEJM, online 2025; issue 2026. https://doi.org/10.1056/NEJMoa2510308
- **S04:** Yang et al. *Overall Survival with Amivantamab–Lazertinib in EGFR-Mutated Advanced NSCLC*. NEJM, 2025. https://doi.org/10.1056/NEJMoa2503001
- **S05:** Oren et al. *Cycling cancer persister cells arise from lineages with distinct programs*. Nature, 2021. https://doi.org/10.1038/s41586-021-03796-6
- **S06:** Haderk et al. *Focal adhesion kinase-YAP signaling axis drives drug-tolerant persister cells and residual disease in lung cancer*. Nature Communications, 2024. https://doi.org/10.1038/s41467-024-47423-0
- **S07:** *Genome-wide CRISPR screens identify the YAP/TEAD axis as a driver of persister cells in EGFR mutant lung cancer*. Communications Biology, 2024. https://doi.org/10.1038/s42003-024-06190-w
- **S08:** Wang et al. *Transient IGF-1R inhibition combined with osimertinib eradicates AXL-low expressing EGFR mutated lung cancer*. Nature Communications, 2020. https://doi.org/10.1038/s41467-020-18442-4
- **S09:** Bertino et al. *Phase IB Study of Osimertinib in Combination with Navitoclax in EGFR-mutant NSCLC Following Resistance to Initial EGFR Therapy (ETCTN 9903)*. Clinical Cancer Research, 2021. https://doi.org/10.1158/1078-0432.CCR-20-4084
- **S10:** Hangauer et al. *Drug-tolerant persister cancer cells are vulnerable to GPX4 inhibition*. Nature, 2017. https://doi.org/10.1038/nature24297
- **S11:** Higuchi et al. *FSP1 and histone deacetylases suppress cancer persister cell ferroptosis*. Science Advances, 2026. https://doi.org/10.1126/sciadv.aea8771
- **S12:** Baldacci et al. *Eradicating Drug-tolerant Persister Cells in EGFR-Mutated Non-Small Cell Lung Cancer by Targeting TROP2 with CAR-T Cellular Therapy*. Cancer Discovery, 2025. https://doi.org/10.1158/2159-8290.CD-24-1515
- **S13:** Brea et al. *Systematic Engineering of TROP2-Targeted CAR T-Cell Therapy Overcomes Resistance Pathways in Solid Tumors*. Cancer Immunology Research, 2025. https://doi.org/10.1158/2326-6066.CIR-25-0527
- **S14:** Liao et al. *Targeting TROP2 in drug-tolerant persister cells delays EGFR tyrosine kinase inhibitor resistance in non-small-cell lung cancer*. Cancer Cell, 2026. https://doi.org/10.1016/j.ccell.2026.05.020
- **S15:** Liu et al. *Transfer of Damaged Mitochondria from Cancer Cells to Cancer-Associated Fibroblasts Promotes Tyrosine Kinase Inhibitor Tolerance in EGFR-Mutant Lung Cancer*. Cancer Research, 2026 issue. https://doi.org/10.1158/0008-5472.CAN-25-0433
- **S16:** Isozaki et al. *Therapy-induced APOBEC3A drives evolution of persistent cancer cells*. Nature, 2023. https://doi.org/10.1038/s41586-023-06303-1
- **S17:** *Osimertinib plus consolidative radiotherapy for advanced EGFR mutant non–small cell lung cancer: a multicentre, single-arm, phase 2 trial*. eClinicalMedicine, 2025. https://doi.org/10.1016/j.eclinm.2025.103435
- **S18:** Herbst et al. *Molecular residual disease analysis of adjuvant osimertinib in resected EGFR-mutated stage IB–IIIA non-small-cell lung cancer*. Nature Medicine, 2025. https://doi.org/10.1038/s41591-025-03577-y
- **S19:** Maynard et al. *Therapy-Induced Evolution of Human Lung Cancer Revealed by Single-Cell RNA Sequencing*. Cell, 2020. https://doi.org/10.1016/j.cell.2020.07.017
