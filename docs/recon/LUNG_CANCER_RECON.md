# Lung cancer recon: where the bottleneck actually is, and which leads this project could move

Recon pass, 2026-09-08. Literature checked live against Europe PMC, PubMed, PMC full text and
trial reporting. Every factual claim below carries a source. Nothing here is a validated result and
nothing here has been modelled yet.

## Framing, stated once and honestly

There is no lead here that produces a cure for lung cancer, and anyone offering one is selling
something. What the recon is actually for is finding the places where (a) a specific, named
mechanism is blocking durable disease control, (b) the mechanism is quantified well enough to be
falsifiable, and (c) this project's existing tooling — a spatial agent-based model of a CD3 T-cell
engager with motile, exhaustible T cells — can say something a wet lab has not already said.

Judged on (c), most of the lung landscape is out of reach. One lead is not, and it is unusually
well matched.

---

## 1. The landscape, and where the deaths come from

| Setting | Status | Is there a modelable bottleneck? |
|---|---|---|
| **ES-SCLC** | Near-uniformly fatal. Chemo works then fails fast. Tarlatamab is the first real second-line win in decades | **Yes — see Lead 1** |
| **EGFR-mutant NSCLC** | Osimertinib works, then drug-tolerant persisters seed resistance | Partly — but the bottleneck is cell-intrinsic, not spatial |
| **KRAS-mutant NSCLC** | G12C inhibitors underwhelm; G12D and pan-RAS now delivering | No — bottleneck is chemistry and adaptive signalling |
| **STK11/KEAP1 co-mutant** | Reliably immunotherapy-resistant | Weakly — mechanism is metabolic/myeloid, not geometric |
| **Early-stage / interception** | Where cures actually happen | No — this project has no tooling for it |

### SCLC is the honest target

SCLC is an "aggressive pulmonary neuroendocrine malignancy featured by a cold tumor immune
microenvironment, limited benefit from immunotherapy, and poor survival"
([Cell Discovery 2024, PMC11375181](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11375181/)).
Then tarlatamab, a DLL3xCD3 bispecific engager, changed the second-line picture:

- **DeLLphi-304 (phase 3):** median OS 13.6 vs 8.3 months against chemotherapy, HR for PFS 0.71
  ([NEJM](https://www.nejm.org/doi/abs/10.1056/NEJMoa2502099);
  [OncLive](https://www.onclive.com/view/dellphi-304-data-reinforce-role-of-tarlatamab-as-second-line-soc-in-small-cell-lung-cancer))
- **DeLLphi-303 (phase 1b, first-line maintenance + PD-L1 inhibitor):** median OS 25.3 months with
  atezolizumab, not yet reached with durvalumab, 12-month OS 82%
  ([Lancet Oncol 2025, PMID 40934933](https://pubmed.ncbi.nlm.nih.gov/40934933/))

So: the drug class works in a solid tumour, the effect size is large, and the field is now
scaling it into first line (DeLLphi-305). **That is exactly when schedule questions become
worth asking, and exactly when nobody is asking them.**

---

## 2. Lead 1 — the strongest, and it is a specific mechanistic claim

### The observation

Two resistance axes to tarlatamab are now documented in patients, and one of them is exhaustion.

**Antigen-side escape.** Transcription factor subtype governs response: 80.0% of ASCL1-subtype
patients had clinical benefit versus 20.0% of NEUROD1-subtype and 0% of POU2F3. Median PFS 9.7
months (SCLC-A) versus 2.0 months (SCLC-N and SCLC-P). DLL3 is a direct transcriptional target of
ASCL1. At progression, one patient showed selection for a NEUROD1-high state with DLL3
downregulation, and copy-number analysis showed the resistant clone "did not evolve from the clone
that predominated at treatment initiation... it may have been present before therapy and become
dominant due to the selective pressure of treatment"
([PMC13082099](https://pmc.ncbi.nlm.nih.gov/articles/PMC13082099/), PMID 41993528).

**Execution-side escape.** A separate resistant lesion showed a high proportion of FOXP3+
regulatory cells and T cells with an exhausted phenotype (CTLA4, TIGIT, LAG3). On rechallenge,
tumours that retained DLL3 showed FOXP3+ Treg enrichment and increased TIM3 and LAG3 on CD8+ T
cells — i.e. **failure with the antigen still present**. Independently, a JTO brief report found
resistance to tarlatamab plus anti-PD-1 driven by low DLL3 *and* T-cell exhaustion
([PMID 41903703](https://pubmed.ncbi.nlm.nih.gov/41903703/)).

### The gap nobody has closed

A 2026 review of DLL3-directed redirection states plainly that preclinical evidence supports
treatment-free intervals ameliorating engager-induced exhaustion, but that **this "has not yet been
formally evaluated in tarlatamab-treated SCLC patients"**
([PMC13457367](https://pmc.ncbi.nlm.nih.gov/articles/PMC13457367/)).

### The mechanistic claim, and why it is not obvious

Tarlatamab is a **half-life extended** BiTE. Terminal half-life is **5.8 days** (DeLLphi-300) to a
model-estimated median of **11.2 days** (population PK, 420 patients), dosed **10 mg every 2 weeks**
after 1 mg step-up priming
([Clin Pharmacokinet, PMID 39589690](https://pubmed.ncbi.nlm.nih.gov/39589690/);
[PMID 40261494](https://pubmed.ncbi.nlm.nih.gov/40261494/)).

The entire treatment-free-interval literature descends from **blinatumomab, half-life ~2 hours,
continuous infusion**. For blinatumomab, "stop dosing" means drug concentration goes to zero within
hours, and T cells genuinely disengage.

**With an 11-day half-life dosed every 14 days, there is no treatment-free interval.** Concentration
never approaches zero; the schedule is a sawtooth on top of continuous target engagement. The
half-life extension that made tarlatamab clinically dosable may have simultaneously **removed the
disengagement window that limits exhaustion** — and clinical exhaustion-mediated resistance is now
observed. Nobody has connected these two facts.

### The falsifiable hypothesis

> In a half-life-extended CD3 engager, the PK-imposed floor on trough target engagement, not the
> nominal dosing interval, determines cumulative effector dwell time and therefore exhaustion. A
> Q2W schedule of an 11-day-half-life agent is pharmacodynamically closer to continuous dosing
> than to an intermittent one, and its exhaustion burden should track the trough, not the interval.

This is testable in-model and **has a hard external decision rule**: if simulated trough
concentration is high enough that engaged fraction never meaningfully falls between doses, the
hypothesis predicts exhaustion accrues as under continuous dosing. If engaged fraction does fall
substantially, the hypothesis is dead and should be reported dead.

It also generates a **dose-fractionation** counter-proposal that is clinically reachable: same total
exposure, lower peak, engineered trough — as opposed to the current schedule, which was selected for
**CRS safety and dosing convenience, not for exhaustion**. Step-up dosing exists "to blunt the
initial activation peak," and CRS is "concentrated around the priming and first full doses"
(PMC13457367). Nothing in the schedule's derivation optimised for week-12 T-cell function.

### Prior art risk: moderate, and specifically checked

- Population PK models exist (NONMEM, two-compartment, linear elimination) but are **PK only** — no
  exhaustion, no spatial term, no schedule optimisation (PMID 40261494).
- QSP models of CD3 bispecifics in solid tumour exist ([PMC7293198](https://pmc.ncbi.nlm.nih.gov/articles/PMC7293198/),
  P-cadherin LP DART) but are non-spatial and not tarlatamab.
- The ABM that motivates this project (Obertopp/Basanta) explicitly used **AMG 562's 210-hour
  half-life as a stand-in** and flagged the PK simplification as a limitation.

No one has combined half-life-extended PK, exhaustion, and spatial structure for tarlatamab. That
is a real gap, and it is narrow enough to be defensible.

---

## 3. Lead 2 — SCLC spatial architecture maps onto the follicle model better than follicular lymphoma did

This matters because it decides whether the existing code is reusable.

Ultra-high-plex Digital Spatial Profiling of **132 tissue microarray cores from 44 treatment-naive
limited-stage SCLC tumours** resolved three pathologically defined compartments: **pan-CK+ tumour
nest**, **CD45/CD3+ tumour stroma**, and para-tumour
([Cell Discovery 2024, PMC11375181](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11375181/)).
T cells are in the stroma; the nests are cold. Fibroblast domains are "frequently observed
interposing between and segmenting distinct tumour domains"
([Cell Rep Med 2026](https://www.cell.com/cell-reports-medicine/fulltext/S2666-3791(26)00158-8)).

That is the **same geometry** as the follicular lymphoma question — an effector compartment
separated from a malignant compartment by a barrier — with three advantages over the FL case:

1. **The trafficking knockout threat is weaker.** The FL model's biggest vulnerability was that
   T cells at 11 µm/min in dense lymph-node cortex plausibly squeeze past B cells, and allowing that
   collapsed the architecture effect (`docs/calibration/TRAFFICKING_PREREG.md`, `docs/findings/FINDINGS.md`). A desmoplastic
   fibroblast/matrix barrier around an epithelial tumour nest is a **physically different and much
   better-evidenced** exclusion mechanism than same-size lymphocytes jostling in a follicle.
2. **MHC-I loss is a confound the engager bypasses.** SCLC-enriched regions "consistently express
   lower levels of HLA class I genes and B2M," contributing to T-cell exclusion (PMC11375181). A
   CD3 engager does not care — it bridges independently of MHC. That makes the engager arm a
   *cleaner* probe of geometry than checkpoint blockade would be.
3. **The clinical schedule is fixed and public**, so the model has something concrete to be
   right or wrong about.

**Caveat, and it is a real one:** the compartment data above is TMA-core-level and I could not
retrieve nest diameters or T-cell-to-tumour distances — nature.com is 403-blocked from this
environment, consistent with the constraints already recorded in `docs/plan/HANDOFF.md`. **Calibrating nest
geometry requires a source not yet obtained.** Until it is, any SCLC geometry is a guess, and the
FL project's own history says guessed geometry is how you get a model artefact.

---

## 4. Leads deliberately rejected, with reasons

| Lead | Why rejected |
|---|---|
| **Osimertinib drug-tolerant persisters** | Genuinely important — MEK, AURKB, BRD4, TEAD vulnerabilities are reproducible ([PMC9794691](https://pmc.ncbi.nlm.nih.gov/articles/PMC9794691/)), EP2 is new ([PMID 42647987](https://pubmed.ncbi.nlm.nih.gov/42647987/)). But the bottleneck is cell-intrinsic state-switching, not spatial immune dynamics. No comparative advantage here, and the field is crowded. |
| **KRAS G12D / pan-RAS** | Moving fast without us — zoldonrasib has FDA Breakthrough designation (Jan 2026) in previously-treated G12D NSCLC ([AACR](https://www.aacr.org/about-the-aacr/newsroom/news-releases/investigational-krason-inhibitor-zoldonrasib-showed-effective-and-durable-responses-in-patients-with-advanced-g12d-mutated-lung-cancer/)). Bottleneck is medicinal chemistry. Nothing an ABM adds. |
| **STK11/KEAP1 resistance** | Real and unsolved; ATR inhibition is an interesting angle ([Cancer Cell 2025](https://www.cell.com/cancer-cell/fulltext/S1535-6108(25)00261-2)). But the mechanism is metabolic/myeloid, requiring model machinery this project does not have. |
| **Early detection / interception** | Where cures actually come from, and entirely outside this toolchain. Flagging it as the honest highest-impact area, not as our lead. |

---

## 5. Portability audit — what survives the move from `lymphoid.py`

This is the part that decides whether the lead is cheap or expensive.

**Transfers directly:** 2D lattice with 10 µm sites, motile T cells with volume exclusion, engage-
kill-detach cycle, exhaustion state variable, compartmentalised seeding, the whole pre-registration
and attack-battery discipline, `exp_attacks.py` structure (knockouts, dose-matched controls).

**Must be rebuilt, and this is not optional:**

1. **A PK compartment.** The FL model treats the drug as a binary on/off. That is correct for
   blinatumomab and **wrong for tarlatamab**, and the wrongness is the entire point of Lead 1.
   Needs two-compartment linear elimination, CL 0.649 L/day, Vc 3.44 L for a 73 kg patient
   (PMID 40261494), driving an occupancy term rather than a switch.
2. **Antigen heterogeneity.** DLL3 positivity is defined as expression on ≥25% of tumour cells;
   85–94% of SCLC is positive by that rule, but only 46.9% is high by H-score
   ([Frontiers Immunol 2025](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2025.1592291/full)).
   A model with uniform antigen cannot represent the dominant clinical escape route, which is
   selection of a pre-existing DLL3-low clone. **Adding a heritable per-cell antigen level is the
   single highest-value structural change**, and it makes the model able to address antigen-side and
   execution-side escape in one framework — which is what the clinical papers say is needed.
3. **The exhaustion recalibration that is already the FL project's blocking task.** It does not go
   away by changing disease. See below.

---

## 6. The blocker that follows us, and a correction to `docs/findings/FINDINGS.md`

`docs/plan/HANDOFF.md` names the first task: recalibrate exhaustion, then verify the positive control passes.
**That blocker is inherited by any lung work using this model.** I read the reference model's
published version ([PMC12667981](https://pmc.ncbi.nlm.nih.gov/articles/PMC12667981/)) and it
sharpens the diagnosis considerably.

**The reference model's actual reported rankings:**

- Day 28: TFI_2 ≈ TFI_3 > TFI_7 > CONT
- Day 42: TFI_2/TFI_3 > CONT > TFI_7

**Two things follow, and they point in opposite directions.**

*The positive-control failure is real and is not a floor effect.* At day 28 this project's dispersed
arm gives median burden CONT 250.5 vs TFI2 589.5 — continuous better by 2.4×, where the reference
has TFI2 better. Day 28 is far from the floor, so this cannot be explained away by the day-42
floor effect noted in `README.md`. `docs/findings/FINDINGS.md` is right that the control failed.

*But the model is not uniformly wrong.* `docs/findings/FINDINGS.md` states the dispersed arm "produced the
opposite ranking." That is imprecise. The reference **also** reports CONT beating TFI_7 at day 42,
and this project reproduces that (7.0 vs 75.5). The failure is specific to the
**short-TFI-versus-continuous** comparison, not to the whole ordering. Worth stating precisely,
because it narrows what recalibration has to fix.

**The root cause named in `docs/plan/HANDOFF.md` is now externally confirmed.** The handoff hypothesised that
T-cell influx dilutes the exhaustion pool. The reference model **has no influx or recruitment at
all** — fixed population, 2,400 T cells at 1:4 E:T, changing only by proliferation, death and state
transitions — and it reaches **over 80% exhaustion by day 16**. This project's model has influx
(200 → ~860 T cells) and reaches **12% at day 28**. The discrepancy is structural, it is
identified, and it is fixable. That is a substantially better position than "unexplained
calibration failure."

I have not edited `docs/findings/FINDINGS.md`'s classification. The correction above is recorded here and should
be merged into `docs/findings/FINDINGS.md` by whoever next touches the FL analysis.

---

## 7. Recommended sequence, with decision gates

1. **Fix exhaustion calibration against Philipp's external curve** (unchanged from `docs/plan/HANDOFF.md`).
   Prime suspect now confirmed: influx dilution. Gate: reproduce ~90% functional loss at day 28
   continuous, and the day-14 recovery comparison.
2. **Re-run the FL positive control.** Gate: does short-TFI beat continuous in the dispersed arm?
   If no, stop and report. Do not carry a broken model into a second disease.
3. **Only if 2 passes:** build the PK layer and ask the Lead 1 question, which does not require the
   FL architecture result to have succeeded — it is a schedule question, not a geometry question.
4. **Obtain SCLC nest geometry from a retrievable source** before attempting any spatial SCLC claim.
   Without it, Lead 2 is not yet runnable.
5. Add heritable per-cell antigen level, then ask whether schedule modulates selection for
   DLL3-low clones. This is the highest-value question in the whole recon and also the least
   constrained by existing data — treat with corresponding suspicion.

## What is NOT claimed

- No cure claim, no clinical recommendation, no proposed change to tarlatamab dosing.
- No new biological mechanism. Exhaustion, antigen escape, lineage plasticity and PK are all
  established separately; the only candidate novelty is the interaction between half-life-extended
  PK and exhaustion accrual, and that is a hypothesis, not a result.
- No modelling has been performed for lung cancer. This document is recon only.
- Lead 2's geometry is uncalibrated and must not be simulated until it is not.
