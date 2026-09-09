# What has already been done, and whether this is a rediscovery

Assessed 2026-09-07 from three independent prior-art searches plus my own. The honest summary is at
the bottom; read that first if you read nothing else.

## The claim being tested for novelty

> Tumour **architecture** (densely packed follicles vs dispersed cells) controls what fraction of
> effector T cells are in antigen contact at any instant. Because exhaustion accrues with contact
> **dwell time**, dense packing protects the effector pool from exhaustion. Therefore the optimal
> **treatment-free interval** for a T-cell engager is architecture-dependent, and the leukaemia
> answer does not transfer to follicular lymphoma.

## Every link in that chain is already published

This is the part that matters. I am not discovering a new mechanism.

| Link in the chain | Already established by |
|---|---|
| T cells are restricted to the rim of dense tumour tissue and fail to penetrate the core | Long known experimentally in spheroids; modelled directly in a 3D agent-based model, bioRxiv 2025.10.21.683801 |
| Spatial configuration sets the CTL–target encounter and contact rate | "Tissue Dimensionality Influences the Functional Response of CTL-Mediated Killing", cellular Potts model, PMC5225319 |
| Exhaustion can be driven by cumulative antigen exposure inside a spatial ABM with density structure | TIED-ABM multiscale exhaustion framework, PLOS Comput Biol, PMC13541251 |
| Treatment-free intervals reverse engager-induced T-cell exhaustion | Philipp et al., Blood 2022, PMID 35878001 (experimental); Weber et al., Science 2021, PMID 33795428 (CAR-T rest) |
| Treatment-free intervals can be optimised in a lattice ABM of a T-cell engager | Obertopp, Froid, Pilon-Thomas & Basanta, bioRxiv 2025, doi 10.1101/2025.11.17.688873 |
| Tissue architecture couples to exhaustion and predicts immunotherapy outcome | "Tumor-immune trajectory context connects static tissue architecture to clinical outcomes", PhysiCell ABM, bioRxiv 2026.03.26.714521 |
| Physical and stoichiometric barriers limit bispecific engager efficacy in tissue | Intravital imaging, PMC9177795 |
| Tumour bulk predicts bispecific and CAR-T failure clinically | Metabolic tumour volume >504.7 cm3 gives median PFS 1 vs 3 months, HR 5.6 |

## What the searches did not find

All three searches, run independently against different angles, agreed on the same negative result.

- **No agent-based model varies packing density or clustering as an independent, manipulated
  architecture variable** in the context of a T-cell engager or CAR-T. Where density appears it is
  either a fixed calibration constant or an *emergent output* (in TIED-ABM the causality runs
  backwards: exhaustion drives the tumour to consolidate into dense nodules).
- **No model makes exhaustion a function of contact dwell time** rather than cumulative antigen dose
  or serial-killing-capacity depletion. The CARCADE CAR-T model has an exhaustion state and lists
  precisely this as a limitation.
- **No model derives an architecture-conditional dosing interval.** The schedule-optimisation work is
  architecture-blind; the architecture work is schedule-blind.
- **No agent-based model of follicular lymphoma exists at all.** The germinal-centre ABM literature is
  mature and has been extended once to DLBCL genetics (Front Syst Biol 2022), but never to therapy.
- **Nobody states the distinction explicitly** that the rationale for treatment holidays under
  immunotherapy is *immunological* (effector recovery) rather than *evolutionary* (competitive
  suppression of resistant clones). The adaptive-therapy literature argues the second; the
  exhaustion literature demonstrates the first; no one has said they are different mechanisms that
  transfer to different diseases.
- The PhysiCell/CompuCell3D angle on this specific question returned nothing at all.

## Closest existing work, which must be cited and built on

1. **Obertopp/Basanta bioRxiv 2025** — closest on the schedule half. Lattice ABM, blinatumomab
   treatment-free intervals in B-ALL, concludes short 2-3 day intervals beat the clinical 7-day one
   and continuous is worst. Random seeding at 50% occupancy, T cells interspersed from t=0, in vitro
   calibration. Their own listed limitation is "random rather than clustered cell seeding".
2. **TIED-ABM, PLOS Comput Biol** — closest on the exhaustion half. Spatial voxel ABM, density
   thresholds, delay-integral exhaustion from cumulative antigen exposure. No drug, no schedule.
3. **PhysiCell trajectory paper, bioRxiv 2026** — closest conceptually. Dispersed vs dense contiguous
   tumour states, CD8 exhaustion rate as dominant outcome driver. But checkpoint blockade in breast
   cancer, exhaustion as a rate parameter, density not independently varied, no scheduling question.
   **I have not yet read this in full and must do so before claiming novelty.**

## Honest verdict

**This is not a rediscovery, but neither is it a new mechanism.** Every individual link in the causal
chain is published. What is unpublished is the conjunction, and the specific consequence that falls
out of it: that a dosing-schedule conclusion derived in a well-mixed or leukaemic setting should not
be carried over to follicular lymphoma, because the architecture changes how fast the effector pool
degrades.

That is a modest contribution. It is worth making because it is checkable and because it bears on a
live clinical question with no randomised answer: glofitamab is given for a fixed 12 cycles and
stopped, epcoritamab is given continuously until progression, and no head-to-head trial exists.

**The two objections that could still sink it**, and which the control experiment (`exp_controls.py`,
pre-registered) is designed to answer:

- *It is trivially true.* The contact-fraction difference between a packed disc and a dispersed
  population is close to a surface-to-volume argument, and anyone could predict its sign without a
  model. The defence is not that the sign is surprising but that the **magnitude and its consequence
  for schedule ranking** are not predictable by inspection. If the schedule ranking does not change,
  there is no result.
- *It is a dose effect, not an exhaustion effect.* Treatment-free-interval arms deliver less total
  drug. Dose-matched controls are required, and if a TFI arm does not beat its own dose-matched
  continuous twin, there is no timing effect to report at all.

A third possibility, which would change the claim rather than kill it: the mean exhaustion may be low
in a follicle while the front cells are fully exhausted and continuously replaced from a fresh
reservoir. That is a **different and more interesting** mechanism, and if the control experiment shows
it, the write-up must say so rather than keep the simpler story.

---

# Phase 8: claim matrix, after reading the closest works in full

Both of the two closest papers were fetched and read in full on 2026-09-07, discharging the debt
recorded above.

## What the two closest works actually do

**"Tumor-immune trajectory context connects static tissue architecture to clinical outcomes"**
(PhysiCell ABM, bioRxiv 2026.03.26.714521, PMC13060115, PubMed 41959335). Read in full.

- Architecture **emerges** from cell-level rules; it is not manipulated as an independent input.
  Initial conditions are ~1,100 cells with "configurations constrained to maintain comparable
  spatial relationships".
- CD8 exhaustion is a **fixed rate parameter** sampled by Latin hypercube, explicitly not
  contact-dependent. The paper's own headline is that "the CD8+ T cell exhaustion rate exhibited the
  highest variability across TME states".
- Therapy is checkpoint blockade only, implemented as reducing the exhaustion-rate parameter.
- **No dosing schedules or treatment intervals are varied.** Intervention timing is examined against
  simulation time, independent of spatial configuration.
- It does **not** claim optimal timing depends on architecture.

**TIED / TIED-ABM multiscale exhaustion framework** (PLOS Comput Biol, PMC13541251). Read in full.

- Exhaustion is driven by **cumulative antigen exposure through a time-delayed integral**
  (a Hill function of the integral of tumour antigen over a window), not by contact dwell time.
- Tumour density is an **emergent output** that modulates behaviour probabilities, not an imposed
  architecture variable.
- Therapy is anti-PD-1 and anti-CTLA-4, implemented as modulating exhaustion and activation rate
  parameters. No schedule is optimised.

## Claim matrix

| Claim | Already known | Closest source | What this project adds |
|---|---|---|---|
| Dense tumour tissue restricts T cells to the rim | **Yes**, experimentally and in models | Spheroid imaging; 3D ABM bioRxiv 2025.10.21.683801 | Nothing |
| Spatial configuration sets CTL–target contact rate | **Yes** | Cellular Potts, tissue dimensionality, PMC5225319 | Nothing |
| Exhaustion accrues from antigen exposure in a spatial ABM | **Yes** | TIED-ABM, PMC13541251 | A *dwell-time* rather than cumulative-integral formulation. Minor. |
| Architecture couples to exhaustion and to immunotherapy outcome | **Yes** | PhysiCell trajectory paper, PMC13060115 | Nothing on outcome |
| Treatment-free intervals reverse engager exhaustion | **Yes**, experimentally | Philipp, Blood 2022, PMID 35878001 | Nothing |
| Treatment-free intervals can be optimised in a lattice ABM of an engager | **Yes** | Obertopp/Basanta, bioRxiv 2025 | Nothing |
| Architecture manipulated as an **independent input** in an engager model | **No** | — | This project varies it directly |
| Exhaustion as a function of **contact dwell time** | **No** | CARCADE lists this as a limitation | This project implements it |
| **Architecture-conditional optimal dosing interval** | **No** | — | This project tests it |
| **Non-transfer of a leukaemia scheduling result to follicular architecture** | **No** | — | This project tests it |

## Novelty verdict after full reading

The bottom four rows are genuinely unclaimed. The top six are not, and I add nothing to them.

That is a narrow strip of novelty resting entirely on a conjunction, and it is only worth anything
if the conjunction produces an effect that survives the attacks in Phases 3 to 7. As of writing, the
Phase 5 trafficking test has already produced a result that threatens it badly (see
`docs/findings/FINDINGS.md`), so novelty is not the binding constraint here. Validity is.
