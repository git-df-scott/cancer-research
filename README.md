# Follicular lymphoma and T-cell-engager scheduling: does tissue architecture change the answer?

A spatial agent-based model asking one narrow question, under active adversarial testing.

**Status as of 2026-09-07: UNRESOLVED. Phase 2 in progress.** The model has been recalibrated
against external experimental data and refrozen; the architecture hypothesis is currently
**untested**, pending a positive control. Nothing here is a validated biological claim. Read
`FINDINGS.md` for the live status, `PHASE2_PREREG.md` for what was fixed in advance, and
`PRIOR_ART.md` before believing any novelty claim.

Experiment L1 is preserved unaltered and reproduces bit-identically. It was **uninformative**,
not negative: its positive control failed, for two independent reasons that Phase 2 identified
and fixed — a structurally mis-specified exhaustion submodel, and a schedule family that never
contained the effect being reproduced.

## The question

Obertopp, Froid, Pilon-Thomas & Basanta (bioRxiv 2025, doi 10.1101/2025.11.17.688873) modelled
treatment-free intervals (TFIs) for a CD19xCD3 T-cell engager in B-cell acute lymphoblastic
leukaemia and concluded that short 2-3 day intervals beat the clinical 7-day interval, and that
continuous dosing is worst. Their model seeds tumour cells randomly at 50% occupancy with T cells
already interspersed - their own stated limitation, and a fair model of leukaemia, where an effector
and a blast are already neighbours.

Follicular lymphoma is not like that. Imaging mass cytometry of paired biopsies found that
"peri-follicular regions represented a barrier for immune infiltration into the follicles", with
malignant cells inside follicles "separated spatially from the attack by CD8+ T cells"
(J Hematol Oncol 2022, PMC9396877).

**Does that architectural difference invalidate the leukaemia scheduling conclusion?**

## The candidate causal chain under test

```
follicular architecture
  -> smaller fraction of T cells simultaneously in contact with malignant B cells
  -> lower cumulative contact dwell time per effector
  -> slower exhaustion
  -> less benefit from treatment-free recovery intervals
  -> different optimal schedule than in dispersed / leukaemia-like disease
```

Links 1 and 2 are confirmed in-model and are **not novel** - geometry changing contact frequency is
expected and covered by prior literature. The question is whether the chain reaches the last step
strongly enough to be a falsifiable prediction.

## What is NOT claimed

- No cure claim. No clinical recommendation.
- No new biological mechanism. Every individual link in the chain has prior art; see the claim
  matrix in `PRIOR_ART.md`.
- Simulated days under stated rate calibrations are not clinical predictions.

## Layout

| Path | What it is |
|---|---|
| `lymphoid.py` | the model: 2D lattice, malignant B cells, motile exhaustible T cells, engager schedule |
| `exp_schedule.py` | experiment L1, pre-registered in its docstring (P1-P5) |
| `exp_attacks.py` | experiment L2 against the legacy model: exhaustion knockout, dose-matched controls, trafficking knockout |
| `calib_philipp.py` | Phase B: fits the exhaustion clock to seven measured points from Philipp 2022, by grid search. Reports identifiability, not one best value |
| `calib_abm_check.py` | Phase B stage 2: rebuilds Philipp's assay on the lattice and checks the fit transfers to the real `Lymphoid.step` path |
| `exp_influx_check.py` | Phase C: tests T-cell influx regimes against patient data not used in the fit, and excludes the incompatible ones |
| `PHASE2_PREREG.md` | Phase D: the frozen Phase 2 pre-registration, committed before any Phase 2 result existed |
| `schedules.py` | both schedule families, including the repeating family L1 never tested |
| `exp_poscontrol.py` | Phase E: the positive control, dispersed architecture only |
| `exp_arch2.py` | Phases F and G: architecture comparison and the trafficking attack |
| `exp_attacks2.py` | Phase H: exhaustion knockout, dose matching, exhaustion distribution, geometry continuum |
| `verify_nbrsum.py` | equivalence evidence for the optimised neighbour-sum against the original |
| `TRAFFICKING_PREREG.md` | Phase 5 pre-registration, written before any alternative was run |
| `analyze_L1_full.py` | reports every pre-registered endpoint regardless of outcome |
| `UNDERSTANDING_NHL.md` | the biology, sourced, including what the first attempt got wrong |
| `PRIOR_ART.md` | prior-art audit and claim matrix |
| `FINDINGS.md` | live results and current classification |
| `results/` | raw outputs, all seeds preserved |
| `superseded/` | the failed first attempt, kept deliberately |

## Calibration and its sources

| Quantity | Value | Source |
|---|---|---|
| T-cell speed, lymph node cortex | ~11 um/min | Miller/Cahalan, Science 2002 |
| CTL killing capacity in vivo | 2-16 targets/CTL/day | Halle et al., Immunity 2016, PMID 26872694 |
| Engage-kill-detach cycle | ~25 min | Cazaux et al., J Exp Med 2019 |
| Exhaustion accrual on remaining function | k = 0.159 /day, identified | fitted to Philipp et al., Blood 2022, PMID 35878001 |
| Exhaustion under continuous engager | specific lysis 88.4% (d7) -> 34.9% (d14) -> 8.6% (d28) | same |
| Recovery on a treatment-free interval | d14 93.4% vs 34.9%; **d28 58.7% vs 8.6%** | same |
| Scale of remaining function (non-saturating readout) | granzyme B MFI ratio, d14 144.5/451.8, d28 45.5/196.1 | same |
| Fraction of exhaustion a 7-day break reverses | 4.3%-55%, not further identifiable | same, carried as a sensitivity axis |
| Population-level functional collapse in patients | ex vivo lysis 73.1% -> 17.4% by d14 | same, Figure 1B; used only to exclude influx regimes |
| Neoplastic follicle diameter | 757 um (577-930) | intestinal FL morphometry |
| Malignant B-cell cycle | ~2 days | standard human tumour cell-cycle time |

One lattice site = one cell = 10 um. One step = 1 minute (runs use dt = 5 min with sub-stepped
motility; validated against dt = 1, agreement within 13% on all reported quantities, and the
discretisation is common-mode across arms).

## Reproduce

```bash
python -m venv .venv && .venv/bin/pip install numpy scipy matplotlib

# Phase 1, preserved. Reproduces bit-identically.
.venv/bin/python exp_schedule.py       # experiment L1
.venv/bin/python analyze_L1_full.py    # all pre-registered endpoints

# Phase 2, in order. Each step depends on the one before it.
.venv/bin/python calib_philipp.py      # B   fit the exhaustion clock to external data
.venv/bin/python calib_abm_check.py    # B   verify the fit in the real code path
.venv/bin/python exp_influx_check.py   # C   exclude inadmissible influx regimes
.venv/bin/python verify_nbrsum.py      #     optimisation equivalence evidence
.venv/bin/python exp_poscontrol.py     # E   positive control. Everything after this is
                                       #     gated on it passing.
.venv/bin/python exp_arch2.py F        # F   architecture comparison
.venv/bin/python exp_arch2.py G        # G   trafficking attack
.venv/bin/python exp_attacks2.py K     # H   exhaustion knockout
.venv/bin/python exp_attacks2.py D     # H   dose matching
.venv/bin/python exp_attacks2.py P     # H   exhaustion distribution
.venv/bin/python exp_attacks2.py C     # H   geometry continuum
```

## Known deviations, recorded rather than hidden

- **2026-09-07 out-of-memory stall.** The first L1 launch deadlocked at 125/150 runs when system
  swap filled. The 125 completed runs are preserved verbatim in `results/expL1_salvage125.json`;
  the remaining 25 cells were re-run by `rerun_missing.py` with identical code and seeds and merged
  into `results/expL1.json`. Nothing was discarded or re-rolled.
- **Floor effect on the primary endpoint.** Under continuous dosing both architectures reach
  near-zero burden by day 42 (medians 4 and 7 cells of 5525), so the pre-registered day-42 endpoint
  may lack dynamic range. Day 28 is reported alongside it and is clearly labelled as secondary.
