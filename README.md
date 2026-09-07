# Follicular lymphoma and T-cell-engager scheduling: does tissue architecture change the answer?

A spatial agent-based model asking one narrow question, under active adversarial testing.

**Status as of 2026-09-07: UNRESOLVED, and currently leaning towards a negative result.**
Nothing here is a validated biological claim. Read `FINDINGS.md` for the live status and
`PRIOR_ART.md` before believing any novelty claim.

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
| `exp_attacks.py` | experiment L2: exhaustion knockout, dose-matched controls, trafficking knockout, exhaustion distribution |
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
| Exhaustion under continuous engager | specific lysis 88.4% (d7) -> 8.6% (d28) | Philipp et al., Blood 2022, PMID 35878001 |
| Recovery on a treatment-free interval | d14 lysis 93.4% vs 34.9% continuous | same |
| Neoplastic follicle diameter | 757 um (577-930) | intestinal FL morphometry |
| Malignant B-cell cycle | ~2 days | standard human tumour cell-cycle time |

One lattice site = one cell = 10 um. One step = 1 minute (runs use dt = 5 min with sub-stepped
motility; validated against dt = 1, agreement within 13% on all reported quantities, and the
discretisation is common-mode across arms).

## Reproduce

```bash
python -m venv .venv && .venv/bin/pip install numpy scipy matplotlib
.venv/bin/python exp_schedule.py     # experiment L1
.venv/bin/python analyze_L1_full.py  # all pre-registered endpoints
.venv/bin/python exp_attacks.py      # the attack battery
```

## Known deviations, recorded rather than hidden

- **2026-09-07 out-of-memory stall.** The first L1 launch deadlocked at 125/150 runs when system
  swap filled. The 125 completed runs are preserved verbatim in `results/expL1_salvage125.json`;
  the remaining 25 cells were re-run by `rerun_missing.py` with identical code and seeds and merged
  into `results/expL1.json`. Nothing was discarded or re-rolled.
- **Floor effect on the primary endpoint.** Under continuous dosing both architectures reach
  near-zero burden by day 42 (medians 4 and 7 cells of 5525), so the pre-registered day-42 endpoint
  may lack dynamic range. Day 28 is reported alongside it and is clearly labelled as secondary.
