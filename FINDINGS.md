# Findings — live status

**Status 2026-09-07 (Phase 2): the model has been recalibrated against external data and refrozen.
The architecture hypothesis remains UNTESTED pending the positive control.**

Experiment L1 is preserved below, unaltered and bit-identical reproducible. Nothing in Phase 2
reinterprets it.

---

## Phase 2 calibration (Phases A–D), completed and committed before any Phase 2 result

### What was wrong with L1, and it was two things, not one

**Cause 1 — the exhaustion clock was mis-specified, not mis-tuned.** The handoff's diagnosis was
that the accrual rate was too slow. It was not: the legacy linear rate drives E→1 over 28 days of
*continuous contact*, which is about right for a linear form. The real problems were that no T cell
in the model ever experiences continuous contact, and that a single fully-reversible exhaustion
state cannot fit the measured curve at *any* parameter values.

Retrieving the full text of Philipp et al. (Blood 2022, PMID 35878001) gave the protocol and four
quantitative readouts the handoff did not have. Two of them force structural changes:

| Recovered fact | Why it forces a structural change |
|---|---|
| TFI arm withdraws the engager on days **7–14 and 21–28**, targets still present | fixes the schedule being calibrated against |
| TFI specific lysis at **day 28 = 58.7%**, not just day-14 93.4% | the TFI arm does **not** return to its day-14 level after its second rest. Recovery is partial and the shortfall grows with cumulative exposure. One reversible state cannot do this → **reversible + durable** |
| Granzyme B MFI ratio, CD8+: d14 CONT 144.5 vs TFI 451.8; d28 45.5 vs 196.1 | a per-cell readout that does **not** saturate the way % specific lysis does. Its CONT/TFI ratio pins the **scale** of remaining function |
| CD2+ proliferation fold change: d14 1.1 vs 4.1; d28 0.06 vs 2.8 | **held out** of the fit as an independent check |

Without the granzyme B constraint the fit drives the assay-saturation constant to 10.0 and declares
every arm functionally dead (f between 0.002 and 0.25) while still reproducing all five lysis
points perfectly. That is a fit to an assay, not to biology, and it was caught only because a
second, non-saturating readout exists in the same figure.

**Cause 2 — L1 tested the wrong schedule family.** The reference being reproduced
(Obertopp/Basanta, PPR1121269) reports that a 7-day interval loses its advantage by day 42, but
that *shorter* intervals consistently win, and that a **Monday-to-Friday** regimen — 5 days on,
2 days off, *every week* — performs comparably. L1 only ever tested a single 28-day cycle with the
last k days off, in which a "2-day TFI" is one 2-day break in 28 days: a 7% dose reduction with a
single recovery opportunity. **L1 never sampled the repeating family at all.** Note that L1 was not
wrong about the 7-day interval; it found TFI-7 loses by day 42, and so does the reference.

### The recalibrated exhaustion submodel

```
in contact, engager present:   dEr/dt = (1-rho) k h(C) (1-E)      reversible
                               dEd/dt =    rho  k h(C) (1-E)      durable
                               dC/dt  = 1                          cumulative engaged exposure
engager absent:                dEr/dt = -Er / tau_r,  Ed unchanged
```

Grid search over 244,800 parameter sets, scored against seven measured points and nothing else.
No schedule ranking, architecture difference or tumour burden was used as a target anywhere.

| Quantity | Status | Value |
|---|---|---|
| `k`, accrual on remaining function | **identified** | **0.159 /day**, at every point of both profiles |
| `rho`, durable share | **not identifiable** | SSE flat from 0.00 to 0.93 |
| `tau_r`, reversible recovery | **not identifiable** | spans three orders of magnitude |
| fraction of exhaustion a 7-day break reverses | **identified as a combination** | **4.3% – 55%**, fit prefers the low end |

The Hill lag did not halve SSE and by the pre-registered rule was **not adopted**. Held-out check:
day-14 proliferation ratio, measured 26.8, model 28.6 — not fitted.

The new clock drives function to 0.33 by day 7 and **0.011 by day 28** under continuous contact.
L1 reached 0.88. That is the gap the handoff identified, now closed against external data.

**Recorded deviation.** The pre-registered RMSE ≤ 5 pp acceptance region is empty; best is 7.69 pp.
The cause is in the data: the three independent day-28 readouts of the *same* CONT/TFI functional
ratio span a 10-fold range (proliferation 2.1, lysis-derived 14.7, granzyme B 23.2), and a
one-functional-axis model predicts a single ratio for all three. The class is therefore defined as
`SSE ≤ 2 × min`. Three representatives spanning the admissible recovery range are carried into
every subsequent run as a sensitivity axis; none is preferred.

### Phase C: the influx rate used in L1 is excluded by external patient data

Anchor independent of the in vitro calibration — Philipp Figure 1B, patients on continuous
blinatumomab: ex vivo specific lysis of **peripheral** T cells fell 73.1% → 17.4% by day 14, i.e. to
**0.238 of baseline**, in a compartment where recruitment, redistribution and engager-driven
expansion are all happening. Criterion fixed before running: admissible only if model pool function
falls to ≤ 0.40 of baseline by day 14 (deliberately generous).

| influx | `rec_low` f(d14) | `rec_mid` | `rec_high` | T pool d28 | verdict |
|---|---|---|---|---|---|
| 0 | 0.122 | 0.264 | 0.370 | ~25 | admissible, all three |
| 2.5e-5 | 0.343 | 0.580 | 0.665 | ~160 | admissible for `rec_low` |
| **1.0e-4 (L1's value)** | **0.533** | **0.828** | **0.901** | ~750 | **excluded, all three** |
| 4.0e-4 | 0.873 | 0.990 | 0.997 | — | excluded |
| 1.6e-3 | — | — | — | — | excluded |

L1's own diagnosis — that influx dilution prevented population-level functional collapse — is
confirmed, and is now excluded on an external measurement rather than by judgement.

**Caveat recorded, not hidden.** Zero influx passes the functional criterion but lets the T-cell
pool fall 200 → ~25 by day 28 through background death, which is its own unrealism. The most
defensible admissible combination is `rec_low` at 2.5e-5: the only one that both reproduces the
functional collapse and holds the pool roughly stable.

**Structural limitation recorded.** Within this model, recruitment and population-level functional
collapse cannot both be represented at a realistic influx rate, because arrivals enter with zero
exhaustion. A non-naive-arrivals variant is listed as a robustness check and is **not** permitted as
a rescue if the positive control fails.

### Phase D

`PHASE2_PREREG.md` was committed before any Phase 2 architecture or positive-control result
existed. The legacy exhaustion path is retained and verified **bit-identical**, so L1 reproduces
exactly.

---

## Experiment L1 (legacy exhaustion model) — preserved unaltered

**Classification: the L1 experiment is UNINFORMATIVE about the hypothesis, because its positive
control failed.** This is not "hypothesis falsified" and it is certainly not support. Phase 2 above
has since identified two independent causes of that failure and fixed both.

## Phase 1 results (experiment L1, 150 runs, 10 seeds, 3 architectures x 5 schedules)

### Pre-registered predictions

| | Prediction | Result | Verdict |
|---|---|---|---|
| P1 | engaged fraction >=3x lower in follicle than dispersed, multi intermediate | 0.098 / 0.135 / 0.463, ratio 4.7 | **HELD** |
| P2 | exhaustion at d28 >=2x lower in follicle | 0.0242 / 0.0352 / 0.1164, ratio 4.8 | **HELD** |
| P3 | some TFI beats continuous in dispersed, none in follicle | **no TFI beat continuous in any architecture** | **FAILED** |
| P4 | TFI benefit ordered dispersed > multi > follicle | no benefit anywhere to order | **FAILED** |
| P5 | TFI benefit correlates with engaged fraction, rho < -0.6 | rho = +0.12, p = 0.71 | **FAILED** |

### Schedule ranking, all three architectures

```
follicle    TFI0(4)  > TFI2(10) > TFI4(18) > TFI7(56)  > TFI14(222)
multi       TFI0(0)  > TFI2(2)  > TFI4(4)  > TFI7(16)  > TFI14(84)
dispersed   TFI0(7)  > TFI2(14) > TFI4(22) > TFI7(76)  > TFI14(414)
```

Continuous dosing wins everywhere, and every interruption is monotonically worse than a shorter one.

## Why this experiment cannot answer the question

**The dispersed arm is the positive control, and it failed.** It was supposed to reproduce the
reference result it is being compared against (Obertopp/Basanta: short TFIs beat continuous in
well-mixed leukaemia-like disease). It produced the opposite ranking. A model that does not
reproduce the effect in the setting where the effect is known cannot be used to argue that the
effect fails to transfer to a different setting.

**The cause is a calibration failure in exhaustion, and it is quantified:**

| Architecture | drug-time delivered | mean exhaustion at d28 | functional loss |
|---|---|---|---|
| dispersed, continuous | 100% | 0.116 | 11.6% |
| dispersed, TFI 7 | 75% | 0.105 | 10.5% |
| follicle, continuous | 100% | 0.024 | 2.4% |
| follicle, TFI 7 | 75% | 0.020 | 2.0% |

Philipp et al. (Blood 2022, PMID 35878001) measured specific lysis falling from 88.4% to 8.6% over
28 days of continuous engager exposure, i.e. roughly **90% functional loss**. This model reaches at
most **12%**. Exhaustion therefore never becomes the binding constraint. A treatment-free interval
costs 25% of drug-time and buys back at most a few percent of function, so it can never win, in any
geometry. The experiment was structurally incapable of detecting the hypothesised interaction.

**The most likely mechanical cause** is T-cell influx diluting the exhaustion pool. The T-cell
population grows from 200 to ~860 over the run, so most effectors present at day 28 are recent
arrivals that have accrued little contact time. Mean exhaustion is therefore held low regardless of
architecture. Whether this dilution is biologically right or a modelling artefact is itself an open
question and should be settled before re-running.

## The one signal that did survive

There is a statistically significant architecture x schedule interaction at the longest interval:

| TFI | dispersed effect | follicle effect | interaction | Mann-Whitney p | Cliff's d |
|---|---|---|---|---|---|
| 2 | 0.0013 | 0.0012 | 0.0001 | 0.97 | -0.02 |
| 4 | 0.0028 | 0.0016 | 0.0012 | 0.79 | 0.08 |
| 7 | 0.0126 | 0.0090 | 0.0036 | 0.43 | 0.22 |
| 14 | 0.0737 | 0.0396 | 0.0341 | **0.021** | 0.62 |

Read honestly: the **penalty** for interrupting treatment is about half as large in follicular
architecture as in dispersed disease. That is a real interaction in the model, but it is a statement
about sensitivity to lost drug exposure, not about exhaustion recovery, and it is one comparison out
of four with no multiplicity correction. It does not support the pre-registered chain.

## Independent threat, already established before these results

The Phase 5 trafficking test (pre-registered in `TRAFFICKING_PREREG.md` before any alternative was
run) shows that allowing T cells to squeeze past malignant B cells — which the 11 um/min lymph-node
cortex measurement arguably requires, since that speed was measured in dense tissue — raises the
follicle engaged fraction from 0.048 to 0.576, essentially the dispersed value, with T cells reaching
93% of the way to the follicle centre **without any killing at all**. If that trafficking rule is the
biologically correct one, links 1 and 2 of the chain collapse too, and P1/P2 holding is an artefact
of absolute volume exclusion.

## Current classification

Between **A (hypothesis falsified)** and **B (model artefact)**, but not yet assignable to either,
because the run that would discriminate them was not a valid test. The honest current label is
**experiment invalid, hypothesis untested**.

## What has to happen before any claim

1. Fix the exhaustion calibration so the model reproduces Philipp's 28-day curve, and verify the
   dispersed arm then reproduces the Obertopp ranking. **That is a positive control and it must pass
   before the architecture comparison means anything.**
2. Only then re-run the architecture comparison.
3. Run the full attack battery (`exp_attacks.py`): exhaustion knockout, dose-matched controls,
   trafficking knockout.
4. Do not tune parameters toward a positive result. The recalibration target is Philipp's measured
   curve, which is fixed and external, not the schedule ranking.
