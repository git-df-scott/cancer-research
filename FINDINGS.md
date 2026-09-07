# Findings — final status

**CLASSIFICATION A — MODEL NOT VALIDATED FOR THE SCHEDULING QUESTION.**

Two positive controls have failed under two pre-registrations. Per `PHASE3_PREREG.md` §4 this line
of attack ends here and is written up as a negative. The model is **not** modified a third time.

**The architecture hypothesis remains UNTESTED. It is not falsified, and it is not supported.**
The architecture comparison, trafficking attack and causal tests were never run.

Experiments L1, Phase 2 and Phase 3 are all preserved, with every seed and raw output.

---

## The result that this project actually establishes

> **In this model, tumour control and binding T-cell exhaustion never co-occur.** Across a 14-fold
> range of effector-to-target ratio, in 24 of 24 (N_T × parameter) cells, there is no regime in
> which the drug controls the disease *and* exhaustion is the limiting constraint. The window in
> which an exhaustion-driven scheduling effect could exist is **empty.**

| N_T | E:T | tumour burden at d42 / n0 | terminal T-cell function | controlled? | exhausted? |
|---|---|---|---|---|---|
| 200 | 1:28 | 2.04 | 0.004 | no | **yes** |
| 400 | 1:14 | 2.03 | 0.008 | no | **yes** |
| 800 | 1:7 | 1.37 | 0.258 | no | **yes** |
| 1400 | 1:4 | 0.01 | 0.947 | **yes** | no |
| 2000 | 1:3 | 0.00 | 0.922 | **yes** | no |
| 2800 | 1:2 | 0.00 | 0.943 | **yes** | no |

This is one finding, not two accidents. It explains both failures at once:

- **Phase 2** (E:T 1:28) — exhaustion binds hard, terminal function 0.001–0.36, but the tumour is
  never controlled. It outgrows its starting burden and saturates the lattice at 78–79% occupancy.
  Arms differ by 0.6%. **The endpoint could not measure anything.**
- **Phase 3** (E:T 1:4) — the endpoint has full dynamic range (continuous leaves 22–1296 cells, the
  worst arm 7,546, a 6–129× spread, nowhere near the ceiling). The tumour is controlled, but killed
  so fast that **antigen disappears before exhaustion accrues**. Terminal function is 0.53–1.00.
  A rest period has nothing to restore and only costs drug time.

In Phase 3 every non-continuous schedule was worse than continuous, monotonically in duty cycle, in
all four combinations. That is not a marginal miss; there was nothing there to find.

## What the model is missing, stated precisely

Real patients occupy the regime this model cannot reach. Philipp's Figure 1B measured it directly:
patients on continuous blinatumomab with **persistent disease and** peripheral T-cell function down
to 0.238 of baseline. Partial control with persistent antigen and a progressively exhausting
effector pool is the clinically relevant state, and it is exactly the state that is unreachable
here.

It is unreachable for one identifiable reason. Sustaining an engaged pool against persistent
antigen requires ongoing recruitment — but recruited T cells arrive with **zero** exhaustion, so
any influx high enough to sustain the pool also resets the population mean and prevents functional
collapse. That is what Phase C measured and excluded.

Every thread converges on one sentence, recorded in `PHASE2_PREREG.md` §2 **before any of these
runs**:

> *"within this model's structure, recruitment and population-level functional collapse cannot both
> be represented at a realistic influx rate, because arrivals enter with zero exhaustion"*

## What a future attempt would need

Recruited T cells that are **not naive** — an exhaustion state that is partly systemic rather than
purely per-cell, so a sustained pool can still lose function. That is a different model, needing
its own external calibration and its own pre-registration. It is **not** attempted here;
`PHASE3_PREREG.md` §4 forbids it.

## What was established along the way, and is worth keeping

1. **An externally calibrated exhaustion submodel.** Reversible + durable, fitted to seven measured
   points from Philipp et al. by grid search over 244,800 parameter sets, scored against no
   schedule ranking of any kind. `k = 0.159/day` is identified at every point of both profiles;
   `rho` and `tau_r` are individually unidentifiable and the identified combination — the fraction
   of exhaustion a 7-day break reverses — is 4.3–55%. Verified in the real code path.
2. **A caught calibration trap.** Fitting to % specific lysis alone reproduces all five points
   perfectly while driving the assay-saturation constant to 10.0 and declaring every arm
   functionally dead. Only a second, non-saturating readout (granzyme B) exposes it.
3. **L1's influx rate is excluded by external patient data**, confirming L1's own suspicion on an
   independent measurement rather than by judgement.
4. **At least 43% of Philipp's measured TFI benefit is reduced cumulative exposure, not
   reinvigoration** (24–57% reinvigoration across the admissible class). The clinical rationale for
   treatment-free intervals is reinvigoration specifically, so this is worth stating — though the
   split is not identifiable and this is a diagnostic, not a finding.
5. **Two named failure modes for this class of experiment**: the floor (everything clears) and the
   ceiling (everything saturates). L1 hit the first, Phase 2 the second.

## Recorded but not acted on

`p_div = 1/2880` per minute is a generic "standard human tumour cell-cycle time". It is defensible
for the leukaemia-like arm the positive controls use, but **too fast for follicular lymphoma**,
which is indolent. Any future attempt must face it — and it confounds a design that deliberately
holds all biology fixed while varying only geometry.

Raw output: `results/expE3_poscontrol.json`, `results/expE3_diagnosis.txt`,
`results/expNT_sweep.json`, `results/expE_poscontrol.json`, `results/expE_diagnosis.txt`.

---

## Phase 3: the positive control at E:T 1:4

## Phase 2: the positive control failed, and why

The pre-registered criterion (PC-1) required at least one non-continuous schedule to reduce day-42
burden versus continuous dosing at p < 0.05 with a median reduction of at least 10%. Across all
four admissible (exhaustion representative × influx) combinations, ten schedules and eight paired
seeds — 320 runs — **no combination passed.**

### The endpoint saturated. The experiment had no dynamic range.

| combination | median nB42 | % of lattice | spread across all 10 arms |
|---|---|---|---|
| `rec_high` × 0 | 11,285 | 78.4% | **0.9%** |
| `rec_low` × 0 | 11,343 | 78.8% | **0.6%** |
| `rec_low` × 2.5e-5 | 9,838 | 68.3% | 11.9% |
| `rec_mid` × 0 | 11,320 | 78.6% | **0.7%** |

Initial burden 5,542; lattice capacity 14,400. Continuous-arm trajectories flatten by day 21–28
(`d21=11000, d28=11298, d35=11302, d42=11313`) **because the domain is full, not because the
disease is controlled.** A 10% minimum effect is unreachable when arms differ by under 1%.

### This is the mirror image of L1's failure

| | failure mode | why burden could not separate the arms |
|---|---|---|
| **L1** | **floor** | every arm cleared the tumour (medians 4 and 7 cells of 5,525) |
| **Phase 2** | **ceiling** | every arm saturates the lattice at 68–79% occupancy |

L1 anticipated the floor, and `PHASE2_PREREG.md` §5 carried a contingency for it. Neither
anticipated the ceiling. That is the gap, and it is now on the record.

### The causal chain, quantified

```
Phase C's external criterion (patient T-cell function must collapse at the population level)
  -> only near-zero influx is admissible
  -> with a 14-day background T-cell half-life and no replenishment, the pool falls 200 -> 10-142
  -> 10-142 T cells cannot control 5,542 malignant cells with a 2-day doubling time
  -> the tumour reaches lattice carrying capacity by day 21-28
  -> the endpoint stops responding to anything
```

**This was recorded in advance.** `PHASE2_PREREG.md` §2: *"within this model's structure,
recruitment and population-level functional collapse cannot both be represented at a realistic
influx rate, because arrivals enter with zero exhaustion."* That is exactly what bit.

### What is NOT the cause

- **Not the kill rate.** 1.9–3.0 targets/T-cell/day, inside the measured range of 2–16
  (Halle et al., Immunity 2016). Not miscalibrated.
- **Not the exhaustion mechanism.** It works exactly as calibrated.

### The exhaustion mechanism works. It just cannot pay for itself here.

In `rec_low` × 2.5e-5 — the only combination with real dynamic range, and the one flagged in
advance as physiologically sensible — end-of-run T-cell function rises **monotonically** with
off-fraction, and so does tumour burden:

| schedule | duty | T-cell function at day 42 | median nB42 |
|---|---|---|---|
| `B_7on7off` | 0.50 | **0.646** | 10,807 |
| `B_4on3off` | 0.57 | 0.545 | 10,408 |
| `B_MO_FR` | 0.71 | 0.459 | 10,174 |
| `B_6on1off` | 0.86 | 0.395 | 9,929 |
| `A_cont` | 1.00 | **0.315** | 9,678 |

Breaks restore function precisely as the external calibration says they should — 0.315 → 0.646 —
and the tumour is **larger** every time. Rank correlation between restored function and worse
control is perfect across the repeating family.

Read honestly: **recovery is real and does not repay the lost drug time.** But this is *not* a
valid test of the hypothesis, because the regime is one of uncontrolled growth against a
saturating boundary. A model that cannot control disease under continuous dosing cannot be used to
rank schedules.

### Consequence

Classification **A**. The architecture comparison, trafficking attack and causal tests were not
run. Per `PHASE2_PREREG.md` §3, any continuation requires a **new externally justified calibration
and a new pre-registration**. The pre-registered non-naive-arrivals variant (§8) is explicitly
**not permitted as a rescue** here, and has not been run.

Raw output: `results/expE_poscontrol.json` (320 runs, all seeds preserved),
`results/expE_poscontrol_verdict.json`, `results/expE_diagnosis.txt`.

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
