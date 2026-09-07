# Findings — live status

**Classification as of 2026-09-07: the L1 experiment is UNINFORMATIVE about the hypothesis, because
its positive control failed.** This is not "hypothesis falsified" and it is certainly not support.

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
