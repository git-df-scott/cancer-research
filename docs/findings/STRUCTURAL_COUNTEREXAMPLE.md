# Phase 5: a structural counterexample to M1's transferability

M1 reproduces Philipp's three continuous-arm numbers. This document records a structural failure
that is independent of fit quality, found by attacking limiting behaviour rather than by fitting
anything.

## The counterexample

**M1's chronic culture has no T-cell proliferation. Philipp's does, and it differs about four-fold
between the two arms M1 is fitted to.**

Model, calibration configuration (L=60, E:T 1:4, `t_div = 0`):

| day | T cells | % of start |
|---|---|---|
| 0 | 360 | 100 |
| 7 | 237 | 66 |
| 14 | 159 | 44 |
| **28** | **64** | **18** |
| 56 | 10 | 3 |
| 70 | 0 | **extinct** |

Philipp, same window (Fig 3, CD2+ fold change after 3 days of assay):

| | continuous | TFI |
|---|---|---|
| day 14 | 1.1 | **4.1** |
| day 28 | 0.06 | **2.8** |

The real culture sustains and expands its T cells, and expands them roughly four times more in the
rested arm. This model's culture only loses them, monotonically, to extinction.

## Why this is a counterexample rather than a nuisance

M1 attributes **all** of the TFI arm's recovered cytotoxicity to per-cell exhaustion decay through
`recover_tau`. The source experiment shows the rested arm also has a large proliferative response
the continuous arm lacks.

So the observable M1 is fitted to — 93.4% specific lysis at day 14 after a rest — has at least two
contributions in reality:

1. per-cell reinvigoration of exhausted T cells,
2. arm-dependent expansion changing the population's composition, plausibly diluting exhausted
   cells with fresh daughters.

M1 can only represent the first, so it necessarily assigns the whole effect to it. **It fits the
right number with a mechanism decomposition the source experiment contradicts.** A model that gets
the right answer by attributing an effect to the wrong cause does not transfer: in R and T the two
contributions would scale differently with schedule, geometry and exposure.

This is exactly the failure mode the phase was designed to find — excellent calibration fit,
structurally wrong.

## Secondary findings from the same audit

**A test passed for the wrong reason.** "Very long stimulation stays bounded in [0,1]" originally
passed with E(56d) = 0.000. That was not boundedness; the population had gone extinct and the mean
of an empty array was being reported as 0.0. The helper now returns NaN on extinction and the test
fails honestly. A test that passes on an empty population is worse than no test, and this one had
been reporting PASS.

**Nine limiting behaviours are sound.** Zero drug gives exactly zero drug-attributable exhaustion;
occupancy 1e-6 gives 3.1e-7; exhaustion is monotone in occupancy with half occupancy properly
intermediate (0.187 vs 0.477); drug without antigen gives no exhaustion; withdrawal with antigen
present permits recovery while sustained engagement does not; lysis stays within [0,100] and is
monotone in exhaustion; the timestep dependence is 0.9% between dt=5 and dt=1; a 1% parameter
perturbation moves the observable 1.8%.

So the graded-occupancy correction works as intended and the model is not numerically pathological.
The failure is biological structure, not arithmetic.

## Scope of the claim

This is a counterexample to the **chronic culture implementation shared by M0 and M1**, not to the
Hill mapping specifically. Fixing it would require modelling proliferation in the chronic culture
and then re-deciding how much of the TFI recovery is per-cell versus compositional — which is a
new mechanism, and this block forbids adding one.

It also means the earlier density-standardisation fix, while correct, papered over a related
problem: resampling 250 cells from 64 survivors preserves the exhaustion distribution but cannot
recover the population dynamics that were never simulated.

## What it does not show

It does not show M1's Hill mapping is wrong. It shows the calibration cannot distinguish a
per-cell recovery mechanism from a compositional one, because the model omits the compositional
channel entirely. That is an identifiability failure with a specific, named, measured cause.

---

**SUPERSEDED IN PART — see `docs/findings/COUNTEREXAMPLE_VERDICT.md`.** The central claim here, that
Philipp's chronic culture is measured to sustain and expand T cells while this model's does not,
is RETRACTED. The 4.1-vs-1.1 fold change is a 3-day assay readout of proliferative capacity, not a
chronic-culture population trajectory, and no absolute chronic-culture counts are reported. What
survives is narrower: two functional axes were measured to diverge between arms, and M1 has one
latent state.
