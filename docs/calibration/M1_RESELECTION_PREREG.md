# Phase 1 pre-registration: clean M1 reselection

**Written and committed before the search is run.** Nothing below may be changed after seeing
results.

## Why this is necessary

M1's current `theta = 0.5, hill = 4.0` were the best cell of the threshold scan run against the
**destructive probe** — an instrument since shown invalid. M1's shape was therefore selected using
contaminated information, and its SSE of 9.0 is inflated by a selection effect of unknown size.

The purpose here is to **kill M1**, not to rescue it. If an independently selected M1 cannot
reproduce the calibration under the predeclared criterion, M1 is dead and the classification
becomes D.

## Parameter ranges, justified independently of the old optimum

**theta** — the latent signal at which cytotoxic function is halved. `lymphoid.py` clips `E` to
[0, 1], so any threshold above 1 is unreachable and any at 0 means permanent full exhaustion. The
admissible domain is therefore (0, 1]. Sampled across it: **0.15, 0.30, 0.50, 0.70, 0.90**.

This ladder is a regular spread over the reachable domain. 0.50 appears because a five-point
spread over (0,1] contains it, not because it was the previous optimum. The grid is not
concentrated around it, and its neighbours are 0.30 and 0.70 rather than 0.45 and 0.55.

**hill** — steepness of the functional mapping. `hill = 1` is a hyperbolic, near-graded loss;
large `hill` approaches the hard PD-1 threshold the reference model uses. A doubling ladder spans
that structural range without privileging any point: **1, 2, 4, 8, 16**.

4.0 appears as a rung of a doubling ladder from 1, not as a centre. Its neighbours are 2 and 8.

**p_kill** and **exhaust_tonic** are searched jointly, because theta and hill trade against them:
p_kill **5e-4, 1e-3, 2e-3**, tonic **2.5e-5, 5e-5, 1e-4**. Ranges carried from the bracketed refit,
whose optimum was interior.

Total: 5 x 5 x 3 x 3 = **225 cells**. Fixed. No extension, no refinement, no second pass.

## Frozen protocol

| | |
|---|---|
| observations | Philipp d7 = 88.4, d14 = 34.9, d28 = 8.6 (continuous arm) |
| held out | d14_tfi, and it is CONTAMINATED, so a pass counts as nothing |
| weights | uniform, unweighted SSE in percentage points |
| seeds | model seeds 0 and 1; assay seeds 900, 901 |
| search | exhaustive grid, no adaptive sampling |
| stopping | when all 225 cells are evaluated. No early stop, no widening |

## Predeclared acceptance criterion

M1 reproduces the calibration if some cell satisfies **both**:

1. **SSE < 50**, and
2. **every individual |residual| ≤ 10 percentage points**

Justification, fixed before running: seed-to-seed standard deviations measured in the refit were
0.2–4.0 points per target. A model reproducing the data within its own stochastic noise would score
roughly 3 × 4² ≈ 50 by chance alone, so SSE < 50 is the natural noise floor. The per-residual cap
is added because a low SSE can hide one large error offset by two small ones, and a 10-point miss
on any single target is not "reproduces the observation".

Both conditions are required. Neither alone suffices.

## Outcome rule, fixed in advance

- **No cell satisfies both criteria** → M1 is killed. Classification **D** for the M0/M1 class. Stop.
- **At least one cell satisfies both** → M1 survives Phase 1; proceed to identifiability (Phase 2).

The full near-optimal region is reported either way, not only the winner.

## What this cannot establish

Passing means an independently selected M1 reproduces three numbers from one paper. It is
candidate generation, not validation. `provenance.py` shows no clean external observation exists,
so classification A remains unreachable regardless of the outcome here.
