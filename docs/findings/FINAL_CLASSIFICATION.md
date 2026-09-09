# Final classification: **B — FIT BUT NOT VALIDATED**

The model reproduces the calibration observations. No legitimate independent test establishes
predictive performance. **Experiments R and T remain blocked**, enforced by `calibration.py`.

## The evidence

Both optima are interior and bracketed on a grid spanning p_kill 3e-4 to 4e-3 and exhaust_tonic
6e-6 to 2e-4, against the density-corrected observation model.

| | best SSE | d7 | d14 | d28 | residuals | seed sd |
|---|---|---|---|---|---|---|
| M0 linear | 384.2 | 75.6 | 47.0 | 0.0 | −12.8 / +12.1 / −8.6 | 4.0 / 2.6 / 0.0 |
| **M1 threshold** | **9.0** | 87.0 | 33.4 | 10.8 | −1.4 / −1.5 / +2.2 | 0.2 / 3.8 / 1.6 |
| target | | 88.4 | 34.9 | 8.6 | | |

**M0 is rejected by the data.** Its residuals are 3–5× its own seed noise and carry an alternating
(−,+,−) sign pattern, the signature of shape mismatch rather than scatter. No M0 cell scores under
SSE 50 anywhere on the grid.

**M1 fits within simulation noise.** Residuals of −1.4, −1.5, +2.2 against seed standard deviations
of 0.2, 3.8, 1.6. Exactly one cell scores under 50, and the next best is 641.2 — a 70× gap, so on
this grid the optimum is sharply identified rather than lying on a ridge.

## Why not D

Model class rejection requires that no defensible parameterisation reproduces the observations.
M1 does, to within the stochastic uncertainty of the simulation producing it. The earlier move
toward D is withdrawn: it rested on a recovery residual that turned out to be the plating-density
artefact.

## Why not C

C requires multiple externally defensible explanations that remain compatible with the observations
and disagree downstream. Only one survives. M0 is excluded by residual structure, and M1's
admissible region on this grid is a single cell. Underdetermination is not what is wrong here.

## Why not A, and this is the binding constraint

A requires a preregistered, provenance-clean external prediction. **None exists, and none can be
constructed from the available data.** `provenance.py` asserts it permanently: `lymphoid.py:91`
sets `recover_tau = 10080` citing *"TFI reinvigoration (Philipp 2022)"*, the same paper every
calibration value comes from. All five Philipp values — including d28_tfi = 58.7%, found late and
recorded — are contaminated through both `recover_tau` and `exhaust_tonic`.

A pass on a contaminated target is uninterpretable in the favourable direction. M1's held-out
reading of 86.6 against 93.4 is therefore **not evidence of predictive validity**, and is not
counted as such.

## Exactly what is missing

1. **A provenance-clean external observation.** Same or comparable mechanism, quantitatively
   extractable, from a source that informed no model default and did not drive M0-versus-M1
   selection. A targeted search surfaced candidates — a 2025 bioRxiv chronic-antigen-stimulation
   study in solid tumours, and PMID 34258584 on tumour burden and exhaustion — but neither was
   assessed against the seven admissibility criteria, and neither is claimed as usable.
2. **M1's structural parameters, re-selected honestly.** `theta = 0.5` and `hill = 4.0` were the
   best cell of the **superseded** threshold scan, run against the destructive probe. M1's shape
   was chosen using an instrument since shown invalid and has never been re-searched against the
   corrected readout. SSE 9.0 is therefore inflated by a selection effect of unknown size. This is
   the single most likely thing to change the picture.
3. **Identification outside the responsive band.** The readout is compressive: at the fitted kill
   rate it is near-saturated across a wide range of exhaustion and collapses only near the top.
   Where dY/dE ≈ 0 the assay carries no information about latent E. Parameters are identified in a
   narrow band; R and T would operate largely outside it.
4. **Grid resolution.** One cell under SSE 50 with the next at 641.2 is consistent with sharp
   identification, and also with a grid too coarse to resolve the basin. Not distinguished.

## What this does and does not license

Licensed: reporting that a threshold functional mapping reproduces Philipp's chronic-stimulation
curve within simulation noise where a linear mapping cannot, on a bracketed grid, under an
observation model built from the published protocol.

**Not licensed:** any statement about treatment schedules, in follicular lymphoma or in lung. The
inferential chain from AMG 562 / CD19 / in vitro to tarlatamab / DLL3 / in vivo crosses molecule,
target, compartment and in-vitro-to-in-vivo boundaries that this calibration constrains not at all.

## Status of the gate

`results/calibration.json` records status `FIT_NOT_VALIDATED`. `require_validated()` refuses R and
T. Verified live, not merely asserted.

The most defensible result the project has produced remains the PK occupancy analysis, which
required no simulation and has survived every correction round untouched.
