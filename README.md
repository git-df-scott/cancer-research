# T-cell engager scheduling: follicular lymphoma, and an attempt to transfer it to lung

Two linked projects under active adversarial testing. **Nothing here is a validated biological
claim, and no experiment licensed to make treatment predictions has been run.**

## Status

| | |
|---|---|
| **Follicular lymphoma (original)** | Experiment L1 invalid. Its schedule encoding did not match the reference it was compared against, and the exhaustion calibration was independently broken. Not yet re-run. |
| **Lung / SCLC (current)** | Calibration classified **B — FIT BUT NOT VALIDATED**. Experiments R and T are **blocked by a hard gate** and have never been run. |

Start with `docs/findings/FINAL_CLASSIFICATION.md`, then `docs/findings/STRATEGIC_ASSESSMENT.md`.

## Layout

```
docs/recon/         why SCLC, prior art, the biology
docs/plan/          phase structure, groundwork state, task board, handoff
docs/calibration/   the Philipp assay as actually performed, data, pre-registrations
docs/findings/      classification, contamination, counterexamples, retractions
docs/review/        Codex's independent review
results/            raw per-run outputs, every seed preserved
superseded/         invalidated work, kept deliberately
```

Model and instrument: `lymphoid.py`, `exhaustion.py`, `philipp_assay.py`, `pk.py`, `schedules.py`,
`calibration.py`, `provenance.py`.
Experiments and fitting: `exp_*.py`, `fit_calibration.py`, `calibrate_*.py`, `structural_audit.py`.

## The gate

`exp_replicate.py` and `exp_tarlatamab.py` refuse to run unless `results/calibration.json` reads
`VALIDATED`. It currently reads `FIT_NOT_VALIDATED`. A fitted model is not a validated one; the
gate exists because documenting that distinction was not enough.

```bash
python -m venv .venv && .venv/bin/pip install numpy scipy
.venv/bin/python philipp_assay.py      # observation-model invariants
.venv/bin/python pk.py                 # PK self-tests
.venv/bin/python provenance.py         # source-to-parameter contamination audit
.venv/bin/python calibration.py        # gate refusal tests
.venv/bin/python structural_audit.py   # limiting-behaviour attack
```

## What has actually been established

**Positive.** A threshold functional mapping reproduces Philipp's chronic-stimulation cytotoxicity
curve within simulation noise, where a linear mapping cannot, on a bracketed grid. The L1 schedule
encoding is provably wrong by exposure arithmetic. Tarlatamab's trough occupancy drop under
clinical Q2W dosing is 17–22% at an EC50 of ~1 nM and 37–45% at 3 nM — a result that needed no
simulation and has survived every correction.

**Negative, and more useful.** Four instrument bugs were found and fixed, none of which raised an
error and all of which produced plausible numbers. The held-out observation turned out to be
contaminated by the parameter it was meant to test. Two measured functional axes are compressed
into one latent variable, so the mechanism of treatment-free-interval recovery is unidentified.

Full record, including retracted claims, in `docs/findings/`.

## Standing rules

Report every pre-registered endpoint regardless of outcome. Calibrate to external measured curves,
never to a schedule ranking. Preserve failed hypotheses and superseded results rather than deleting
them. No cure claims, no clinical recommendations, no proposed dosing changes. A clean negative is
a success.
