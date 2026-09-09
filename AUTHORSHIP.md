# Authorship

All code, experiments, analyses, corrections and retractions on this branch were produced by
**Claude Opus 5**.

- Commits from this point carry the git author **Opus <noreply@anthropic.com>**.
- Every commit, including earlier ones authored as `Claude`, carries the trailer
  `Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>` and a session link.

Earlier commits on this branch were authored as `Claude <noreply@anthropic.com>`. They are the
same author and are **not** being rewritten: the branch is shared — Codex pushed a review commit to
it — and rewriting pushed history would invalidate other checkouts and risk losing that work. The
identity change is applied going forward and recorded here instead.

## Files authored on this branch

**Model and instrument**
`lymphoid.py` (two behaviour-preserving edits: `kill_efficiency` hook, occupancy-graded
exhaustion/recovery), `exhaustion.py`, `philipp_assay.py`, `pk.py`, `schedules.py`,
`calibration.py`, `provenance.py`

**Experiments and fitting**
`exp_replicate.py`, `exp_tarlatamab.py`, `fit_calibration.py`, `calibrate_exhaustion.py`,
`calibrate_threshold.py`, `run_calib_scan.py`, `structural_audit.py`

**Record**
`docs/recon/LUNG_CANCER_RECON.md`, `docs/plan/PLAN.md`, `docs/plan/GROUNDWORK.md`, `docs/calibration/PHILIPP_ASSAY.md`,
`docs/calibration/PHILIPP_POPULATION_DATA.md` + `.csv`, `docs/calibration/M1_RESELECTION_PREREG.md`,
`docs/findings/HELD_OUT_CONTAMINATION.md`, `docs/findings/STRUCTURAL_COUNTEREXAMPLE.md`, `docs/findings/COUNTEREXAMPLE_VERDICT.md`,
`docs/findings/STRATEGIC_ASSESSMENT.md`, `docs/findings/FINAL_CLASSIFICATION.md`, `docs/plan/TASKS_ASTRA.md`, `AUTHORSHIP.md`

**Superseded, preserved deliberately**
`superseded/calib_probe_bug/`, `superseded/calib2_density_bug/`

Not authored here: `docs/review/CODEX_GROUNDWORK_REVIEW.md`, `docs/review/CODEX_GROUNDWORK_REVIEW_CHECKS.json`,
`docs/review/CODEX_verify_groundwork.py` (Codex), and the pre-existing follicular lymphoma files.
