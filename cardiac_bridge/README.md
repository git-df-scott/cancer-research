# TAH bridge duration after total cardiectomy for primary cardiac sarcoma

A competing-risks decision model asking one narrow question, pre-registered before running.

**Status: classification B — interior optimum, parameter-sensitive. Not a clinical
recommendation.** Read `PREREG.md` first, then `FINDINGS.md`.

## The question

Total cardiectomy with a total artificial heart is documented in primary cardiac sarcoma.
It gives R0 resection by definition, avoids immunosuppression (which accelerates
angiosarcoma — median survival after heart transplant is 9 months for angiosarcoma versus
36 for other cardiac sarcoma histologies), and uses the bridge period as a biological
filter that lets occult metastases unmask before an organ is committed.

The surgical literature states that rationale. It does not say **how long to wait**.
Bridge duration is currently chosen by clinical judgement. This is a first attempt to
compute it.

## What came out

Mean survival prefers a ~6-month bridge. Five-year survival prefers transplanting
immediately, unanimously across 208 parameter draws. The bridge trades roughly 15 donor
organs saved per five-year survivor lost — which makes it a weak individual therapy and a
strong organ-stewardship intervention. Optimal duration is driven mostly by device hazard,
so it is a property of the centre rather than of the tumour.

## Layout

| Path | What it is |
|---|---|
| `PREREG.md` | pre-registration: endpoints E1-E5, controls PC1-PC3, rules R1-R4, fixed in advance |
| `bridge_model.py` | the model: latent occult disease, unmasking, device hazard, transplant branches |
| `run_controls.py` | PC1-PC3; must pass before any result is reportable |
| `run_main.py` | base-case bridge-duration curve |
| `run_robustness.py` | parameter-space sampling with per-draw recalibration, Weibull check |
| `run_verdict.py` | decision rules and the organ-stewardship trade-off |
| `FINDINGS.md` | live results, recorded deviation, and limits |
| `results/` | raw outputs, seeds preserved |
