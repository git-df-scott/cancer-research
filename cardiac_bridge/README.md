# TAH bridge duration after total cardiectomy for primary cardiac sarcoma

A competing-risks decision model asking one narrow question, pre-registered before running.

**Status after Phase 3: classification C — a conditional bridge regime exists, above a
quantified boundary that the cardiac-specific evidence does not reach. The qualifying
region occupies 0.3% of the evidence-consistent parameter space and is operationally
indistinguishable from failure. Not a clinical recommendation.**

Read `PREREG.md` first, then `PROVENANCE.md` (evidence audit), then `FINDINGS.md`.

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

**Phase 1 (untreated bridge).** Mean survival prefers a ~6-month bridge; five-year survival
prefers transplanting immediately, unanimously across 208 draws. The bridge trades roughly
15 donor organs saved per five-year survivor lost.

**Phase 2 (constant eradication hazard).** Appeared to invert Phase 1: an ~11% clearance
chance made waiting beneficial. Superseded — that model made every patient partially
curable.

**Phase 3 (latent responder model).** The strategy turns on ONE parameter, the durable
clearance fraction pi = p_R x p_dur (99.1% of variance; the two are not separately
identifiable). Required pi >= 0.21 at best case; cardiac-specific evidence gives
pi ~ 0.00-0.07. Only 0.3% of the evidence-consistent region reaches a 5 pp survival gain.
At the evidence-level pi the strategy needs a total artificial heart 10-45x better than the
state of the art. The patient-survival rationale is dead; the organ-stewardship rationale
never depended on pi and is untouched.

## Layout

| Path | What it is |
|---|---|
| `PREREG.md` | pre-registration: endpoints E1-E5, controls PC1-PC3, rules R1-R4, fixed in advance |
| `bridge_model.py` | the model: latent occult disease, unmasking, device hazard, transplant branches |
| `run_controls.py` | PC1-PC3; must pass before any result is reportable |
| `run_main.py` | base-case bridge-duration curve |
| `run_robustness.py` | parameter-space sampling with per-draw recalibration, Weibull check |
| `run_verdict.py` | decision rules and the organ-stewardship trade-off |
| `PROVENANCE.md` | Phase 3 evidence audit: every load-bearing number, its population, endpoint, source and applicability |
| `phase3_model.py` | latent responder model; all parameters applied at evaluation time |
| `run_phase3_controls.py` | PC6-PC13 including a closed-form check and an all-parameters-live guard |
| `run_phase3_identifiability.py` | is the model constrainable from external evidence at all |
| `run_phase3_final.py` | phase boundary and admissible-region Monte Carlo |
| `run_phase3_attacks2.py` | G-series attacks; biology vs observation separation |
| `run_phase3_registry.py` | device axis and the prospective dataset sizing |
| `FINDINGS.md` | live results, every superseded interpretation, bug and correction |
| `results/` | raw outputs, seeds preserved |
