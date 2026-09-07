# Findings — optimal TAH bridge duration after total cardiectomy for cardiac sarcoma

**Classification: B — interior optimum, parameter-sensitive.**
Recorded deviation from PREREG below. Nothing here is a clinical recommendation.

## Positive controls (run before any headline number)

| | Control | Result | Verdict |
|---|---|---|---|
| PC1 | median OS = 9 mo when transplanting with no bridge | 8.98 mo at calibrated k = 2.19 | **PASS** |
| PC2 | p_occ = 0 forces T* = 0 | T* = 0 | **PASS** |
| PC3 | zero device hazard forces T* to horizon | T* = 24 (horizon) | **PASS** |

k = 2.19 was the only tuned parameter, solved against the external 9-month anchor, never
against the shape of the T curve. It says immunosuppression roughly doubles the pace of
occult disease — a plausible value that was not assumed, it was recovered.

## E1/E2 — the two pre-registered endpoints disagree, and this is the result

| Bridge (mo) | RMST-60 (mo) | P(alive at 60 mo) | Transplanted | Futile transplants |
|---|---|---|---|---|
| 0 | 21.41 | **0.242** | 100% | 65% |
| 6 | **22.08** | 0.215 | 58% | 28% |
| 12 | 21.81 | 0.190 | 38% | 12% |
| 24 | 20.71 | 0.149 | 22% | 2% |

Mean survival peaks at a ~6-month bridge. Probability of five-year survival — the cure
proxy — **falls monotonically with every month waited, in 100% of 208 parameter draws.**

Mechanism, confirmed analytically against the simulation: waiting kills the curable.
A patient with no occult disease was cured by the cardiectomy itself and gains nothing
from the filter, while paying the device hazard every month. The filter works as designed
(futility 65% -> 28%) and selects against exactly the patients it was meant to protect.

## E3/E4 — robustness

- Interior optimum in **51.4% of 208 draws** consistent with the 9-month anchor.
  95% CI [44.4%, 58.4%]. **The CI straddles 50%, so R1 is not meaningfully cleared.**
- T* median 5 mo, IQR 0–15, full range 0–24. The optimum is wherever the parameters put it.
- Best-worst RMST gap: median 3.17 mo; only 5% of draws fall below the R2 negligibility
  threshold. The differences are real in size, just not stable in location.
- 42 of 250 draws discarded because the 9-month anchor was unreachable for any k in [1,15].

## E5 — what actually drives the answer

| Parameter | T* across its plausible range |
|---|---|
| device event rate 0.010 → 0.045 /mo | **24 → 10 → 2 → 0 months** |
| P(occult disease) 0.5 → 0.8 | 2 → 7 → 12 months |
| median unmasking time 3 → 12 mo | 2 → 7 → 12 → 17 months |
| metastatic survival 4 → 9 mo | 2 → 7 → 11 months |

Device hazard dominates everything. The optimal bridge duration is close to a direct
readout of how good a given centre's TAH programme is, not a property of the tumour.

## R3 — not a memorylessness artefact

Weibull unmasking with shape 1.0 / 1.5 / 2.0 / 3.0 gives T* = 7 / 7 / 8 / 8 months and
never inverts the direction. Passed.

## The finding that was not pre-registered

The model is patient-centric; donor hearts are a scarce shared resource, and the base case
makes the trade explicit. Per 100 patients reaching cardiectomy + TAH, a 6-month bridge
versus immediate transplant:

- **41.8 donor organs not consumed**
- **36.9 futile transplants averted** (futility among organs used falls 65% → 48%)
- **2.8 five-year survivors lost**
- **≈15 donor organs saved per five-year survivor lost**

Read honestly: the bridge is a poor individual therapy and a strong organ-stewardship
intervention. That reframing was not anticipated in PREREG and is reported as a
hypothesis-generating observation, not a tested endpoint.

## Recorded deviation

PREREG rule R1 classifies a result as A when an interior optimum survives in >= 50% of
draws. 51.4% meets that by the letter. It is not meaningfully cleared — the binomial CI
runs from 44.4% to 58.4%. Reporting A would overstate the evidence, so the classification
is **B**. The letter of the rule and the honest reading diverge; the honest reading wins,
and the divergence is recorded here rather than absorbed silently.

## What would have to be true before any of this touched a patient

1. Device hazard is the dominant parameter and is centre-specific. Any real use requires
   that centre's own TAH event rate, not a literature average.
2. p_occ, unmasking kinetics and k have no direct measurement in this disease. They were
   swept, and the anchor constrains their joint values, but they are not observed.
3. Systemic therapy during the bridge is not modelled. KDR-directed or DDR-directed
   treatment during the bridge would change the answer and is the obvious next extension.
4. Quality of life on TAH support is not modelled and is substantial.
5. The 9-month calibration anchor is from a small retrospective transplant series.
6. Device hazard is modelled as constant; the real hazard is front-loaded, which would
   penalise short bridges less and long bridges more than modelled here.

## Reproduce

```
python -m venv .venv && .venv/bin/pip install numpy scipy
.venv/bin/python run_controls.py     # PC1-PC3, must pass first
.venv/bin/python run_main.py         # E1-E3 base case
.venv/bin/python run_robustness.py   # E4, E5, R3
.venv/bin/python run_verdict.py      # decision rules + stewardship trade-off
```
Seeds fixed in each script. Raw outputs preserved in `results/`.
