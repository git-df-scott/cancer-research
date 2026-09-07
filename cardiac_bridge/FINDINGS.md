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

---

# Phase 2 — systemic therapy during the bridge

**Phase 2 inverts the Phase 1 conclusion.** Phase 1's "never wait" holds only if the bridge
is empty time. It is not a property of the bridge; it is a property of an untreated bridge.

## Controls

| | Control | Result | Verdict |
|---|---|---|---|
| PC4 | therapy off reproduces the Phase 1 formula | max element-wise difference 0.0e+00 across 25 bridge durations x 600,000 patients | **PASS (bitwise)** |
| PC5 | eradication machinery actually fires | carried disease at T=12: 0.121 -> 0.066 -> 0.003 -> 0.000 as rate rises | **PASS** |

## Two control failures, found and fixed rather than reasoned past

**PC4, first version.** Compared therapy-off against a *stored* Phase 1 curve that had been
computed under a tighter `calibrate_k` tolerance, so it could never match bitwise. It was
testing pipeline reproducibility across a tolerance change, not the therapy code. Replaced
with a literal Phase 1 reference implementation compared element-wise at identical k.

**Eradication was a silent no-op.** `tx_erad_rate` was baked into the latent draw. Latents
are drawn once and the rate varied afterwards via `dataclasses.replace()`, so every
eradication rate returned byte-identical output. The pre-registered endpoint did not catch
this — it reported "no tested rate moved T* off zero", which reads exactly like a clean
negative result and was a dead code path. What caught it was checking the reported
magnitude against a hand calculation: at 84% clearance the curable fraction should rise
from 0.35 to ~0.89, and the model was returning the untreated value to four decimal places.
Fixed by drawing a unit exponential and applying the rate at evaluation time. PC5 added to
test the mechanism directly rather than trusting an endpoint.

This is worth recording plainly: a pre-registered negative result was one commit away from
being reported, and it was a bug.

## F2 — the efficacy bar, and it is low

| Eradication rate | Clearance by 6 mo | Best bridge | P(alive 5yr) | vs transplant now |
|---|---|---|---|---|
| 0.00 | 0% | 0 mo | 0.2421 | — |
| 0.02 | **11%** | 4 mo | 0.2493 | +3% |
| 0.05 | 26% | 7 mo | 0.2916 | **+20%** |
| 0.10 | 45% | 8 mo | 0.3523 | **+46%** |
| 0.30 | 84% | 6 mo | 0.4781 | **+97%** |

**The bar is an ~11% chance of clearing occult micrometastatic disease over six months.**
Per F-R1 this is REACHABLE, not aspirational: ipilimumab + nivolumab produced a 25%
response rate with durable responses in angiosarcoma (SWOG S1609 / DART).

## F3 — the bar is centre-specific

| Device hazard | Required clearance by 6 mo |
|---|---|
| 0.010 /mo | 3% |
| 0.020 /mo | 6% |
| 0.025 /mo | 11% |
| 0.035 /mo | 11% |
| 0.045 /mo | 16% |

Confirmed as predicted. A strong TAH programme needs a much weaker drug to justify waiting.

## F1 — cytostasis alone: right direction, wrong mechanism of harm

| Stasis | T*(RMST) | T*(5yr) | P(alive 5yr) @T=6 | Carried disease @T=6 |
|---|---|---|---|---|
| 1.0 | 6 | 0 | 0.2146 | 0.280 |
| 2.0 | 14 | 0 | 0.2145 | 0.396 |
| 3.0 | 24 | 0 | 0.2145 | 0.444 |

F1 predicted cytostasis would be harmful. Direction confirmed, but the survival effect is
-0.0001 — negligible. The real damage is elsewhere: cytostasis raises carried occult
disease by 59% relative and drags the RMST-optimal bridge from 6 months to the full
horizon. **A cytostatic agent does not harm the patient, it corrupts the decision rule.**
A centre choosing bridge duration by mean survival, in patients on effective growth
suppression, is driven toward ever-longer bridges that burn device-months and produce more
futile transplants.

## F-E5 — and the interaction reverses the practical reading

| | erad 0.0 | erad 0.05 | erad 0.10 |
|---|---|---|---|
| stasis 1.0 | 0.2416 (T=0) | 0.2907 (T=7) | 0.3509 (T=8) |
| stasis 2.0 | 0.2416 (T=0) | 0.3094 (T=10) | 0.3812 (T=10) |
| stasis 3.0 | 0.2416 (T=0) | 0.3193 (T=11) | 0.3953 (T=11) |

Interaction **+0.0186: SYNERGY, not antagonism.** Cytostasis contributes exactly zero on
its own and potentiates eradication substantially. Mechanism: slowing growth holds occult
disease inside the treatable-and-still-hidden window for longer, giving the eradicating
agent more time to act before the patient is either detected or transplanted.

The naive reading of F1 alone — "cytostatic therapy is the wrong drug" — is wrong, and is
recorded here because it was stated before F-E5 was run. The correct reading is that a
cytostatic agent is **insufficient alone and valuable in combination.**

## What Phase 2 says, in one paragraph

The bridge's value is entirely contingent on whether bridge-period therapy can eradicate
micrometastatic disease rather than merely control it, and the bar is low enough to be
clinically reachable (~11% clearance over six months). Cytostatic therapy contributes
nothing alone but potentiates an eradicating agent, so the indicated regimen is a
combination rather than either component. Notably, the TAH bridge is the only window in
this disease in which the patient carries no allograft and therefore no immunosuppression
opposing an immune-based therapy — and immunotherapy is precisely the eradication-type
modality the model requires. The model was not built to find that; it fell out of
separating two mechanisms the pre-registration forced apart.

## Additional Phase 2 limits

- Therapy is modelled as active only during the bridge and stopping at transplant (F-R3).
- Eradication is modelled as a constant hazard. Real immunotherapy response is not
  exponential in time and is concentrated in a responder subpopulation.
- The 25% angiosarcoma response rate is from a 16-patient cohort dominated by non-cardiac
  primaries; scalp and face tumours (UV-driven, mutation-rich) drove most of it. Cardiac
  angiosarcoma is visceral and likely less immunogenic, so mapping that response rate onto
  this eradication hazard is an assumption, not a measurement.
