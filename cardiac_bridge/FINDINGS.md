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

---

# Phase 3 — realistic response heterogeneity

**Classification: C — conditional bridge regime survives, above a clearly quantified
boundary that the cardiac-specific evidence does not reach.**

Reported per pre-registered rule D3. Stated plainly: the qualifying region occupies
**0.3% of the evidence-consistent parameter space** and is operationally indistinguishable
from failure. C is the honest label under the rules fixed in advance; it is not support.
Not a clinical recommendation.

## Controls (D)

| | Control | Result |
|---|---|---|
| PC6 | p_R = 0 reproduces the no-eradication formula | 0.0e+00 bitwise, 25 T x 300,000 |
| PC7 | perfect clearance drives carried disease to 0 | 0.00e+00 at T = 3, 6, 12 |
| PC8 | zero device hazard monotone in T | min step +1.4e-03 |
| PC9 | very high device hazard collapses T* | T* = 0 at 0.60/mo |
| PC10/11 | disabling each branch reproduces the other | 0.0e+00 bitwise |
| PC12 | **closed-form hand calculation** | max diff 5.7e-04 |
| PC13 | **all 14 parameters reachable** | pass after correction |

## The bug class, third instance, and the structural fix

Phase 2's eradication rate was inert because it was baked into the latent draw. Phase 3
repeated it with `dev_rate` and `med_unmask`: their exactly-zero sensitivity was initially
read as a "flat dead zone" finding before it turned out to be a dead code path.

Fixed structurally: `draw_latents3` now stores only unit/uniform draws and every parameter
applies at evaluation time. PC13 asserts all 14 parameters are reachable, so this cannot
recur silently.

PC13's first version demanded each parameter move the PRIMARY endpoint and flagged
`s_detect`. That was a false positive, and chasing it produced a real result: `s_detect`
moves carried disease 0.29 -> 0.45 and transplant rate 0.53 -> 0.67 while leaving five-year
survival flat. **The observation process reshuffles who is transplanted without changing
who survives** — detection-before-transplant and carriage-through-transplant are both
uniformly fatal inside the horizon. PC13 corrected to require reachability across all
outputs.

## C — identifiability

**p_R and p_dur are not separately identifiable.** 99.1% of the variance in benefit is
explained by their product alone (within-bin SD 0.15 pp vs total SD 1.82 pp). The model has
ONE effective parameter, the durable clearance fraction

    pi = p_R x p_dur

Consequence for any future trial: **response rate without matched durability data cannot
constrain this strategy at all.** An ORR readout is uninformative here by construction.

Sensitivity at the generous corner, per 10% of admissible range:

| parameter | effect |
|---|---|
| dev_rate | **-0.470 pp** |
| m_clear | -0.424 pp |
| p_R | +0.369 pp |
| p_dur | +0.275 pp |
| a_clear | +0.134 pp |
| med_unmask | +0.054 pp |
| tx_stasis | **+0.019 pp** |

Device hazard dominates; cytostasis is the weakest of the seven.

## F — the phase boundary (the deliverable)

Minimum pi reaching the 5 pp MCID. Admissible ceiling is pi = 0.30 x 0.70 = **0.21**.

| device hazard | stasis 1.0 | stasis 2.0 | stasis 3.0 |
|---|---|---|---|
| 0.010 /mo (best real) | 0.26 | 0.22 | 0.22 |
| 0.025 /mo (typical) | 0.40 | 0.36 | 0.34 |
| 0.045 /mo (poor) | 0.56 | 0.52 | 0.50 |

The clearance-versus-unmasking race, at the most favourable device and stasis:

| clearance | unmasking | required pi |
|---|---|---|
| 2 mo | 12 mo | 0.16 |
| 3 mo | 9 mo | 0.18 |
| 4 mo | 6 mo | 0.22 |
| **6 mo** | **3 mo (IGR-measured for cardiac AS)** | **0.30** |

The bottom row is the evidence-anchored one: IGR measured median time to progression of
**3 months** in cardiac angiosarcoma. At that unmasking rate the requirement is pi >= 0.30
against an available 0.00-0.07 — a factor of four or more.

## The device axis — the cleanest kill

At the cardiac-evidence pi = 0.07, with every other parameter at its most favourable
admissible value:

| device hazard /mo | benefit |
|---|---|
| 0.045 | 0.00 pp |
| 0.025 | 0.11 pp |
| 0.010 | 1.35 pp |
| 0.002 | 4.81 pp |
| **0.001** | **5.53 pp (meets MCID)** |

Clearing the bar at the evidence-level pi needs a device hazard of ~0.001/month, about
**1.2% fatal or transplant-precluding events per year**. Real TAH programmes run
0.010-0.045/month. That is a device **10-45x better than the state of the art**.

This closes the question without needing to resolve 0.22 against 0.21 — two numbers
separated by less than the resolution of the pi grid.

## Admissible-region Monte Carlo — the decision rule

| | |
|---|---|
| usable draws | 300 |
| **reaching the 5 pp MCID** | **0.3%** (threshold 50%) -> **FAILS** |
| any benefit at all (T* > 0) | 21.7% |
| median benefit | **0.00 pp** |
| 90th percentile benefit | 1.05 pp |
| max anywhere in the region | 5.31 pp (at pi = 0.183) |
| admissible pi | median 0.058, max 0.204 |

The median draw from the evidence-consistent region produces exactly zero benefit at any
bridge duration. The median admissible pi (0.058) sits almost exactly on the independent
cardiac-specific estimate (0.07).

## G — attacks

**G7. A pre-registered prediction of mine, falsified.** I predicted the Phase 2 synergy was
an artefact of cytostasis delaying detection. It is not:

| | benefit | carried disease |
|---|---|---|
| neither | 1.70 pp | 0.331 |
| BIOLOGY only (s=3, detect=1) | 2.23 pp | 0.286 |
| OBSERVATION only (s=1, detect=3) | 1.70 pp | 0.450 |
| both | 2.23 pp | 0.411 |

Biology **+0.54 pp**, observation **-0.00 pp**. The synergy is genuine biology — cytostasis
gives clearance more time to act. The observation channel moves carried disease by 12
points and survival by nothing.

**Stacked attacks at the generous corner:**

| attack | benefit | |
|---|---|---|
| no attack | 6.19 pp | meets |
| G3 non-exponential clearance (a=2.5) | 7.48 pp | *improves* |
| G6 IGR-anchored unmasking (3 mo) | 5.22 pp | meets |
| **no cytostasis at all** | **5.22 pp** | **meets** |
| G4 acquired escape (p_dur 0.70 -> 0.35) | 2.56 pp | below |
| G5 typical centre (dev 0.025) | 3.22 pp | below |
| G5 poor centre (dev 0.045) | 1.29 pp | below |
| **all together** | **0.85 pp** | below |

**Cytostasis is dispensable** — removing it entirely still clears MCID, contributing only
0.97 pp of 6.19. Phase 2 presented combination therapy as the mechanism making the bridge
viable; under responder heterogeneity it is a rounding correction. The Phase 2 synergy is
real, biological, and nearly irrelevant.

The attacks that bite are **durability and device quality**, not drug potency. Either alone
drops the corner below MCID, and neither is a pessimistic assumption.

## H — required versus observed

| pi | benefit |
|---|---|
| 0.07 (cardiac evidence) | 1.43 pp |
| 0.15 | 4.14 pp |
| **0.21 (admissible ceiling)** | **6.21 pp, meets MCID** |

- **Required pi**, single most favourable admissible setting: **>= 0.21**
- **Admissible ceiling**: 0.21 — a bound set deliberately generously in the pre-registration,
  from non-cardiac data, to give the hypothesis its best case
- **Cardiac-specific evidence**: DART cardiac 0/1, visceral 1/7; 1/7 x 2/4 durable = **0.07**;
  adjuvant chemotherapy in resected cardiac sarcoma **~0**

MCID sensitivity, reported because the verdict depends on a threshold I chose (secondary,
NOT the pre-registered rule):

| MCID | required pi |
|---|---|
| 2 pp | 0.10 |
| 3 pp | 0.15 |
| **5 pp (pre-registered)** | **0.21** |
| 8 pp | 0.30 |

Even at a 2 pp MCID the requirement (0.10) exceeds the cardiac evidence ceiling (0.07), so
the direction survives the threshold choice even where the specific verdict would not.

## I — POT1 / KDR branch, closed on evidence

No ATR-inhibitor efficacy was modelled. There is no direct preclinical evidence that POT1
deficiency confers ATR-inhibitor sensitivity in angiosarcoma; the nearest adjacent claim
(ALT-mediated ATRi hypersensitivity, Flynn *Science* 2015) was directly contradicted by
Deeg et al., who attributed the signal to cell-line background. KDR alteration in 9/11 is a
real observation from a tiny cohort and is not universal.

**Minimal discriminating experiment**, in place of simulating an invented drug effect:
cardiac angiosarcoma PDX or patient-derived lines stratified by POT1 status, dosed with a
clinical ATR inhibitor, with ATR-pathway engagement confirmed pharmacodynamically and an
isogenic POT1-restored control. Endpoint: differential surviving fraction. Without that,
any ATR arm here would be manufactured.

## The model is biased in FAVOUR of the bridge and still says no

Two structural assumptions maximise the filter's value: a patient detected during the
bridge receives no transplant and dies, and a patient carrying occult disease through
transplant is essentially always killed by it. Relaxing either makes the bridge look worse.
The negative result is obtained despite the bias, not because of it.

## What survives Phase 3

Phase 1's organ-stewardship result **does not depend on pi at all**. It requires only that
waiting reveals occult disease before an organ is committed — confirmed here, since the
observation channel moves carried disease from 0.29 to 0.45 while changing no one's
survival. The two rationales separate cleanly and Phase 3 kills exactly one:

- **Patient survival** — needs durable clearance 3-4x the cardiac evidence, or a device an
  order of magnitude beyond the state of the art. **Dead.**
- **Organ stewardship** — ~15 donor hearts saved per five-year survivor lost; futile
  transplants 65% -> 28% at a six-month bridge. **Untouched.**

The bridge is not a treatment strategy. It may still be an allocation strategy, and those
are decided by different people against different endpoints.

## Translational deliverable — the minimum prospective dataset

Because pi is the single effective parameter and is not decomposable, the registry must
measure **durable clearance**, not response rate.

**Arm 1 — the survival question.** Proportion of patients undergoing R0 resection or total
cardiectomy for primary cardiac sarcoma who are free of radiographically detectable disease
at 24 months on systemic therapy.
- **n = 45** gives 82% power to distinguish pi = 0.07 from pi = 0.20 (single-arm exact
  binomial, one-sided alpha 0.05; reject H0 if >= 7 of 45 achieve durable clearance).
- n = 30 gives 57%, n = 60 gives 87%.

**Arm 2 — the question that survived, and it is cheaper.** Futile-transplant rate against
bridge duration. Requires only observation of existing practice, no new therapy, no
randomisation.

**Mandatory covariates**, in sensitivity order, since these drive the boundary more than
the drug does:
1. centre-specific device event rate (fatal or transplant-precluding), per month
2. time from resection to radiographic recurrence (unmasking kinetics)
3. time from therapy start to best response (clearance kinetics)
4. duration of response, not just its occurrence

At roughly 300-500 primary cardiac sarcomas per year worldwide, n = 45 is accruable by an
international consortium in a few years. Nothing smaller can separate the hypothesis from
its null.

## Superseded interpretations, preserved

1. Phase 2's "the bar is an 11% chance of clearing occult disease" — an artefact of a model
   in which every patient was partially curable. Under a responder model the requirement is
   pi >= 0.21 durable clearance.
2. Phase 2's "combination therapy is the indicated regimen" — the synergy is real and
   biological but contributes ~1 pp of ~6; removing cytostasis entirely still clears MCID.
3. My Phase 3 prediction that the synergy was a decision-process artefact — falsified.
4. My reading of exactly-zero dev_rate/med_unmask sensitivity as a "flat dead zone" — that
   instance was a dead code path. The flatness at the evidence centre is real, but was only
   verified after the fix.
