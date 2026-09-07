# Phase 3 evidence audit — parameter provenance

Every load-bearing number verified against a primary source where one exists. Confidence
is graded on applicability to **occult micrometastatic disease from a CARDIAC primary**,
which is the setting the bridge model actually describes.

## The extrapolation problem, quantified up front

Phase 2 mapped "25% ORR in angiosarcoma" onto a constant eradication hazard. Reading the
DART primary source makes the weakness of that mapping explicit:

| DART / SWOG S1609 cohort 51 | value |
|---|---|
| registered / evaluable | 18 / 16 |
| confirmed responses | 4 (**1 CR**, 3 PR) — ORR 25% |
| CR rate | **1/16 = 6.3%** |
| cutaneous scalp/face | 3 responses / 5 patients (60%) |
| radiation-associated breast | 1 response |
| non-cutaneous total | 7 (breast 4, liver 2, **heart 1**, spleen 1) |
| **cardiac primary** | **n = 1, did not respond** |
| response durations | 5, 7, 12+, 13+ months (2 ongoing) |
| 6-month PFS | 38% |
| TMB-high | 1 of 7 assessed; that patient responded |

Four separate degradations sit between that headline and this model:
1. **Site.** The signal is concentrated in UV-driven scalp/face and a radiation-associated
   breast primary — both mutation-rich contexts. Cardiac angiosarcoma is visceral and
   neither UV- nor (usually) radiation-associated. n = 1 cardiac, and it did not respond.
2. **Endpoint.** ORR is RECIST shrinkage of bulky disease. Only 1/16 reached CR.
3. **Durability.** Two of four responses had ended by 7 months.
4. **Disease state.** Measurable metastases, not occult micrometastatic disease.

## The closest real analogue to the bridge setting, and it is negative

Adjuvant doxorubicin-based chemotherapy after resection of primary cardiac sarcoma
(IGR series, PMC2063231) is the nearest published approximation to "resected cardiac
sarcoma plus immediate systemic therapy against micrometastatic disease":

| | value |
|---|---|
| patients | 15 (6 cardiac angiosarcoma) |
| relapsed | **13 / 15** |
| relapsed *during* therapy | **5** |
| median OS | 12 months; 2-year survival 26% |
| median TTP, angiosarcoma | **3 months** (vs 14 for other histologies, P<0.01) |
| authors' conclusion | "Post-operative conventional doxorubicin-based chemotherapy **failed to modify the natural history** of patients with resected cardiac sarcomas." |

This is direct, on-target, negative evidence for chemotherapy-mediated eradication in
exactly the population the bridge model describes.

## Provenance table

| Parameter | Estimate / range | Population | Endpoint | Source | Confidence | Directly applicable? |
|---|---|---|---|---|---|---|
| `p_R` immunotherapy, cardiac/visceral | 0.00–0.14 (1/7 non-cutaneous confirmed; 0/1 cardiac) | non-cutaneous angiosarcoma | RECIST ORR | DART S1609 c51 | **very low** (n=7, n=1 cardiac) | no — bulky disease, wrong sites |
| `p_R` immunotherapy, scalp/face | 0.60 (3/5) | UV-driven cutaneous | RECIST ORR | DART S1609 c51 | low (n=5) | **no** — different mutational aetiology |
| CR rate, immunotherapy | 0.063 (1/16) | mixed angiosarcoma | confirmed CR | DART S1609 c51 | low | partial — CR ≠ eradication |
| CR rate, weekly paclitaxel | 0.13 (9/68) | metastatic angiosarcoma | RECIST CR | Italiano 2012, Cancer | moderate | partial |
| CR rate, doxorubicin | 0.06 (2/34) | metastatic angiosarcoma | RECIST CR | Italiano 2012, Cancer | moderate | partial |
| chemo eradication, adjuvant cardiac | **≈ 0** | resected cardiac sarcoma | relapse-free survival | IGR, PMC2063231 | moderate (n=15, n=6 AS) | **yes — closest analogue** |
| durable fraction among responders | 0.50 (2/4 ongoing at ~12 mo) | mixed angiosarcoma | duration of response | DART S1609 c51 | **very low** (n=4) | no |
| cytostasis strength `s` | 2.0–2.9 (PFS 4.6 vs 1.6 mo) | soft tissue sarcoma | PFS ratio vs placebo | PALETTE | moderate | partial — not cardiac AS |
| cytostasis, cutaneous AS | ORR 31.8%, DCR 63.6%, PFS 2.8 mo | cutaneous angiosarcoma | PFS | JCOG pazopanib | moderate | partial |
| unmasking kinetics, cardiac AS | median TTP **3 months** | resected cardiac sarcoma | time to progression | IGR, PMC2063231 | moderate | **yes** |
| post-transplant OS, cardiac AS | 9 months (vs 36 other histologies) | transplanted cardiac sarcoma | median OS | transplant series | low (small retrospective) | yes — PC1 anchor |

## POT1 / ATR branch — insufficient evidence, stated plainly

- POT1 germline alterations in 45.5% of cardiac angiosarcoma and KDR alterations in 9/11
  are real findings, but n = 11 is a tiny cohort. KDR involvement must not be called
  universal on that basis.
- **There is no direct preclinical evidence that POT1 deficiency confers ATR-inhibitor
  sensitivity in angiosarcoma.** The nearest adjacent claim is ALT-mediated ATRi
  hypersensitivity (Flynn, Science 2015), and that was **directly contradicted** by Deeg
  et al. (PMC4993795), who could not confirm general hypersensitivity and attributed the
  original signal to cell-line background rather than ALT.
- **Therefore Phase 3 models no ATR-inhibitor efficacy at all.** Simulating an invented
  drug effect here would manufacture the result. The minimal discriminating experiment is
  specified in FINDINGS.md instead.

## Admissible region carried into Phase 3 (fixed before modelling)

Derived from the table above, deliberately generous on the upper side so the hypothesis
gets its best defensible case:

| Parameter | Admissible range | Basis |
|---|---|---|
| `p_R` responder fraction | **0.00 – 0.30** | 0/1 cardiac and 1/7 visceral in DART; ~0 for adjuvant chemo; upper bound stretched well above any observation to be generous |
| `p_dur` durable among cleared | **0.20 – 0.70** | 2/4 DART responses ongoing; enormous uncertainty |
| median clearance time | **2 – 8 months** | time-to-first-response in solid tumours plus consolidation |
| clearance Weibull shape | **1.0 – 2.5** | exponential to growth-like |
| cytostasis `s` | **1.0 – 3.0** | PALETTE PFS ratio 2.9; JCOG angiosarcoma weaker |
| median unmasking | **3 – 12 months** | IGR cardiac AS median TTP 3 months anchors the fast end |
| device hazard | **0.010 – 0.045 /mo** | Phase 1 range, centre-specific |
