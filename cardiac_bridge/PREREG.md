# Pre-registration — optimal bridge duration on a total artificial heart after
# total cardiectomy for primary cardiac sarcoma

Written before any model code was executed. Decision rules below are fixed in advance.
Deviations, if any, are recorded in FINDINGS.md and not silently absorbed.

## The clinical decision

A patient with primary cardiac sarcoma (predominantly angiosarcoma) undergoes total
cardiectomy with implantation of a total artificial heart (TAH). This is already done —
documented in case reports and small series (SynCardia TAH after sarcoma excision;
HeartMate 6 in TAH configuration after total cardiectomy; RV intimal spindle cell sarcoma).

Total cardiectomy plus TAH gives three things simultaneously:
  1. R0 resection by definition — the tumour bed is removed.
  2. No immunosuppression — there is no allograft to reject, so occult micrometastatic
     disease is not being pharmacologically accelerated.
  3. Time as a biological filter — during the bridge, occult metastases can become
     radiographically detectable, identifying patients for whom transplant would be futile.

The surgical literature states this rationale explicitly. What it does not state is
**how long to wait**. Bridge duration is currently chosen by clinical judgement.

## The question

Is there an interior optimum bridge duration T* > 0, and if so where, and over what
fraction of plausible parameter space does it survive?

## The tension being modelled

Longer bridge  -> more occult disease unmasks -> fewer futile transplants
Longer bridge  -> more device attrition (thromboembolism, infection, device failure,
                  deconditioning) -> more patients never reach transplant at all

If both arms are real, the objective is non-monotonic in T and an interior optimum exists.
If either arm dominates across plausible parameters, the answer is a corner solution
(transplant immediately, or never transplant) and that is a finding, not a failure.

## Model structure (fixed before running)

Discrete-event Monte Carlo, per patient, clock starting at cardiectomy.

Latent state at t=0:
  occult ~ Bernoulli(p_occ)     micrometastatic disease present but radiographically silent

Competing processes:
  U  ~ time to radiographic unmasking of occult disease (occult patients only)
  D  ~ time to a fatal or transplant-precluding device event
  Surveillance imaging every q months; detection occurs at the first imaging time >= U.

Decision: transplant at time T if no metastases have been detected and no device event
has occurred.

Outcome branches:
  A. Device event first            -> death on device at D
  B. Metastases detected first     -> no transplant; remains on device with metastatic
                                      disease; survival governed by whichever of
                                      metastatic progression or device hazard comes first
  C. Transplant at T, truly free   -> perioperative mortality, then long-term graft survival
  D. Transplant at T, occult miss  -> immunosuppression accelerates disease by factor k;
                                      residual unmasking time and post-detection survival
                                      both divided by k

## Endpoints (pre-specified, all reported regardless of outcome)

  E1. RMST-60: restricted mean survival time over a 60-month horizon, as a function of T.
  E2. P(alive at 60 months) — proxy for cure, as a function of T.
  E3. T* = argmax RMST-60, and the width of the region within 1 month of maximal RMST.
  E4. Fraction of sampled plausible parameter space in which T* is interior (T* not 0
      and not the maximum tested horizon).
  E5. Sensitivity of T* to each parameter individually.

## Positive controls — these must pass before any headline number is reported

  PC1 (external calibration anchor). With T at or near 0 — i.e. transplant everyone
      immediately, no filtering — the model must reproduce the observed median overall
      survival after heart transplantation for cardiac angiosarcoma of approximately
      9 months (versus 36 months for other cardiac sarcoma histologies).
      Calibration target is this external published number. It is fixed.
      The single free parameter tuned to hit it is k, the immunosuppression acceleration
      factor, which is the least directly measured quantity in the model.
      DO NOT tune any parameter toward producing an interior optimum.

  PC2 (degenerate case, no disease). With p_occ = 0, waiting buys nothing and costs
      device hazard. T* must equal 0. If it does not, the model is wrong.

  PC3 (degenerate case, no device hazard). With device hazard = 0, waiting is free and
      strictly informative. T* must go to the horizon. If it does not, the model is wrong.

  If PC1 fails, recalibrate against PC1 and report that it was recalibrated.
  If PC2 or PC3 fail, the model is structurally broken and no result is reportable.

## Decision rules fixed in advance

  R1. If an interior T* exists in the base case but survives in fewer than 50% of sampled
      plausible parameter draws, classify as PARAMETER-SENSITIVE and do not present it as
      a clinical recommendation.
  R2. If the RMST-60 difference between the best and worst bridge duration is smaller than
      1.0 month, classify as CLINICALLY NEGLIGIBLE regardless of statistical clarity —
      a mathematically real optimum that moves survival by two weeks is not actionable.
  R3. If the exponential (memoryless) unmasking assumption is doing the work — i.e. the
      result inverts under a Weibull unmasking distribution with shape > 1, which is the
      more biologically plausible growth model — classify as MODEL ARTEFACT.
  R4. Report every endpoint above regardless of whether it supports an interior optimum.

## What is NOT claimed

  - Not a clinical recommendation. No patient decision should be made from this.
  - Parameters are drawn from small retrospective series and case reports. Several
    (p_occ, k, unmasking kinetics) have no direct measurement at all and are stated
    as assumptions with swept ranges.
  - Quality of life on TAH support is not modelled and is substantial.
  - Systemic therapy during the bridge (KDR-directed, DDR-directed) is not modelled in
    the base case. It would change the answer and is a stated limitation.
  - The 9-month calibration anchor comes from a small retrospective transplant series.

## Final classification must be exactly one of

  A. Interior optimum, robust           (interior T*, survives >=50% of parameter space,
                                         RMST gap >= 1.0 month, survives Weibull check)
  B. Interior optimum, parameter-sensitive
  C. Corner solution (transplant immediately, or bridge indefinitely)
  D. Clinically negligible
  E. Model artefact

---

# Phase 2 pre-registration — systemic therapy during the bridge

Written after Phase 1 completed and before any Phase 2 code was executed.

## Why this phase exists

Phase 1 found the bridge to be pure attrition for curable patients: they wait, they pay
device hazard, they get nothing back. But the bridge window is the one period in this
disease when a patient has no allograft and therefore **no immunosuppression suppressing
their systemic therapy**. Cardiac angiosarcoma carries KDR alterations in 9/11 profiled
cases with a mechanistic account of their origin (POT1 loss -> unrepressed ATR-dependent
damage signalling -> somatic activating VEGF-pathway mutations), so KDR-directed and
DDR-directed therapy during the bridge is mechanistically indicated.

If therapy during the bridge works, waiting stops being dead time and becomes treatment
time. That is the only route by which a bridge could win on cure fraction rather than only
on organ arithmetic.

## The distinction being tested, which Phase 1 could not see

Therapy can act two ways, and the model separates them deliberately:

  CYTOSTATIC (parameter `tx_stasis` = s > 1). Slows growth of occult disease. Disease
      progress accrues at rate 1/s during the bridge. Does not remove disease.
  ERADICATING (parameter `tx_erad_rate` = per-month hazard). Some probability per unit
      time of actually clearing occult micrometastatic disease.

These have opposite effects on the *filter*, and that is the point of the phase. The bridge
is a diagnostic intervention. Slowing disease growth makes the diagnostic window less
informative — occult disease stays hidden longer and is more likely to be carried through
transplant — while the patient still pays full device hazard.

## Pre-registered predictions (direction stated in advance)

  F1. Pure cytostasis (s > 1, eradication = 0) will REDUCE five-year survival relative to
      s = 1 at any fixed bridge duration T > 0. Direction: HARMFUL. If cytostasis is
      neutral or beneficial, F1 is falsified and the "treating during a diagnostic window
      degrades the diagnostic" reasoning is wrong.
  F2. There exists a threshold eradication rate above which the five-year-survival-
      maximising bridge duration becomes > 0 — i.e. an efficacy bar that a bridge-period
      regimen must clear before waiting is justified at all. Report the bar as a number.
  F3. The F2 threshold rises with device hazard. Higher device hazard demands more
      effective therapy before waiting pays.

## Endpoints

  F-E1. T* by RMST-60 and by P(alive at 60 mo), over a grid of (eradication rate, stasis).
  F-E2. The eradication-rate threshold at which argmax P(alive60) leaves 0, in the base case.
  F-E3. That threshold as a function of device hazard.
  F-E4. Sign and size of the pure-cytostasis effect on P(alive60) at fixed T.
  F-E5. Whether cytostasis and eradication interact, or are separable.

## Positive control

  PC4. With eradication rate = 0 and stasis = 1, Phase 2 must reproduce the Phase 1 base
       case exactly (same T*, same RMST curve to Monte Carlo error). If it does not, the
       therapy machinery has changed the model rather than extended it, and Phase 2 is
       void.

## Decision rules

  F-R1. If the eradication rate required to move T* off zero exceeds what any systemic
        therapy plausibly achieves in metastatic angiosarcoma, report that the bar is
        UNREACHABLE and say so plainly. Do not soften it.
  F-R2. Report F1 direction as found, including if it is the opposite of predicted.
  F-R3. Therapy is modelled as active only during the bridge and stopping at transplant.
        This is a stated simplification, not a claim about practice.

---

# Phase 3 pre-registration — realistic response heterogeneity

Written after the Phase 3 literature audit (PROVENANCE.md) and **before any Phase 3 model
code was executed**. Nothing below is changed after seeing results.

## Why this phase exists

Phase 2's eradication model gave every patient a continuous per-month probability of cure.
That is biologically dubious: it implies universal partial curability. Real systemic
therapy in angiosarcoma works in a minority and fails outright in the rest. Phase 3 exists
to determine whether the Phase 2 result survives when that assumption is made realistic.

The objective is NOT to make Phase 2 survive. It is to find out whether it deserves to.

## The replacement model (fixed before running)

For each patient with occult micrometastatic disease:

    R_i ~ Bernoulli(p_R)                        latent responder status
    R_i = 0  ->  eradication probability is exactly ZERO, forever
    R_i = 1  ->  clearance achieved at C_i ~ Weibull(median m_clear, shape a_clear)
                 eradication occurs only if C_i < min(unmasking time, T)
    among patients who achieve clearance:
        durable with probability p_dur   -> cured
        otherwise                        -> disease re-emerges after G_i ~ Exp(median m_regrow)
                                            (acquired escape / residual disease)

Cytostasis is unchanged from Phase 2: disease progress accrues at rate 1/s during therapy.

Terminology held strictly distinct throughout, per the audit:
**tumour shrinkage != radiographic response != clearance of occult disease != cure.**

## Primary endpoint

  **P(alive at 60 months)**, as a function of bridge duration T.

Chosen because the strategy is curative-intent and because Phase 1 established that RMST
and the cure proxy diverge and can recommend opposite actions. The cure proxy is the
clinically meaningful one here.

## Primary question

Across externally admissible responder fractions and clearance kinetics (PROVENANCE.md),
does there remain a nonzero bridge duration T* > 0 that improves P(alive at 60 months)
relative to immediate transplantation, by a clinically meaningful margin?

## Secondary endpoints

  S1. RMST-60 as a function of T (reported for continuity with Phases 1-2).
  S2. Fraction of the cohort transplanted while still carrying occult disease.
  S3. Donor organs consumed per 100 patients.
  S4. The T* = 0 / T* > 0 phase boundary in (p_R, p_dur, m_clear, s, h_device) space.
  S5. Which parameter moves that boundary most strongly.
  S6. Whether the Phase 2 cytostasis x eradication synergy survives responder
      heterogeneity, and whether any surviving interaction is biological or an artefact
      of cytostasis changing WHEN occult disease becomes observable.
  S7. Required (p_R x p_dur) for benefit, compared against the observed/plausible region.

## Minimum clinically meaningful improvement (MCID)

  **5 percentage points absolute** in P(alive at 60 months) versus T = 0.

Fixed in advance. Justification: a bridge commits a patient to months of total artificial
heart support with substantial morbidity and quality-of-life cost that this model does not
capture. An absolute five-year survival gain below 5 points does not justify that.

## Robustness threshold

  The benefit must reach MCID in **at least 50%** of draws from the externally admissible
  parameter region defined in PROVENANCE.md.

## Decision rule for whether a finite bridge is supported

  D1. If the (p_R x p_dur) required to reach MCID lies ENTIRELY ABOVE the externally
      admissible region -> classification **A, bridge hypothesis fails**.
  D2. If MCID is reached in >=50% of admissible draws AND survives every Phase 3 attack
      (G-series below) -> classification **D, robust**.
  D3. If MCID is reached only above a quantifiable efficacy / device-performance boundary
      that lies partly inside the admissible region -> classification **C, conditional
      regime**, and the boundary must be reported as the deliverable, not a point estimate.
  D4. If the admissible region cannot be constrained enough by external evidence to
      distinguish A from C -> classification **B, model-indeterminate**. Report that
      rather than manufacturing precision.

## Stopping rules

  ST1. If any Phase 3 positive control fails, the main experiment is BLOCKED until fixed.
       Controls are asserted automatically and abort the run.
  ST2. If identifiability analysis shows the parameters are not constrainable from
       available external evidence, stop and report classification B. Do not select a
       favourite point in an unidentified region.
  ST3. If the required efficacy region and the evidence-supported region do not overlap,
       stop and report classification A. Do not soften.

## Positive controls (D) — asserted automatically, failure blocks the experiment

  PC6.  p_R = 0 must reproduce the Phase 2 no-eradication result exactly.
  PC7.  p_R = 1, immediate clearance (m_clear -> 0), p_dur = 1 must drive the carried-
        disease fraction to ~0 and P(alive60) toward the analytic ceiling
        (1 - periop_mort) x P(graft survives 60 mo) x P(no device event before T).
  PC8.  Device hazard = 0 must be monotone non-decreasing in T for the primary endpoint.
  PC9.  Very high device hazard must collapse T* to 0.
  PC10. Disabling cytostasis (s = 1) must reproduce the eradication-only branch exactly.
  PC11. Disabling eradication (p_R = 0) must reproduce the cytostasis-only branch exactly.
  PC12. Hand-calculation check: for p_R = 1, p_dur = 1, m_clear small, the cured fraction
        must match a closed-form expression to within Monte Carlo error.

PC12 exists specifically because the Phase 2 dead-code bug was caught by a hand
calculation and NOT by any pre-registered endpoint. Endpoint agreement is not evidence
that a mechanism fired.

## Attacks the surviving result must withstand (G)

  G1. Responder subpopulation only (the core Phase 3 change).
  G2. Cytostasis delays radiographic detectability (already in model; isolate its effect).
  G3. Non-exponential clearance kinetics (Weibull shape > 1).
  G4. Acquired escape / non-durable clearance (p_dur < 1).
  G5. Device hazard varying by centre quality.
  G6. Heterogeneous metastatic burden.
  G7. Separation of true biological synergy from decision-process interaction, by holding
      the observation process fixed while varying only the biology.

## What Phase 3 will NOT do

  - Will not model ATR-inhibitor efficacy. No defensible preclinical basis exists
    (PROVENANCE.md); simulating it would manufacture the result.
  - Will not treat KDR alteration as universal on the basis of an 11-patient cohort.
  - Will not tune any parameter toward a positive result.
  - Will not issue a clinical recommendation under any classification.

## Final classification must be exactly one of

  A — Bridge hypothesis fails. Realistic response heterogeneity removes the benefit.
  B — Model-indeterminate. Available data cannot constrain the parameters enough.
  C — Conditional bridge regime survives, above a quantified boundary.
  D — Robust bridge signal across the admissible region and all attacks.
