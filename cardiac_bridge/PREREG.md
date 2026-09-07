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
