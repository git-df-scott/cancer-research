"""
Part 2: does the Step-3 conclusion survive (a) non-exponential growth and
(b) an out-of-sample test that was NOT used to derive it?

The Step-3 conclusion under test:
    "The 0.49 yr ctDNA sojourn and the 3 yr MSCE sojourn are not in conflict.
     They imply detection thresholds ~2.3 cm and ~0.5 cm respectively -- i.e.
     they were never measuring the same state."

Threat 1: it assumed exponential growth. PDAC may be Gompertzian.
Threat 2: it is so far only self-consistent, not validated.
"""
import numpy as np
from scipy import integrate, optimize

DPY = 365.25
V0 = 1e-9
vol = lambda d: (np.pi / 6.0) * d ** 3
dia = lambda v: (6.0 * v / np.pi) ** (1.0 / 3.0)

SIGMA = 0.8326
MU = -4.9027
E_INV_G = np.exp(-MU + SIGMA ** 2 / 2)          # 190.4 d
D_PRES, V_PRES = 3.1, vol(3.1)

print("=" * 78)
print("STEP 6  Growth-law robustness: exponential vs Gompertz")
print("=" * 78)
print("""  Gompertz:  dV/dt = a*V*ln(K/V).  Growth decelerates near carrying
  capacity K, so the time to cross the LAST few doublings is inflated
  relative to exponential. If PDAC is Gompertzian, a fixed observed sojourn
  implies a detection threshold CLOSER to presentation size than the
  exponential calculation says -- which would strengthen, not weaken, the
  claim about ctDNA. Check the magnitude.\n""")


def gompertz_time(v_a, v_b, a, K):
    """Time to grow v_a -> v_b under dV/dt = a V ln(K/V)."""
    # solution: ln(K/V(t)) = ln(K/V_a) exp(-a t)  =>  t = ln( ln(K/v_a)/ln(K/v_b) ) / a
    return np.log(np.log(K / v_a) / np.log(K / v_b)) / a


def implied_threshold_gompertz(sojourn_yr, K_mult):
    """K = K_mult * V_pres. Calibrate 'a' so the Gompertz curve has the same
    mean VDT as the exponential fit at the midpoint size, then invert."""
    K = K_mult * V_PRES
    # calibrate a: match instantaneous doubling time at V = V_PRES/8 (~1.55 cm)
    v_ref = V_PRES / 8.0
    g_ref = np.log(2) / (np.log(2) * E_INV_G / 1.0)   # 1/E_inv_g -> per day
    g_ref = 1.0 / E_INV_G
    a = g_ref / np.log(K / v_ref)
    f = lambda logv: gompertz_time(np.exp(logv), V_PRES, a, K) - sojourn_yr * DPY
    lo, hi = np.log(V0), np.log(V_PRES * 0.999)
    try:
        logv = optimize.brentq(f, lo, hi)
    except ValueError:
        return np.nan
    return dia(np.exp(logv))


def implied_threshold_exp(sojourn_yr):
    return dia(V_PRES / np.exp(sojourn_yr * DPY / E_INV_G))


print(f"  {'sojourn':>10}{'exponential':>16}{'Gompertz K=2x':>16}"
      f"{'K=10x':>12}{'K=100x':>12}")
print("  " + "-" * 66)
for s in [0.49, 3.00, 8.93]:
    row = [implied_threshold_exp(s)] + [implied_threshold_gompertz(s, k)
                                        for k in (2.0, 10.0, 100.0)]
    print(f"  {s:>8.2f} yr" + "".join(f"{v:>13.2f} cm" if np.isfinite(v)
                                      else f"{'--':>16}" for v in row[:1])
          + "".join(f"{v:>13.2f} cm" if np.isfinite(v) else f"{'--':>16}"
                    for v in row[1:]))
print()
print("  Conclusion is ORDER-ROBUST: under every growth law tested, the ctDNA")
print("  anchor implies a threshold within ~1 cm of presentation size, and the")
print("  MSCE anchor implies one far below imaging resolution.")
print()


# ----------------------------------------------------------------------------
print("=" * 78)
print("STEP 7  OUT-OF-SAMPLE TEST (anchor never used in the derivation)")
print("=" * 78)
print("""  Derivation used A1 (ctDNA sojourn), A4 (VDT), A5 (size at dx) only.
  It PREDICTS a ctDNA limit of detection at ~2.3 cm diameter.
  That prediction implies a stage-specific ctDNA sensitivity, which has been
  measured independently and was NOT used above. This is the falsification
  test.\n""")

V_CTDNA = V_PRES / np.exp(0.49 * DPY / E_INV_G)
D_CTDNA = dia(V_CTDNA)

# Size distribution of clinically-presenting stage I/II PDAC.
# AJCC 8th: T1 <=2cm, T2 >2-4cm, T3 >4cm. Stage I/II resected series report
# median ~2.5-3.0 cm. Model as lognormal, median 2.7 cm, GSD 1.45.
rng = np.random.default_rng(1)
d_stage12 = np.exp(rng.normal(np.log(2.7), np.log(1.45), 200000))
d_stage12 = d_stage12[d_stage12 <= 4.0]          # stage I/II size ceiling
d_stage1 = d_stage12[d_stage12 <= 2.0]           # T1 = stage IA

for label, arr in [("stage I/II (<=4 cm)", d_stage12),
                   ("stage IA  (<=2 cm)", d_stage1)]:
    pred = np.mean(arr >= D_CTDNA)
    print(f"  predicted ctDNA sensitivity, {label:<22} {100*pred:5.1f} %")

print()
print("  MEASURED, independently:")
print("    KRAS-only digital NGS, stage I/II          31 %")
print("    multi-analyte ctDNA+protein, stage I       62 %")
print("    multi-analyte ctDNA+protein, stage II      56 %")
print("""
  READ THIS HONESTLY -- it is a SPLIT result, not a clean pass.

  PASS at stage I/II: predicted 62.7% vs measured 31% (weak single-locus
  assay) to 56-62% (strong multi-analyte assay). The prediction lands at the
  top of the measured range. For a number derived from three unrelated
  anchors with zero tuning, that is a real hit.

  FAIL at stage IA: predicted 0.0% vs measured non-zero. This is a genuine
  falsification of the HARD-THRESHOLD detection assumption, not of the
  underlying claim. A step function at 2.27 cm makes every <=2 cm tumour
  invisible by construction, which is plainly wrong -- ctDNA shedding varies
  per tumour and the assay has a soft ROC, so detection must be probabilistic
  in volume, e.g. P(detect) = 1 - exp(-lambda * V), not a step.

  CONSEQUENCE: the model needs a soft detection function, which adds ONE
  parameter (lambda) per modality. That changes the identifiability budget --
  re-run Step 5 with the soft version before building anything. The stage-IA
  sensitivity then becomes a second usable anchor rather than a degenerate
  prediction, so it may pay for itself. This must be settled first.
""")

# ----------------------------------------------------------------------------
print("=" * 78)
print("STEP 8  Closing the one flat direction")
print("=" * 78)
print("""  Step 5 found rank 5 of 6, with the unconstrained direction dominated by
  V_img (loading -0.84). That is expected: A2 (MSCE) was used to set V_img,
  but the MSCE 'preclinical screen-detectable' compartment is an incidence-fit
  abstraction, not an imaging state -- it implies a 0.46 cm threshold, well
  below what CT resolves for PDAC.

  So V_img needs its own external anchor. The natural one:

    A6  EUS/MRI detection limit for solid pancreatic lesions in surveillance,
        ~0.5-1.0 cm, plus the size distribution of CAPS screen-detected
        lesions (reported in the CAPS series; needs extraction from the
        paper, not available from abstract-level search).

  With A6 supplied, the parameter count (6) equals the anchor count (6) and
  the flat direction closes. Condition number at the trial point was 87 --
  well conditioned, not marginal. Nothing else in the spectrum is near zero.
""")
print("=" * 78)
print("VERDICT")
print("=" * 78)
print("""  IDENTIFIABLE, conditional on one more anchor (A6), and only because the
  two DIRECT measurements (A4 volume doubling times, A5 size at diagnosis)
  are brought in. With the three original anchors alone (A1, A2, A3) the
  model would have been badly under-determined -- 6 parameters, 3 constraints.

  The reconciliation project is worth running.

  But note what already fell out of the algebra, before any simulation:
  the paradox is substantially DISSOLVED rather than pending. See REPORT.""")
