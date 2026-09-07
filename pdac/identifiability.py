"""
Identifiability analysis for a minimal PDAC natural-history model.

QUESTION THIS SCRIPT ANSWERS (and only this):
    Do the published anchors actually CONSTRAIN the parameters of a minimal
    PDAC natural-history model, or would fitting them be an unfalsifiable
    curve-fit with more freedom than data?

Run BEFORE building any simulation. If the answer is "under-identified",
the reconciliation project as scoped is not worth running.

------------------------------------------------------------------------------
MODEL M0 (deliberately minimal)
------------------------------------------------------------------------------
A single invasive PDAC clone, founder cell at t=0.

    V(t) = v0 * exp(g*t)                     volume, cm^3
    v0   = 1e-9 cm^3                         one cell (~1000 um^3)
    g ~ LogNormal(mu, sigma)                 per-tumour growth rate, /day

Clinical presentation when V crosses V_pres.
Modality k detects when V crosses V_k  (k in {img, ctdna}).

    S_k = time from V_k to V_pres = ln(V_pres / V_k) / g      "sojourn"

Metastatic seeding: hazard  m * V(t)  (seeding rate proportional to burden).

Parameters theta = (mu, sigma, V_pres, V_img, V_ctdna, m).

------------------------------------------------------------------------------
ANCHORS (all external, all published, none tunable)
------------------------------------------------------------------------------
A1  ctDNA sojourn                 0.49 yr  [0.26, 0.88]
        Hubbell et al., Cancer 2026 (CPS-3 / CCGA3), shortest of all cancers.
A2  MSCE preclinical sojourn      ~3 yr
        Luebeck MSCE fit to SEER incidence.
A3  CAPS5 stage-I fraction        7/9 = 0.778 under ANNUAL imaging surveillance
        Multicenter CAPS study, JCO 2022.
A4  volume doubling time          mean 132 d, SD 132 d  (range 20-977)
        Serial-CT series, J Gastrointest Cancer 2016.
A5  tumour size at diagnosis      mean ~3.1 cm diameter
        SEER / standard PDAC presentation.

A4 and A5 are DIRECT measurements, not model fits. That distinction is the
whole ballgame -- see the verdict at the bottom.
"""

import numpy as np
from scipy import stats, optimize

DAYS_PER_YEAR = 365.25
V0 = 1e-9                       # cm^3, one cell


def vol_from_diam(d_cm):
    return (np.pi / 6.0) * d_cm ** 3


def diam_from_vol(v):
    return (6.0 * v / np.pi) ** (1.0 / 3.0)


# ----------------------------------------------------------------------------
# STEP 1. Pin the growth distribution from A4 alone (direct measurement).
# ----------------------------------------------------------------------------
# VDT = ln2 / g.  If g ~ LogNormal(mu, sigma) then VDT ~ LogNormal(ln(ln2)-mu, sigma).
#   mean(VDT) = ln2 * exp(-mu + sigma^2/2)
#   CV(VDT)   = sqrt(exp(sigma^2) - 1)      <- depends on sigma ONLY

VDT_MEAN, VDT_SD = 132.0, 132.0
cv = VDT_SD / VDT_MEAN
sigma_hat = np.sqrt(np.log(1.0 + cv ** 2))
mu_hat = -(np.log(VDT_MEAN / np.log(2.0)) - sigma_hat ** 2 / 2.0)

# Moments of 1/g, which is what every sojourn depends on.
#   1/g ~ LogNormal(-mu, sigma)
E_inv_g = np.exp(-mu_hat + sigma_hat ** 2 / 2.0)                    # days
E_inv_g2 = np.exp(-2 * mu_hat + 2 * sigma_hat ** 2)                 # E[(1/g)^2]

print("=" * 78)
print("STEP 1  Growth distribution, from A4 (direct VDT measurement) alone")
print("=" * 78)
print(f"  CV of VDT                    {cv:.3f}   (50-fold spread, 20-977 d)")
print(f"  sigma  (log-scale spread)    {sigma_hat:.4f}   <- fixed by CV alone")
print(f"  mu     (log-scale location)  {mu_hat:.4f}   <- fixed by mean given sigma")
print(f"  median VDT                   {np.log(2)*np.exp(-mu_hat):.1f} d")
print(f"  E[1/g]                       {E_inv_g:.1f} d  = {E_inv_g/DAYS_PER_YEAR:.4f} yr")
print(f"  1/E[g]                       {1.0/np.exp(mu_hat+sigma_hat**2/2)/1:.1f} d"
      f"   <- NOT the same; using it biases sojourn by "
      f"{100*(E_inv_g/(1.0/np.exp(mu_hat+sigma_hat**2/2))-1):.0f}%")
print()
print("  => mu and sigma are IDENTIFIED by A4 without touching any sojourn anchor.")
print()


# ----------------------------------------------------------------------------
# STEP 2. Structural result: what the sojourn anchors can and cannot see.
# ----------------------------------------------------------------------------
# S_k = ln(V_pres/V_k) / g.  Therefore
#     E[S_k] = ln(V_pres/V_k) * E[1/g]
# Two consequences, both structural (no numbers needed):
#
#   (a) Thresholds enter ONLY as the ratio V_pres/V_k. Absolute volumes are
#       NOT identifiable from sojourn data. Need A5 to set the scale.
#
#   (b) The RATIO of two sojourns is completely growth-free:
#           E[S_img] / E[S_ctdna] = ln(V_pres/V_img) / ln(V_pres/V_ctdna)
#       E[1/g] cancels. So A1 and A2 together impose a constraint that CANNOT
#       be absorbed by retuning growth. That makes it a TEST, not a fit.

print("=" * 78)
print("STEP 2  Structural identifiability (analytic)")
print("=" * 78)
print("  E[S_k] = ln(V_pres/V_k) * E[1/g]")
print("  (a) thresholds enter only as ratios  -> absolute volumes need A5")
print("  (b) E[S_img]/E[S_ctdna] is growth-INDEPENDENT -> A1+A2 is a test, not a fit")
print()

# Set the volume scale from A5.
D_PRES = 3.1                                    # cm, mean size at diagnosis
V_PRES = vol_from_diam(D_PRES)

# What detection threshold does each sojourn anchor IMPLY?
def implied_threshold(sojourn_years):
    """Invert E[S] = ln(V_pres/V_k) * E[1/g] for V_k, then report diameter."""
    log_ratio = sojourn_years * DAYS_PER_YEAR / E_inv_g
    v_k = V_PRES / np.exp(log_ratio)
    return v_k, diam_from_vol(v_k), log_ratio

print("=" * 78)
print("STEP 3  Detection thresholds implied by each sojourn anchor")
print("=" * 78)
print(f"  volume scale set by A5: presentation at {D_PRES} cm = {V_PRES:.2f} cm^3")
print()
print(f"  {'anchor':<34}{'sojourn':>10}{'implied threshold':>22}")
print("  " + "-" * 66)
for name, s in [("A1  ctDNA (Hubbell 2026)", 0.49),
                ("A1  ctDNA, lower CI", 0.26),
                ("A1  ctDNA, upper CI", 0.88),
                ("A2  MSCE preclinical (Luebeck)", 3.00),
                ("    proteomic relay CTHRC1", 8.93)]:
    v, d, lr = implied_threshold(s)
    print(f"  {name:<34}{s:>8.2f} yr{d:>16.2f} cm dia")
print()
ratio_obs = 3.00 / 0.49
print(f"  Growth-free structural constraint from A1+A2:")
print(f"    E[S_MSCE]/E[S_ctDNA] = {ratio_obs:.2f}")
print(f"    must equal ln(V_pres/V_MSCE)/ln(V_pres/V_ctDNA)  -- no growth freedom")
print()


# ----------------------------------------------------------------------------
# STEP 4. Length bias. With CV=1 this is NOT a footnote.
# ----------------------------------------------------------------------------
# A screen-detected prevalent case is sampled with probability proportional to
# its sojourn length. The length-biased mean is E[S^2]/E[S].
#   E[S^2] = ln(ratio)^2 * E[(1/g)^2]

def sojourn_moments(log_ratio):
    e1 = log_ratio * E_inv_g
    e2 = log_ratio ** 2 * E_inv_g2
    return e1, e2

print("=" * 78)
print("STEP 4  Length-biased sampling (matters a lot at CV=1)")
print("=" * 78)
for name, s in [("ctDNA (A1)", 0.49), ("MSCE (A2)", 3.00)]:
    _, _, lr = implied_threshold(s)
    e1, e2 = sojourn_moments(lr)
    lb = e2 / e1
    print(f"  {name:<14} incident mean {e1/DAYS_PER_YEAR:5.2f} yr   "
          f"length-biased mean {lb/DAYS_PER_YEAR:5.2f} yr   "
          f"inflation x{lb/e1:.2f}")
print()
print(f"  Inflation factor E[S^2]/E[S]^2 = exp(sigma^2) = {np.exp(sigma_hat**2):.2f}")
print("  => If A1 and A2 are not BOTH incident-mean or BOTH length-biased,")
print("     they differ by 2x for purely methodological reasons.")
print()


# ----------------------------------------------------------------------------
# STEP 5. Practical identifiability: Jacobian rank / singular values.
# ----------------------------------------------------------------------------
# Predict the anchor vector from theta, take the Jacobian in log-parameters,
# and look at the singular value spectrum. Near-zero singular values = flat
# directions = parameters the anchors cannot separate.

def predict(theta_log, annual_interval_days=365.25, n_mc=40000, seed=0):
    """theta_log = log of (mu_neg, sigma, V_pres, V_img, V_ctdna, m).
    mu is passed as -mu so it can live in log space."""
    neg_mu, sigma, v_pres, v_img, v_ctdna, m = np.exp(theta_log)
    mu = -neg_mu
    rng = np.random.default_rng(seed)
    g = np.exp(rng.normal(mu, sigma, n_mc))

    t_of = lambda v: np.log(v / V0) / g
    t_pres, t_img, t_ct = t_of(v_pres), t_of(v_img), t_of(v_ctdna)

    s_ct = np.mean(t_pres - t_ct) / DAYS_PER_YEAR
    s_img = np.mean(t_pres - t_img) / DAYS_PER_YEAR

    # metastatic seeding: hazard m*V(t); cumulative hazard to volume v is
    # integral_0^{t(v)} m*V0*exp(g t) dt = m*(v - V0)/g  ~= m*v/g
    # Annual imaging surveillance: detection uniformly within a year of
    # crossing V_img (first scan after the threshold is crossed).
    delay = rng.uniform(0.0, annual_interval_days, n_mc)
    v_at_detect = np.minimum(v_img * np.exp(g * delay), v_pres)
    cumhaz = m * v_at_detect / g
    stage1_frac = np.mean(np.exp(-cumhaz))          # P(no seeding yet)

    d_pres = diam_from_vol(v_pres)
    vdt_mean = np.mean(np.log(2) / g)
    return np.array([s_ct, s_img, stage1_frac, vdt_mean, d_pres])


ANCHOR_NAMES = ["A1 ctDNA sojourn (yr)", "A2 MSCE sojourn (yr)",
                "A3 CAPS stage-I frac", "A4 mean VDT (d)", "A5 size at dx (cm)"]
PARAM_NAMES = ["-mu", "sigma", "V_pres", "V_img", "V_ctdna", "m"]

# A plausible point that roughly hits the anchors (starting point, not a fit).
v_ct_star = implied_threshold(0.49)[0]
v_img_star = implied_threshold(3.00)[0]
theta0 = np.log(np.array([-mu_hat, sigma_hat, V_PRES, v_img_star, v_ct_star, 0.02]))

# tune m to land near 0.778 stage-I
def _f(logm):
    th = theta0.copy(); th[5] = logm
    return predict(th)[2] - 0.778
try:
    theta0[5] = optimize.brentq(_f, np.log(1e-4), np.log(1e3))
except ValueError:
    pass

pred0 = predict(theta0)
print("=" * 78)
print("STEP 5  Practical identifiability: Jacobian singular spectrum")
print("=" * 78)
print("  anchor predictions at the trial point:")
for n, p, t in zip(ANCHOR_NAMES, pred0, [0.49, 3.00, 0.778, 132.0, 3.1]):
    print(f"    {n:<26} pred {p:9.3f}   target {t:8.3f}")
print()

# Jacobian by central differences in log-parameter space (same seed = common
# random numbers, so the derivative is not swamped by MC noise).
eps = 1e-3
J = np.zeros((len(pred0), len(theta0)))
for j in range(len(theta0)):
    tp, tm = theta0.copy(), theta0.copy()
    tp[j] += eps; tm[j] -= eps
    J[:, j] = (predict(tp) - predict(tm)) / (2 * eps)

# scale rows so each anchor contributes on a comparable relative footing
Jn = J / np.maximum(np.abs(pred0)[:, None], 1e-12)
U, S, Vt = np.linalg.svd(Jn)
print("  singular values of the scaled Jacobian (6 params, 5 anchors):")
for i, s in enumerate(S):
    print(f"    s{i+1} = {s:10.4f}")
print(f"  condition number (s1/s_last) = {S[0]/S[-1]:.3e}")
print(f"  numerical rank (tol 1e-8 * s1) = {np.sum(S > 1e-8*S[0])} of {len(theta0)}")
print()
print("  worst-constrained direction (right singular vector of smallest s):")
worst = Vt[np.argmin(S)]
for n, w in sorted(zip(PARAM_NAMES, worst), key=lambda x: -abs(x[1])):
    print(f"    {n:<10}{w:+.4f}")
print()
