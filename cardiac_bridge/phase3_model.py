"""Phase 3: latent responder model replacing Phase 2's constant eradication hazard.

Pre-registered in PREREG.md (Phase 3 section) before this file was executed.

Phase 2 gave every patient a continuous per-month cure probability. Phase 3 replaces it:
a latent responder fraction, with non-responders having exactly zero eradication
probability forever, responders clearing on their own kinetics, and clearance being
durable only some of the time.

Held strictly distinct: shrinkage != radiographic response != clearance != cure.
"""
import numpy as np
from dataclasses import dataclass, replace
from bridge_model import HORIZON, CAL_TARGET, weibull_median


@dataclass(frozen=True)
class P3:
    # --- inherited Phase 1 structure ---
    p_occ: float = 0.65
    med_unmask: float = 6.0
    unmask_shape: float = 1.0
    dev_rate: float = 0.025
    med_met: float = 6.0
    k: float = 2.8
    periop_mort: float = 0.08
    med_graft: float = 144.0
    q: float = 1.0
    tx_stasis: float = 1.0          # cytostasis, unchanged from Phase 2
    # --- Phase 3 responder model ---
    p_R: float = 0.0                # latent responder fraction; 0 -> no eradication at all
    m_clear: float = 4.0            # median months to clearance, RESPONDERS ONLY
    a_clear: float = 1.0            # Weibull shape of clearance time
    p_dur: float = 0.5              # P(clearance is durable | cleared)
    m_regrow: float = 6.0           # median months to re-emergence if clearance not durable
    # G7: separates the BIOLOGICAL effect of cytostasis (more time for clearance to act)
    # from the DECISION-PROCESS effect (disease stays radiographically hidden, so more
    # patients remain transplant-eligible). In reality these are the same drug effect;
    # splitting them is diagnostic, not a claim about biology. None -> tied to tx_stasis.
    s_detect: float = None          # stasis applied to the DETECTION clock only


def draw_latents3(p: P3, n, rng):
    occult = rng.random(n) < p.p_occ
    U = np.where(occult, weibull_median(rng, p.med_unmask, p.unmask_shape, n), np.inf)
    D = rng.exponential(1.0 / p.dev_rate, n)
    met = rng.exponential(p.med_met / np.log(2.0), n)
    met_after_tx = rng.exponential(p.med_met / np.log(2.0), n)
    dev_resid = rng.exponential(1.0 / p.dev_rate, n)
    graft = rng.exponential(p.med_graft / np.log(2.0), n)
    periop_death = rng.random(n) < p.periop_mort
    # Phase 3 latents drawn as UNIFORMS so that p_R, m_clear, a_clear and p_dur can all be
    # varied at evaluation time without redrawing. Baking a rate into the draw is exactly
    # the Phase 2 dead-code bug; PC6/PC12 exist to catch a recurrence.
    R_u = rng.random(n)             # -> responder if R_u < p_R
    C_u = rng.random(n)             # -> clearance time by inverse-CDF at evaluation time
    dur_u = rng.random(n)           # -> durable if dur_u < p_dur
    G_unit = rng.exponential(1.0, n)
    return dict(occult=occult, U=U, D=D, met=met, met_after_tx=met_after_tx,
                dev_resid=dev_resid, graft=graft, periop_death=periop_death,
                R_u=R_u, C_u=C_u, dur_u=dur_u, G_unit=G_unit)


def _clearance_time(p: P3, L):
    """Weibull inverse-CDF from stored uniforms, so median and shape vary at eval time."""
    scale = p.m_clear / (np.log(2.0) ** (1.0 / p.a_clear))
    return scale * (-np.log1p(-L["C_u"])) ** (1.0 / p.a_clear)


def disease_state(T, p: P3, L):
    """Resolve each patient's occult-disease trajectory under bridge therapy.

    Returns (U_eff, cured, s_eff):
      U_eff  time at which occult disease becomes radiographically detectable (inf if none)
      cured  durable clearance achieved during the bridge
      s_eff  therapy-progress divisor applying to residual disease after transplant
    """
    s = p.tx_stasis
    s_det = s if p.s_detect is None else p.s_detect
    U_bio = L["U"] * s                          # biological clock: when clearance must beat
    U_obs = L["U"] * s_det                      # observation clock: when it becomes visible
    responder = L["occult"] & (L["R_u"] < p.p_R)
    C = _clearance_time(p, L)
    cleared = responder & (C < np.minimum(U_bio, T))        # must clear before unmasking/tx
    durable = cleared & (L["dur_u"] < p.p_dur)
    transient = cleared & ~durable

    # non-durable clearance buys a delay, not a cure: disease re-emerges at C + G
    G = L["G_unit"] * (p.m_regrow / np.log(2.0))
    U_eff = np.where(durable, np.inf,
             np.where(transient, C + G, U_obs))
    s_eff = np.where(transient, 1.0, s_det)     # re-emergence delay already encodes therapy
    return U_eff, durable, s_eff


def survival3(T, p: P3, L):
    D = L["D"]
    U_eff, cured, s_eff = disease_state(T, p, L)
    has_disease = L["occult"] & ~cured

    U_det = np.where(has_disease, np.ceil(U_eff / p.q) * p.q, np.inf)
    transplanted = (D > T) & (U_det > T)
    dev_first = (~transplanted) & (D <= U_det)
    det_first = (~transplanted) & (D > U_det)

    surv = np.empty_like(D)
    surv[dev_first] = D[dev_first]
    surv[det_first] = U_det[det_first] + np.minimum(L["met"][det_first],
                                                    L["dev_resid"][det_first])

    tx = transplanted
    occ_tx = tx & has_disease & (U_eff > T)     # transplanted still carrying disease
    free_tx = tx & ~occ_tx
    surv[free_tx] = T + L["graft"][free_tx]

    resid = (U_eff[occ_tx] - T) / (s_eff[occ_tx] * p.k)
    post = L["met_after_tx"][occ_tx] / p.k
    surv[occ_tx] = T + np.minimum(resid + post, L["graft"][occ_tx])

    surv[tx & L["periop_death"]] = T
    return surv


def carried3(T, p: P3, L):
    D = L["D"]
    U_eff, cured, _ = disease_state(T, p, L)
    has_disease = L["occult"] & ~cured
    U_det = np.where(has_disease, np.ceil(U_eff / p.q) * p.q, np.inf)
    tx = (D > T) & (U_det > T)
    return float(np.mean(tx & has_disease & (U_eff > T)))


def tx_rate3(T, p: P3, L):
    U_eff, cured, _ = disease_state(T, p, L)
    has_disease = L["occult"] & ~cured
    U_det = np.where(has_disease, np.ceil(U_eff / p.q) * p.q, np.inf)
    return float(np.mean((L["D"] > T) & (U_det > T)))


def metrics3(T, p: P3, L):
    s = survival3(T, p, L)
    return dict(rmst=float(np.mean(np.minimum(s, HORIZON))),
                alive60=float(np.mean(s >= HORIZON)),
                median=float(np.median(s)),
                tx_rate=tx_rate3(T, p, L),
                carried=carried3(T, p, L))


def calibrate_k3(p: P3, L, target=CAL_TARGET, lo=1.0, hi=15.0, tol=0.02, iters=22):
    """Same external anchor as Phases 1-2: median OS ~9 mo transplanting with no bridge.
    At T=0 there is no bridge therapy exposure, so the anchor is independent of p_R."""
    f = lambda kk: float(np.median(survival3(0.0, replace(p, k=kk), L))) - target
    flo, fhi = f(lo), f(hi)
    if flo * fhi > 0:
        return None
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if abs(fm) < tol:
            return mid
        if flo * fm <= 0:
            hi, fhi = mid, fm
        else:
            lo, flo = mid, fm
    return 0.5 * (lo + hi)


def analytic_ceiling(T, p: P3):
    """Closed form for PC12: if every occult patient is cured instantly at bridge start,
    P(alive at 60 mo) = P(no device event before T) x P(survive transplant) x P(graft>60-T).
    No Monte Carlo involved."""
    return (np.exp(-p.dev_rate * T) * (1 - p.periop_mort)
            * 0.5 ** ((HORIZON - T) / p.med_graft))
