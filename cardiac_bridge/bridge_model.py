"""
Optimal bridge duration on a total artificial heart after total cardiectomy for
primary cardiac sarcoma.

Pre-registered in PREREG.md. Endpoints E1-E5, positive controls PC1-PC3 and
decision rules R1-R4 were fixed before this file was executed.

Clock starts at cardiectomy. All times in months.
"""
import numpy as np
from dataclasses import dataclass, replace

HORIZON = 60.0          # months, restricted-mean horizon (E1)
CAL_TARGET = 9.0        # PC1: median OS after transplant for cardiac angiosarcoma


@dataclass(frozen=True)
class Params:
    p_occ: float = 0.65        # P(occult micrometastatic disease at cardiectomy)
    med_unmask: float = 6.0    # median months to radiographic unmasking, no immunosuppression
    unmask_shape: float = 1.0  # Weibull shape; 1.0 = exponential (memoryless)
    dev_rate: float = 0.025    # /month, fatal or transplant-precluding device event
    med_met: float = 6.0       # median survival from detected metastatic disease, no IS
    k: float = 2.8             # immunosuppression acceleration factor (calibrated to PC1)
    periop_mort: float = 0.08  # perioperative mortality at transplant
    med_graft: float = 144.0   # median graft survival if truly disease-free
    q: float = 1.0             # surveillance imaging interval
    # --- Phase 2: systemic therapy during the bridge (inactive at these defaults) ---
    tx_erad_rate: float = 0.0  # /month hazard of clearing occult micrometastatic disease
    tx_stasis: float = 1.0     # >1 slows occult disease growth; disease progress accrues at 1/s


def weibull_median(rng, median, shape, size):
    """Weibull draws with a specified median. shape=1 gives the exponential."""
    scale = median / (np.log(2.0) ** (1.0 / shape))
    return scale * rng.weibull(shape, size)


def draw_latents(p: Params, n, rng):
    """Latent patient state, drawn once and reused across every candidate T.

    Common random numbers across T are what make the RMST(T) curve smooth enough
    for argmax to mean anything.
    """
    occult = rng.random(n) < p.p_occ
    U = np.where(occult, weibull_median(rng, p.med_unmask, p.unmask_shape, n), np.inf)
    U_det = np.where(np.isfinite(U), np.ceil(U / p.q) * p.q, np.inf)   # detected at next scan
    D = rng.exponential(1.0 / p.dev_rate, n)                           # device event
    met = rng.exponential(p.med_met / np.log(2.0), n)                  # metastatic survival, no IS
    met_after_tx = rng.exponential(p.med_met / np.log(2.0), n)         # independent draw post-tx
    dev_resid = rng.exponential(1.0 / p.dev_rate, n)                   # residual device life, branch B
    graft = rng.exponential(p.med_graft / np.log(2.0), n)
    periop_death = rng.random(n) < p.periop_mort
    # Drawn LAST so every Phase 1 variable keeps its exact random stream (PC4).
    # UNIT exponential: the eradication rate is applied in _eradicated() at evaluation
    # time. Baking the rate in here made tx_erad_rate a silent no-op, because latents are
    # drawn once and the rate is then varied via dataclasses.replace(). Caught by PC5.
    E_unit = rng.exponential(1.0, n)
    return dict(occult=occult, U=U, U_det=U_det, D=D, met=met,
                met_after_tx=met_after_tx, dev_resid=dev_resid,
                graft=graft, periop_death=periop_death, E_unit=E_unit)


def _eradicated(T, p: Params, L, U_bridge):
    """Occult disease cleared by bridge therapy before it unmasks or the bridge ends."""
    if p.tx_erad_rate <= 0:
        return np.zeros(len(L["D"]), dtype=bool)
    E = L["E_unit"] / p.tx_erad_rate
    return L["occult"] & (E < np.minimum(U_bridge, T))


def survival(T, p: Params, L):
    """Overall survival from cardiectomy under bridge duration T.

    Bridge therapy (Phase 2) acts only during [0, T] and stops at transplant.
    Cytostasis stretches occult disease progress: it accrues at rate 1/s while on
    therapy, so unmasking that would occur at U occurs at U*s instead.
    Eradication clears occult disease at hazard tx_erad_rate while on therapy.
    At s = 1 and rate = 0 this is identical to Phase 1.
    """
    U, D, s = L["U"], L["D"], p.tx_stasis

    U_bridge = U * s                                    # unmasking time under therapy
    eradicated = _eradicated(T, p, L, U_bridge)
    U_det = np.where(L["occult"] & ~eradicated,
                     np.ceil(U_bridge / p.q) * p.q, np.inf)

    transplanted = (D > T) & (U_det > T)
    dev_first    = (~transplanted) & (D <= U_det)
    det_first    = (~transplanted) & (D > U_det)

    surv = np.empty_like(D)
    surv[dev_first] = D[dev_first]
    surv[det_first] = U_det[det_first] + np.minimum(L["met"][det_first],
                                                    L["dev_resid"][det_first])

    tx = transplanted
    occ_tx = tx & L["occult"] & ~eradicated & (U_bridge > T)   # disease carried through
    free_tx = tx & ~occ_tx

    surv[free_tx] = T + L["graft"][free_tx]

    # residual natural progress left at transplant, then accelerated by immunosuppression
    resid = (U[occ_tx] - T / s) / p.k
    post = L["met_after_tx"][occ_tx] / p.k
    surv[occ_tx] = T + np.minimum(resid + post, L["graft"][occ_tx])

    surv[tx & L["periop_death"]] = T
    return surv


def carried_disease(T, p: Params, L):
    """Fraction of the whole cohort transplanted while still carrying occult disease."""
    U, s = L["U"], p.tx_stasis
    U_bridge = U * s
    eradicated = _eradicated(T, p, L, U_bridge)
    U_det = np.where(L["occult"] & ~eradicated, np.ceil(U_bridge / p.q) * p.q, np.inf)
    tx = (L["D"] > T) & (U_det > T)
    return float(np.mean(tx & L["occult"] & ~eradicated & (U_bridge > T)))


def transplant_rate(T, p: Params, L):
    U_bridge = L["U"] * p.tx_stasis
    eradicated = _eradicated(T, p, L, U_bridge)
    U_det = np.where(L["occult"] & ~eradicated, np.ceil(U_bridge / p.q) * p.q, np.inf)
    return float(np.mean((L["D"] > T) & (U_det > T)))


def metrics(T, p: Params, L):
    s = survival(T, p, L)
    return dict(
        rmst=float(np.mean(np.minimum(s, HORIZON))),     # E1
        alive60=float(np.mean(s >= HORIZON)),            # E2
        median=float(np.median(s)),
        tx_rate=transplant_rate(T, p, L),
        futile_tx=carried_disease(T, p, L),
    )


def median_os_at_T0(p: Params, L):
    """PC1 observable: median OS when everyone is transplanted immediately."""
    return float(np.median(survival(0.0, p, L)))


def calibrate_k(p: Params, L, target=CAL_TARGET, lo=1.0, hi=15.0, tol=0.02, iters=22):
    """Solve for the immunosuppression acceleration factor k that reproduces the
    external anchor: median OS ~= 9 months when transplanting without a bridge.

    k is the only parameter tuned, and it is tuned to a fixed published number,
    never toward the shape of the T curve. Returns None if the anchor is
    unreachable anywhere in [lo, hi] for this parameter draw.
    """
    f = lambda k: median_os_at_T0(replace(p, k=k), L) - target
    flo, fhi = f(lo), f(hi)
    if flo * fhi > 0:
        return None                      # anchor not attainable -> draw discarded
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


def sweep(p: Params, L, T_grid):
    return [dict(T=float(T), **metrics(T, p, L)) for T in T_grid]
