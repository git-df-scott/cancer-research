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
    return dict(occult=occult, U=U, U_det=U_det, D=D, met=met,
                met_after_tx=met_after_tx, dev_resid=dev_resid,
                graft=graft, periop_death=periop_death)


def survival(T, p: Params, L):
    """Overall survival from cardiectomy under bridge duration T."""
    U, U_det, D = L["U"], L["U_det"], L["D"]

    transplanted = (D > T) & (U_det > T)          # neither event happened before transplant
    dev_first    = (~transplanted) & (D <= U_det)  # branch A: died on device
    det_first    = (~transplanted) & (D > U_det)   # branch B: mets found, transplant cancelled

    surv = np.empty_like(D)

    # A. device event before anything else
    surv[dev_first] = D[dev_first]

    # B. metastases detected during the bridge -> no transplant, stays on device
    surv[det_first] = U_det[det_first] + np.minimum(L["met"][det_first],
                                                    L["dev_resid"][det_first])

    # C/D. transplanted at T
    tx = transplanted
    occ_tx = tx & L["occult"] & (U > T)   # D: occult disease carried through transplant
    free_tx = tx & ~occ_tx                # C: truly disease-free (or occult already excluded)

    # C: disease-free graft survival
    surv[free_tx] = T + L["graft"][free_tx]

    # D: immunosuppression accelerates residual unmasking AND post-detection survival by k
    resid = (U[occ_tx] - T) / p.k
    post = L["met_after_tx"][occ_tx] / p.k
    surv[occ_tx] = T + np.minimum(resid + post, L["graft"][occ_tx])

    # perioperative mortality applies to everyone actually transplanted
    surv[tx & L["periop_death"]] = T
    return surv


def metrics(T, p: Params, L):
    s = survival(T, p, L)
    return dict(
        rmst=float(np.mean(np.minimum(s, HORIZON))),     # E1
        alive60=float(np.mean(s >= HORIZON)),            # E2
        median=float(np.median(s)),
        tx_rate=float(np.mean((L["D"] > T) & (L["U_det"] > T))),
        futile_tx=float(np.mean((L["D"] > T) & (L["U_det"] > T) & L["occult"] & (L["U"] > T))),
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
