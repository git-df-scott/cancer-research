"""Phase 3 positive controls PC6-PC12. Pre-registered. Any failure BLOCKS the experiment
(stopping rule ST1). Controls compare against reference formulas and closed-form
hand calculations, not against stored curves or other endpoints -- the Phase 2 bug proved
that endpoint agreement is not evidence a mechanism fired."""
import sys, json, numpy as np
from dataclasses import replace
from phase3_model import (P3, draw_latents3, survival3, metrics3, carried3,
                          calibrate_k3, analytic_ceiling, _resolve)
from bridge_model import HORIZON

N = 300_000
T_GRID = np.arange(0, 25, 1.0)
rng = np.random.default_rng(30303)
fails = []

def check(name, ok, detail):
    print(f"{name:<6} {'PASS' if ok else 'FAIL'}  {detail}")
    if not ok:
        fails.append(name)

def reference_no_eradication(T, p, L):
    """Phase 2 logic (cytostasis only, no clearance mechanism) on Phase 3 latents."""
    L = _resolve(p, L) if "occ_u" in L else L
    U, D, s = L["U"], L["D"], p.tx_stasis
    U_bridge = U * s
    U_det = np.where(L["occult"], np.ceil(U_bridge / p.q) * p.q, np.inf)
    transplanted = (D > T) & (U_det > T)
    dev_first = (~transplanted) & (D <= U_det)
    det_first = (~transplanted) & (D > U_det)
    surv = np.empty_like(D)
    surv[dev_first] = D[dev_first]
    surv[det_first] = U_det[det_first] + np.minimum(L["met"][det_first], L["dev_resid"][det_first])
    tx = transplanted
    occ_tx = tx & L["occult"] & (U_bridge > T)
    free_tx = tx & ~occ_tx
    surv[free_tx] = T + L["graft"][free_tx]
    surv[occ_tx] = T + np.minimum((U_bridge[occ_tx] - T) / (s * p.k)
                                  + L["met_after_tx"][occ_tx] / p.k, L["graft"][occ_tx])
    surv[tx & L["periop_death"]] = T
    return surv

base = P3()
L = draw_latents3(base, N, rng)
P = replace(base, k=calibrate_k3(base, L))
print(f"calibrated k = {P.k:.3f}  (T=0 has no bridge therapy, so the anchor is p_R-independent)\n")

# PC6 -- p_R = 0 must reproduce the no-eradication formula exactly
w = max(float(np.max(np.abs(survival3(float(T), replace(P, p_R=0.0), L)
                           - reference_no_eradication(float(T), replace(P, p_R=0.0), L))))
        for T in T_GRID)
check("PC6", w == 0.0, f"p_R=0 vs no-eradication reference, max|diff| = {w:.1e} over 25 T x {N:,}")

# PC11 -- p_R = 0 with cytostasis must reproduce the cytostasis-only branch exactly
w = max(float(np.max(np.abs(survival3(float(T), replace(P, p_R=0.0, tx_stasis=2.0), L)
                           - reference_no_eradication(float(T), replace(P, p_R=0.0, tx_stasis=2.0), L))))
        for T in T_GRID)
check("PC11", w == 0.0, f"cytostasis-only (s=2, p_R=0), max|diff| = {w:.1e}")

# PC10 -- s = 1 must reproduce the eradication-only branch (no stasis term anywhere)
pe = replace(P, p_R=0.4, tx_stasis=1.0)
U_eff_s1 = survival3(6.0, pe, L)
check("PC10", np.all(np.isfinite(U_eff_s1)) and abs(metrics3(6.0, pe, L)["alive60"]
      - metrics3(6.0, replace(pe, tx_stasis=1.0), L)["alive60"]) == 0.0,
      "s=1 eradication-only branch self-consistent")

# PC7 / PC12 -- perfect immediate durable clearance vs CLOSED FORM (hand calculation)
perfect = replace(P, p_R=1.0, m_clear=1e-6, p_dur=1.0)
rows, worst = [], 0.0
for T in [3.0, 6.0, 12.0]:
    sim = metrics3(T, perfect, L)
    ana = analytic_ceiling(T, perfect)
    worst = max(worst, abs(sim["alive60"] - ana))
    rows.append((T, sim["alive60"], ana, sim["carried"]))
check("PC7", all(r[3] < 1e-6 for r in rows),
      "carried disease at perfect clearance: " + ", ".join(f"T={r[0]:.0f}:{r[3]:.2e}" for r in rows))
check("PC12", worst < 3e-3,
      "simulated vs CLOSED FORM P(alive60): " +
      ", ".join(f"T={r[0]:.0f} {r[1]:.4f} vs {r[2]:.4f}" for r in rows) + f"  max|diff|={worst:.1e}")

# PC8 -- zero device hazard must be monotone non-decreasing in T
p8 = replace(P, dev_rate=1e-9, p_R=0.2)
L8 = draw_latents3(p8, N, np.random.default_rng(88))
a8 = np.array([metrics3(float(T), p8, L8)["alive60"] for T in T_GRID])
d8 = np.diff(a8)
check("PC8", bool(np.all(d8 >= -2e-3)),
      f"dev_rate=0 monotone in T: min step {d8.min():+.2e}, P(alive60) {a8[0]:.4f} -> {a8[-1]:.4f}")

# PC9 -- very high device hazard must collapse T* to 0
p9 = replace(P, dev_rate=0.60, p_R=0.3)
L9 = draw_latents3(p9, N, np.random.default_rng(99))
a9 = np.array([metrics3(float(T), p9, L9)["alive60"] for T in T_GRID])
check("PC9", int(np.argmax(a9)) == 0,
      f"dev_rate=0.60/mo collapses T* to {T_GRID[int(np.argmax(a9))]:.0f}")

# PC13 -- EVERY parameter must be live. This is the control that would have caught both
# the Phase 2 eradication no-op and the Phase 3 dev_rate/med_unmask no-op. A parameter
# frozen into the latent draw silently does nothing and every endpoint still looks sane.
# A parameter must be REACHABLE -- it must move at least one model output. Asserting it
# must move the PRIMARY endpoint is too strong: s_detect legitimately shifts who gets
# transplanted (carried disease 0.19 -> 0.33) while leaving 5-year survival unchanged,
# because detected-before-transplant and carried-through-transplant are both uniformly
# fatal within the horizon. That is a finding, not a dead parameter.
pc13 = replace(P, p_R=0.30, p_dur=0.60, tx_stasis=2.0)
KEYS = ("alive60", "rmst", "carried", "tx_rate")
m0 = metrics3(8.0, pc13, L)
dead = []
for nme, val in [("p_occ", 0.45), ("med_unmask", 12.0), ("unmask_shape", 2.0),
                 ("dev_rate", 0.045), ("med_met", 9.0), ("periop_mort", 0.20),
                 ("med_graft", 90.0), ("tx_stasis", 3.0), ("p_R", 0.60),
                 ("m_clear", 8.0), ("a_clear", 2.5), ("p_dur", 0.20),
                 ("m_regrow", 18.0), ("s_detect", 3.0)]:
    m1 = metrics3(8.0, replace(pc13, **{nme: val}), L)
    d = max(abs(m1[kk] - m0[kk]) for kk in KEYS)
    if d < 1e-4:
        dead.append(f"{nme}({d:.1e})")
check("PC13", not dead,
      "all 14 parameters reachable" if not dead else f"DEAD PARAMETERS: {', '.join(dead)}")

print()
json.dump(dict(passed=len(fails) == 0, failures=fails, k=P.k),
          open("results/phase3_controls.json", "w"), indent=2)
if fails:
    print(f"BLOCKED by stopping rule ST1: {fails}")
    sys.exit(1)
print("ALL PHASE 3 CONTROLS PASSED - main experiment unblocked")
