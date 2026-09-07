"""Phase 2: systemic therapy during the bridge. PC4 + endpoints F-E1..F-E5.
Pre-registered in PREREG.md (Phase 2 section) before this file was executed."""
import json, numpy as np
from dataclasses import replace
from bridge_model import Params, draw_latents, metrics, calibrate_k

N = 600_000   # must match run_main.py exactly for PC4 to be an exact reproduction
T_GRID = np.arange(0, 25, 1.0)
rng = np.random.default_rng(11071905)          # same seed as Phase 1 run_main.py
base = Params()
L = draw_latents(base, N, rng)
k = calibrate_k(base, L)
P = replace(base, k=k)
print(f"calibrated k = {k:.3f}   (therapy has no effect at T=0, so the PC1 anchor is unchanged)\n")

def curve(p):
    m = [metrics(T, p, L) for T in T_GRID]
    rmst = np.array([x["rmst"] for x in m]); a60 = np.array([x["alive60"] for x in m])
    return rmst, a60, np.array([x["futile_tx"] for x in m]), m

# ---------------- PC4 (revised): therapy-off must reproduce the Phase 1 FORMULA exactly
# Original PC4 compared against a stored Phase 1 curve, but that curve was computed under a
# tighter calibrate_k tolerance, so it could never match bitwise. This is the stronger test:
# a literal reference implementation of the Phase 1 survival function, compared element-wise
# at identical k. Deviation recorded in FINDINGS.md.
def survival_phase1_reference(T, p, L):
    """Verbatim Phase 1 logic, before any therapy machinery existed."""
    U, U_det, D = L["U"], L["U_det"], L["D"]
    transplanted = (D > T) & (U_det > T)
    dev_first = (~transplanted) & (D <= U_det)
    det_first = (~transplanted) & (D > U_det)
    surv = np.empty_like(D)
    surv[dev_first] = D[dev_first]
    surv[det_first] = U_det[det_first] + np.minimum(L["met"][det_first], L["dev_resid"][det_first])
    tx = transplanted
    occ_tx = tx & L["occult"] & (U > T)
    free_tx = tx & ~occ_tx
    surv[free_tx] = T + L["graft"][free_tx]
    surv[occ_tx] = T + np.minimum((U[occ_tx] - T) / p.k + L["met_after_tx"][occ_tx] / p.k,
                                  L["graft"][occ_tx])
    surv[tx & L["periop_death"]] = T
    return surv

from bridge_model import survival as survival_phase2
worst = 0.0
for T in T_GRID:
    a = survival_phase2(float(T), P, L)      # therapy off: stasis=1, erad rate=0
    b = survival_phase1_reference(float(T), P, L)
    worst = max(worst, float(np.max(np.abs(a - b))))
pc4 = worst == 0.0
print(f"PC4  therapy-off vs Phase 1 reference formula, all {len(T_GRID)} bridge durations,"
      f" {N:,} patients each")
print(f"     max element-wise |difference| = {worst:.1e}   "
      f"{'PASS (bitwise identical)' if pc4 else 'FAIL'}")
if not pc4:
    raise SystemExit("PC4 failed - Phase 2 void")

# ---------------- PC5: the eradication machinery must actually fire ----------------
# Added after a bug: tx_erad_rate was baked into the latent draw, so varying it via
# replace() silently did nothing and every eradication rate returned identical results.
# This control tests the mechanism directly instead of trusting the endpoint.
from bridge_model import carried_disease
pc5_row = [(r, carried_disease(12.0, replace(P, tx_erad_rate=r), L))
           for r in [0.0, 0.05, 0.3, 1.0]]
pc5 = (pc5_row[0][1] > 0.05 and pc5_row[-1][1] < 0.001
       and all(pc5_row[i][1] > pc5_row[i+1][1] for i in range(len(pc5_row)-1)))
print("PC5  eradication fires: carried disease at T=12 by rate  " +
      "  ".join(f"{r}/mo:{v:.4f}" for r, v in pc5_row) +
      f"   {'PASS' if pc5 else 'FAIL'}")
if not pc5:
    raise SystemExit("PC5 failed - eradication machinery inert, Phase 2 void")

r0, a0, f0, _ = curve(P)

# ---------------- F-E4 / F1: pure cytostasis, no eradication ----------------
print("\nF1  PURE CYTOSTASIS (slows occult disease, does not clear it), eradication = 0")
print(f"{'stasis s':>9} {'T*(RMST)':>9} {'T*(alive60)':>12} {'P(alive60)@T=6':>15} "
      f"{'carried disease@T=6':>21}")
cyto = {}
for s in [1.0, 1.25, 1.5, 2.0, 3.0]:
    p = replace(P, tx_stasis=s)
    r, a, f, _ = curve(p)
    cyto[s] = dict(T_rmst=float(T_GRID[int(np.argmax(r))]),
                   T_a60=float(T_GRID[int(np.argmax(a))]),
                   a60_at6=float(a[6]), carried6=float(f[6]), a60_at0=float(a[0]))
    print(f"{s:>9.2f} {cyto[s]['T_rmst']:>9.0f} {cyto[s]['T_a60']:>12.0f} "
          f"{a[6]:>15.4f} {f[6]:>21.4f}")
d = cyto[3.0]["a60_at6"] - cyto[1.0]["a60_at6"]
print(f"    effect of s=1 -> s=3 on P(alive60) at T=6: {d:+.4f}   "
      f"-> F1 direction {'HARMFUL (as predicted)' if d < 0 else 'NOT harmful (F1 FALSIFIED)'}")

# ---------------- F-E2: the eradication efficacy bar ----------------
print("\nF2  ERADICATION during the bridge: what is the efficacy bar?")
print(f"{'rate/mo':>8} {'P(clear by 6mo)':>16} {'T*(alive60)':>12} {'P(alive60)@T*':>14} "
      f"{'P(alive60)@T=0':>15}")
bar = None; erad = {}
for rate in [0.0, 0.01, 0.02, 0.03, 0.05, 0.075, 0.10, 0.15, 0.20, 0.30]:
    p = replace(P, tx_erad_rate=rate)
    r, a, f, _ = curve(p)
    i = int(np.argmax(a)); Ta = float(T_GRID[i])
    erad[rate] = dict(T_a60=Ta, a60_max=float(a[i]), a60_0=float(a[0]),
                      T_rmst=float(T_GRID[int(np.argmax(r))]))
    if bar is None and Ta > 0: bar = rate
    print(f"{rate:>8.3f} {1-np.exp(-rate*6):>16.3f} {Ta:>12.0f} {a[i]:>14.4f} {a[0]:>15.4f}")
print(f"    F2 efficacy bar: eradication rate ~{bar}/month "
      f"(= {1-np.exp(-bar*6):.0%} chance of clearing occult disease over 6 months)"
      if bar else "    F2: no tested rate moved T* off zero")

# ---------------- F-E3: bar vs device hazard ----------------
print("\nF3  efficacy bar as a function of device hazard")
bars = {}
for dev in [0.010, 0.020, 0.025, 0.035, 0.045]:
    pb = replace(P, dev_rate=dev)
    Lb = draw_latents(pb, N, np.random.default_rng(4242))
    kb = calibrate_k(pb, Lb); pb = replace(pb, k=kb)
    b_ = None
    for rate in [0.0,0.005,0.01,0.02,0.03,0.05,0.075,0.10,0.15,0.20,0.30,0.50]:
        a = np.array([metrics(T, replace(pb, tx_erad_rate=rate), Lb)["alive60"] for T in T_GRID])
        if np.argmax(a) > 0: b_ = rate; break
    bars[dev] = b_
    txt = (f"{b_:.3f}/mo ({1-np.exp(-b_*6):.0%} clearance by 6mo)") if b_ is not None else "unreachable"
    print(f"    device hazard {dev:.3f}/mo -> bar {txt}")

json.dump(dict(pc4_pass=bool(pc4), k=k, cytostasis=cyto, eradication=erad,
               efficacy_bar=bar, bar_vs_device={str(x): y for x, y in bars.items()}),
          open("results/therapy.json","w"), indent=2)
