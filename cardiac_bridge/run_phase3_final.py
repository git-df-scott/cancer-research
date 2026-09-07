"""Phase 3 definitive run on the bug-fixed model (all parameters evaluation-time,
PC6-PC13 passing). Sections C2, F, G, H plus the admissible-region Monte Carlo
needed for pre-registered decision rules D1-D4."""
import json, numpy as np
from dataclasses import replace, asdict
from phase3_model import P3, draw_latents3, metrics3, calibrate_k3

N = 200_000
T_GRID = np.arange(0, 25, 1.0)
MCID = 0.05
rng = np.random.default_rng(31415)
base = P3()
L = draw_latents3(base, N, rng)          # unit draws; every parameter applies at eval time
P = replace(base, k=calibrate_k3(base, L))
print(f"calibrated k = {P.k:.3f}   MCID = {MCID*100:.0f} pp (pre-registered)\n")

def benefit(p):
    a = np.array([metrics3(float(T), p, L)["alive60"] for T in T_GRID])
    i = int(np.argmax(a))
    return float(a[i] - a[0]), float(T_GRID[i])

# ---------------- C2 (now with dev_rate and med_unmask actually live) ----------------
print("C2  absolute sensitivity: pp change in benefit per 10% of admissible range")
RANGES = dict(p_R=(0.0,0.30), p_dur=(0.20,0.70), m_clear=(2.0,8.0), a_clear=(1.0,2.5),
              tx_stasis=(1.0,3.0), med_unmask=(3.0,12.0), dev_rate=(0.010,0.045))
for label, pt in [("evidence centre", dict(p_R=0.15, p_dur=0.45, m_clear=5.0, a_clear=1.75,
                                           tx_stasis=2.0, med_unmask=7.5, dev_rate=0.0275)),
                  ("generous corner", dict(p_R=0.30, p_dur=0.70, m_clear=3.0, a_clear=1.0,
                                           tx_stasis=3.0, med_unmask=9.0, dev_rate=0.010))]:
    c = replace(P, **pt)
    kk = calibrate_k3(c, L); c = replace(c, k=kk) if kk else c
    g0, T0 = benefit(c)
    print(f"\n  {label}: benefit {g0*100:.2f} pp at T*={T0:.0f}  (pi = {pt['p_R']*pt['p_dur']:.3f})")
    sens = {}
    for nme, (lo, hi) in RANGES.items():
        st = 0.10 * (hi - lo); v = pt[nme]
        sens[nme] = (benefit(replace(c, **{nme: min(v+st, hi)}))[0]
                     - benefit(replace(c, **{nme: max(v-st, lo)}))[0]) / 2 * 100
    for nme, v in sorted(sens.items(), key=lambda x: -abs(x[1])):
        print(f"      {nme:<12} {v:+.3f} pp")

# ---------------- F: phase boundary ----------------
PI = np.round(np.arange(0.0, 0.92, 0.02), 3)
def required_pi(**kw):
    p = replace(P, p_dur=1.0, **kw)
    kk = calibrate_k3(p, L)
    if kk is None: return None, None
    p = replace(p, k=kk)
    for pi in PI:
        g, T = benefit(replace(p, p_R=float(pi)))
        if g >= MCID: return float(pi), T
    return None, None

print("\n\nF  PHASE BOUNDARY: minimum pi = p_R x p_dur reaching the 5 pp MCID")
print("   admissible ceiling is pi = 0.30 x 0.70 = 0.21\n")
STAS, DEVS = [1.0, 2.0, 3.0], [0.010, 0.025, 0.045]
print(f"{'device hazard':>14}" + "".join(f"{'stasis '+str(s):>18}" for s in STAS))
diagram = {}
for dev in DEVS:
    cells = []
    for s in STAS:
        pi, T = required_pi(dev_rate=dev, tx_stasis=s, m_clear=4.0, med_unmask=6.0)
        diagram[f"dev{dev}_s{s}"] = dict(required_pi=pi, T_star=T)
        cells.append(f"pi>={pi:.2f} (T={T:.0f})" if pi is not None else "unreachable")
    print(f"{dev:>14.3f}" + "".join(f"{c:>18}" for c in cells))

print("\n   best case for the hypothesis (dev=0.010, s=3), varying the clearance/unmasking race:")
for mc, mu in [(2.0, 12.0), (3.0, 9.0), (4.0, 6.0), (6.0, 3.0)]:
    pi, T = required_pi(dev_rate=0.010, tx_stasis=3.0, m_clear=mc, med_unmask=mu)
    diagram[f"race_mc{mc}_mu{mu}"] = dict(required_pi=pi, T_star=T)
    print(f"      clearance {mc:.0f} mo vs unmasking {mu:>4.0f} mo -> "
          + (f"pi >= {pi:.2f} (T*={T:.0f})" if pi is not None else "unreachable"))

# ---------------- admissible-region Monte Carlo (decision rules D1-D4) --------------
print("\n\nADMISSIBLE-REGION MONTE CARLO (robustness threshold: MCID in >=50% of draws)")
draws, r2 = [], np.random.default_rng(2718)
for _ in range(300):
    d = dict(p_R=r2.uniform(0.0, 0.30), p_dur=r2.uniform(0.20, 0.70),
             m_clear=r2.uniform(2.0, 8.0), a_clear=r2.uniform(1.0, 2.5),
             tx_stasis=r2.uniform(1.0, 3.0), med_unmask=r2.uniform(3.0, 12.0),
             dev_rate=r2.uniform(0.010, 0.045))
    p = replace(P, **d)
    kk = calibrate_k3(p, L)
    if kk is None: continue
    g, T = benefit(replace(p, k=kk))
    draws.append(dict(**d, pi=d["p_R"]*d["p_dur"], gain=g, T=T))
g = np.array([x["gain"] for x in draws]); Ts = np.array([x["T"] for x in draws])
pis = np.array([x["pi"] for x in draws])
print(f"    usable draws                      : {len(draws)}")
print(f"    reaching the 5 pp MCID            : {(g>=MCID).mean()*100:.1f}%   "
      f"(threshold 50%) -> {'MEETS' if (g>=MCID).mean()>=0.5 else 'FAILS'}")
print(f"    with any benefit at all (T* > 0)  : {(Ts>0).mean()*100:.1f}%")
print(f"    median benefit                    : {np.median(g)*100:.2f} pp")
print(f"    90th percentile benefit           : {np.percentile(g,90)*100:.2f} pp")
print(f"    max benefit anywhere in the region: {g.max()*100:.2f} pp (at pi={pis[g.argmax()]:.3f})")
print(f"    admissible pi: median {np.median(pis):.3f}, max {pis.max():.3f}")

json.dump(dict(mcid=MCID, diagram=diagram, draws=draws,
               frac_mcid=float((g>=MCID).mean()), frac_positive_T=float((Ts>0).mean()),
               median_gain=float(np.median(g)), max_gain=float(g.max())),
          open("results/phase3_final.json","w"), indent=2)
