"""Phase 3 main experiment: sections C2 (corrected), F (phase diagram), G (synergy
attack), H (immunotherapy window). Pre-registered; MCID = 5.0 pp fixed in advance."""
import json, numpy as np
from dataclasses import replace
from phase3_model import P3, draw_latents3, metrics3, calibrate_k3

N = 250_000
T_GRID = np.arange(0, 25, 1.0)
MCID = 0.05                      # 5 percentage points absolute, pre-registered
rng = np.random.default_rng(717)
base = P3()
L = draw_latents3(base, N, rng)
P = replace(base, k=calibrate_k3(base, L))
print(f"calibrated k = {P.k:.3f}   MCID = {MCID*100:.0f} pp absolute (pre-registered)\n")

def benefit(p, Lx=None):
    Lx = L if Lx is None else Lx
    a = np.array([metrics3(float(T), p, Lx)["alive60"] for T in T_GRID])
    i = int(np.argmax(a))
    return float(a[i] - a[0]), float(T_GRID[i])

# ---------------- C2 corrected: absolute sensitivity, elasticity is undefined at g=0 ----
print("C2 (corrected)  absolute sensitivity: pp change in benefit per 10% of admissible range")
print("   The original C2 divided by the benefit, which is ZERO at the centre of the")
print("   admissible region, producing meaningless elasticities. Recorded as a correction.\n")
RANGES = dict(p_R=(0.0,0.30), p_dur=(0.20,0.70), m_clear=(2.0,8.0), a_clear=(1.0,2.5),
              tx_stasis=(1.0,3.0), med_unmask=(3.0,12.0), dev_rate=(0.010,0.045))
for label, pt in [("centre of admissible region",
                   dict(p_R=0.15, p_dur=0.45, m_clear=5.0, a_clear=1.75, tx_stasis=2.0,
                        med_unmask=7.5, dev_rate=0.0275)),
                  ("generous corner (p_R=0.30, p_dur=0.70)",
                   dict(p_R=0.30, p_dur=0.70, m_clear=3.0, a_clear=1.0, tx_stasis=3.0,
                        med_unmask=9.0, dev_rate=0.012))]:
    c = replace(P, **pt); g0, T0 = benefit(c)
    print(f"  {label}: benefit {g0*100:.2f} pp at T*={T0:.0f}")
    sens = {}
    for nme, (lo, hi) in RANGES.items():
        step = 0.10 * (hi - lo)
        v = pt[nme]
        gp, _ = benefit(replace(c, **{nme: min(v + step, hi)}))
        gm, _ = benefit(replace(c, **{nme: max(v - step, lo)}))
        sens[nme] = (gp - gm) / 2 * 100
    for nme, v in sorted(sens.items(), key=lambda x: -abs(x[1])):
        print(f"      {nme:<12} {v:+.3f} pp")
    print()

# ---------------- F: the phase diagram -- required pi = p_R x p_dur to reach MCID -------
print("F  PHASE BOUNDARY: minimum effective durable clearance fraction pi = p_R x p_dur")
print("   required to reach the 5 pp MCID. '-' = MCID unreachable at any pi <= 0.90.\n")
PI = np.round(np.arange(0.0, 0.92, 0.02), 3)
def required_pi(s, dev, m_clear=4.0, med_unmask=6.0, a_clear=1.0):
    Lx = draw_latents3(replace(P, dev_rate=dev, med_unmask=med_unmask), N,
                       np.random.default_rng(4321))
    pb = replace(P, dev_rate=dev, med_unmask=med_unmask, tx_stasis=s,
                 m_clear=m_clear, a_clear=a_clear, p_dur=1.0)
    kk = calibrate_k3(pb, Lx)
    if kk is None: return None, None
    pb = replace(pb, k=kk)
    for pi in PI:
        g, T = benefit(replace(pb, p_R=float(pi)), Lx)
        if g >= MCID:
            return float(pi), T
    return None, None

DEVS = [0.010, 0.025, 0.045]
STAS = [1.0, 2.0, 3.0]
print(f"{'device hazard':>14}" + "".join(f"{'stasis '+str(s):>16}" for s in STAS))
diagram = {}
for dev in DEVS:
    cells = []
    for s in STAS:
        pi, T = required_pi(s, dev)
        diagram[f"dev{dev}_s{s}"] = dict(required_pi=pi, T_star=T)
        cells.append(f"pi>={pi:.2f} (T={T:.0f})" if pi is not None else "unreachable")
    print(f"{dev:>14.3f}" + "".join(f"{c:>16}" for c in cells))

print("\n   effect of faster clearance and slower unmasking (both favourable), dev=0.010, s=3:")
for mc, mu in [(2.0, 12.0), (3.0, 9.0), (4.0, 6.0), (6.0, 3.0)]:
    pi, T = required_pi(3.0, 0.010, m_clear=mc, med_unmask=mu)
    print(f"      m_clear={mc:.0f} mo, med_unmask={mu:>4.0f} mo -> "
          + (f"pi >= {pi:.2f} (T*={T:.0f})" if pi is not None else "unreachable"))

json.dump(dict(mcid=MCID, diagram=diagram), open("results/phase3_diagram.json","w"), indent=2)
