"""Translational output: the minimum prospective dataset that would settle this.

C1 established that p_R and p_dur are not separately identifiable and only their product
pi = durable clearance fraction drives the strategy. So the registry must measure DURABLE
clearance, not response rate. This computes how many patients that takes, and the other
axis of the boundary: how good a device programme would have to be for the cardiac-
evidence value of pi to reach the MCID.
"""
import json, numpy as np
from scipy import stats
from dataclasses import replace
from phase3_model import P3, draw_latents3, metrics3, calibrate_k3

MCID = 0.05
PI_EVIDENCE = 0.07      # cardiac-specific: 1/7 visceral responders x 2/4 durable
PI_BOUNDARY = 0.20      # required at best device + max cytostasis (section F)

# ---------- 1. exact binomial sample size to distinguish the two ----------
print("REGISTRY SIZE: distinguishing pi = 0.07 (cardiac evidence) from pi = 0.20 (boundary)")
print("  single-arm, one-sided alpha = 0.05, exact binomial\n")
def power_exact(n, p0, p1, alpha=0.05):
    # smallest critical count c with P(X >= c | p0) <= alpha
    c = next((c for c in range(n + 1) if stats.binom.sf(c - 1, n, p0) <= alpha), n + 1)
    return (stats.binom.sf(c - 1, n, p1) if c <= n else 0.0), c
for n in range(10, 121, 5):
    pw, c = power_exact(n, PI_EVIDENCE, PI_BOUNDARY)
    if pw >= 0.80:
        print(f"  n = {n:>3}  reject H0 if >= {c} of {n} achieve durable clearance   power = {pw:.3f}")
        n80 = n
        break
for n in [20, 30, 40, 50, 60, 80]:
    pw, c = power_exact(n, PI_EVIDENCE, PI_BOUNDARY)
    print(f"    n={n:>3}: power {pw:.2f} (reject if >= {c} events)")

# ---------- 2. the other axis: how good must the device programme be? ----------
print("\n\nDEVICE AXIS: what device hazard would let the cardiac-evidence pi reach the MCID?")
N = 200_000
T_GRID = np.arange(0, 25, 1.0)
L = draw_latents3(P3(), N, np.random.default_rng(99))
P = replace(P3(), k=calibrate_k3(P3(), L))
def benefit(p):
    kk = calibrate_k3(p, L); p = replace(p, k=kk) if kk else p
    a = np.array([metrics3(float(T), p, L)["alive60"] for T in T_GRID])
    i = int(np.argmax(a))
    return float(a[i] - a[0]), float(T_GRID[i])

fav = dict(p_dur=1.0, m_clear=3.0, a_clear=1.0, tx_stasis=3.0, med_unmask=9.0)
print(f"  at pi = {PI_EVIDENCE} with every other parameter at its most favourable value:\n")
print(f"  {'device hazard /mo':>18}{'benefit (pp)':>15}{'T*':>5}")
dev_rows = []
for dev in [0.045, 0.025, 0.010, 0.005, 0.002, 0.001, 0.0001]:
    g, T = benefit(replace(P, p_R=PI_EVIDENCE, dev_rate=dev, **fav))
    dev_rows.append(dict(dev=dev, gain=g, T=T))
    print(f"  {dev:>18.4f}{g*100:>15.2f}{T:>5.0f}   {'MEETS MCID' if g>=MCID else ''}")
best = max(r["gain"] for r in dev_rows)
print(f"\n  ceiling as device hazard -> 0 : {best*100:.2f} pp")
print(f"  -> the cardiac-evidence pi {'CAN' if best>=MCID else 'CANNOT'} reach the MCID even")
print(f"     with a physically perfect device and every other parameter at its best value.")

json.dump(dict(n80=n80, pi_evidence=PI_EVIDENCE, pi_boundary=PI_BOUNDARY,
               device_axis=dev_rows, ceiling_at_zero_device=best),
          open("results/phase3_registry.json","w"), indent=2)
