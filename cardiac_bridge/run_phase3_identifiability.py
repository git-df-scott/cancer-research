"""Section C: is the Phase 3 model identifiable from available external evidence?

Run BEFORE the bridge experiment. If parameters are not constrainable, stopping rule ST2
says report classification B rather than picking a favourite point.
"""
import json, itertools, numpy as np
from dataclasses import replace
from phase3_model import P3, draw_latents3, metrics3, calibrate_k3

N = 200_000
T_GRID = np.arange(0, 25, 1.0)
rng = np.random.default_rng(5150)
base = P3()
L = draw_latents3(base, N, rng)
P = replace(base, k=calibrate_k3(base, L))

def benefit(p):
    """Primary endpoint gain: max_T P(alive60) - P(alive60) at T=0, and the arg."""
    a = np.array([metrics3(float(T), p, L)["alive60"] for T in T_GRID])
    i = int(np.argmax(a))
    return float(a[i] - a[0]), float(T_GRID[i]), float(a[i])

# ---------- C1: are p_R and p_dur separately identifiable, or only their product? ----------
print("C1  p_R x p_dur degeneracy: benefit (pp gain in 5-yr survival) on a grid\n")
pRs = [0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.60]
pDs = [0.10, 0.20, 0.30, 0.50, 0.70, 1.00]
cells = {}
print("  pR\\pdur" + "".join(f"{d:>9.2f}" for d in pDs))
for r in pRs:
    row = []
    for d in pDs:
        g, T, _ = benefit(replace(P, p_R=r, p_dur=d))
        cells[(r, d)] = dict(gain=g, T=T, product=r * d)
        row.append(g)
    print(f"{r:>10.2f}" + "".join(f"{g*100:>9.2f}" for g in row))

# do level sets follow constant product?
prods = np.array([v["product"] for v in cells.values()])
gains = np.array([v["gain"] for v in cells.values()])
order = np.argsort(prods)
# bin by product and measure within-bin spread of the gain
bins = np.linspace(prods.min(), prods.max(), 9)
idx = np.digitize(prods, bins)
spreads, ranges = [], []
for b in np.unique(idx):
    g = gains[idx == b]
    if len(g) >= 3:
        spreads.append(float(np.std(g)))
        ranges.append(float(g.max() - g.min()))
tot = float(np.std(gains))
print(f"\n    total SD of gain across the grid       : {tot*100:.2f} pp")
print(f"    mean SD WITHIN constant-product bins   : {np.mean(spreads)*100:.2f} pp")
print(f"    variance explained by the product alone: {1 - np.mean(np.array(spreads)**2)/tot**2:.1%}")
degenerate = np.mean(spreads) / tot < 0.35
print(f"    -> p_R and p_dur are {'NOT separately identifiable (product only)' if degenerate else 'separately identifiable'}")

# ---------- C2: normalised sensitivity (Jacobian) at a central admissible point ----------
print("\nC2  normalised sensitivity of the benefit to each parameter (central point)")
centre = replace(P, p_R=0.15, p_dur=0.45, m_clear=4.0, a_clear=1.0, tx_stasis=2.0,
                 med_unmask=6.0, dev_rate=0.025)
g0, T0, _ = benefit(centre)
jac = {}
for name, val, step in [("p_R", 0.15, 0.03), ("p_dur", 0.45, 0.09), ("m_clear", 4.0, 0.8),
                        ("a_clear", 1.0, 0.25), ("tx_stasis", 2.0, 0.4),
                        ("med_unmask", 6.0, 1.2), ("dev_rate", 0.025, 0.005)]:
    gp, _, _ = benefit(replace(centre, **{name: val + step}))
    gm, _, _ = benefit(replace(centre, **{name: max(val - step, 1e-6)}))
    # elasticity: relative change in benefit per relative change in parameter
    jac[name] = float((gp - gm) / (2 * step) * val / (abs(g0) + 1e-9))
for nme, v in sorted(jac.items(), key=lambda x: -abs(x[1])):
    print(f"    {nme:<12} elasticity {v:+.2f}")

# ---------- C3: what does external evidence actually constrain? ----------
print("\nC3  external constraint on each parameter (from PROVENANCE.md)")
constraint = {
    "p_R":        ("0.00-0.30", "0/1 cardiac, 1/7 visceral in DART; ~0 adjuvant chemo", "WEAK"),
    "p_dur":      ("0.20-0.70", "2/4 DART responses ongoing at ~12 mo", "VERY WEAK"),
    "m_clear":    ("2-8 mo",    "time-to-response in solid tumours, indirect", "WEAK"),
    "a_clear":    ("1.0-2.5",   "no angiosarcoma data at all", "NONE"),
    "tx_stasis":  ("1.0-3.0",   "PALETTE PFS ratio 4.6/1.6; JCOG AS weaker", "MODERATE"),
    "med_unmask": ("3-12 mo",   "IGR cardiac AS median TTP 3 mo", "MODERATE"),
    "dev_rate":   ("0.010-0.045", "TAH programme data, centre-specific", "MODERATE"),
}
for nme, (rng_, basis, conf) in constraint.items():
    print(f"    {nme:<12} {rng_:<14} {conf:<10} {basis}")

json.dump(dict(grid={f"{r}_{d}": v for (r, d), v in cells.items()},
               product_degenerate=bool(degenerate),
               within_bin_sd=float(np.mean(spreads)), total_sd=tot,
               jacobian=jac, centre_gain=g0, centre_T=T0),
          open("results/phase3_identifiability.json", "w"), indent=2)
