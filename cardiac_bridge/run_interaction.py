"""F-E5: do cytostasis and eradication interact, or are they separable?

This is the practical question. A VEGFR/KDR inhibitor is predominantly cytostatic.
An eradicating modality would be given alongside it, not instead of it. F1 showed
cytostasis corrupts the filter; F2 showed eradication is what makes waiting pay.
So: does the cytostatic component cancel the eradication benefit?
"""
import json, numpy as np
from dataclasses import replace
from bridge_model import Params, draw_latents, metrics, calibrate_k, carried_disease

N = 400_000
T_GRID = np.arange(0, 25, 1.0)
rng = np.random.default_rng(11071905)
base = Params()
L = draw_latents(base, N, rng)
P = replace(base, k=calibrate_k(base, L))

STASIS = [1.0, 1.5, 2.0, 3.0]
ERAD = [0.0, 0.02, 0.05, 0.10, 0.20]

print("P(alive at 60 months) at the best bridge duration, and that duration in brackets")
print(f"{'':>12}" + "".join(f"{'erad ' + str(e):>18}" for e in ERAD))
grid = {}
for s in STASIS:
    cells = []
    for e in ERAD:
        p = replace(P, tx_stasis=s, tx_erad_rate=e)
        a = np.array([metrics(T, p, L)["alive60"] for T in T_GRID])
        i = int(np.argmax(a))
        grid[(s, e)] = dict(a60=float(a[i]), T=float(T_GRID[i]))
        cells.append(f"{a[i]:.4f} (T={T_GRID[i]:.0f})")
    print(f"stasis {s:>4.1f}  " + "".join(f"{c:>18}" for c in cells))

print("\nEffect of adding cytostasis on top of a fixed eradication rate:")
for e in ERAD:
    b = grid[(1.0, e)]["a60"]
    row = "  ".join(f"s={s}:{grid[(s,e)]['a60']-b:+.4f}" for s in STASIS[1:])
    print(f"    erad {e:>5}/mo (vs no stasis, P={b:.4f}):  {row}")

# is the combination better than either alone?
lone_e = grid[(1.0, 0.05)]["a60"]; lone_s = grid[(2.0, 0.0)]["a60"]
both = grid[(2.0, 0.05)]["a60"]; neither = grid[(1.0, 0.0)]["a60"]
print(f"\n  neither: {neither:.4f}   cytostasis only: {lone_s:.4f}   "
      f"eradication only: {lone_e:.4f}   both: {both:.4f}")
print(f"  additive prediction {neither + (lone_s-neither) + (lone_e-neither):.4f} vs "
      f"observed {both:.4f}  -> interaction {both - (lone_s + lone_e - neither):+.4f} "
      f"({'SYNERGY' if both > lone_s + lone_e - neither else 'ANTAGONISM'})")

json.dump({f"s{k[0]}_e{k[1]}": v for k, v in grid.items()},
          open("results/interaction.json", "w"), indent=2)
