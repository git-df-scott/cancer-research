"""Endpoints E1-E3: base-case bridge-duration curve, with k calibrated to the
external PC1 anchor. Reports every pre-registered endpoint regardless of outcome."""
import json, numpy as np
from dataclasses import replace, asdict
from bridge_model import Params, draw_latents, sweep, calibrate_k

N = 600_000
rng = np.random.default_rng(11071905)
T_GRID = np.arange(0, 25, 1.0)

base = Params()
L = draw_latents(base, N, rng)
k = calibrate_k(base, L)
p = replace(base, k=k)
print(f"calibrated k = {k:.3f}\n")

rows = sweep(p, L, T_GRID)
rmst = np.array([r["rmst"] for r in rows])
a60  = np.array([r["alive60"] for r in rows])
i    = int(np.argmax(rmst))
T_star = float(T_GRID[i])

print(f"{'T':>4} {'RMST-60':>9} {'P(alive60)':>11} {'transplanted':>13} {'futile tx':>10}")
for r in rows:
    mark = "  <-- T*" if r["T"] == T_star else ""
    print(f"{r['T']:>4.0f} {r['rmst']:>9.2f} {r['alive60']:>11.3f} "
          f"{r['tx_rate']:>13.3f} {r['futile_tx']:>10.3f}{mark}")

# E3: plateau width — every T within 1.0 month of the maximum
plateau = [float(t) for t, v in zip(T_GRID, rmst) if v >= rmst[i] - 1.0]
gap = float(rmst.max() - rmst.min())

# pre-registered decision rule R2
interior = 0.0 < T_star < float(T_GRID[-1])
negligible = gap < 1.0

print(f"\nE3  T* = {T_star:.0f} months   RMST-60 = {rmst[i]:.2f} mo")
print(f"    plateau within 1 mo of max: T = {plateau[0]:.0f}..{plateau[-1]:.0f}")
print(f"    RMST-60 at T=0  : {rmst[0]:.2f} mo   (transplant immediately)")
print(f"    RMST-60 at T=24 : {rmst[-1]:.2f} mo   (bridge two years)")
print(f"    best-worst gap  : {gap:.2f} mo")
print(f"    P(alive at 60mo): {a60[0]:.3f} at T=0  ->  {a60[i]:.3f} at T*   "
      f"(absolute gain {a60[i]-a60[0]:+.3f})")
print(f"\n    interior optimum: {interior}    R2 clinically negligible: {negligible}")

json.dump(dict(params=asdict(p), T_grid=T_GRID.tolist(), rows=rows,
               T_star=T_star, plateau=[plateau[0], plateau[-1]], rmst_gap=gap,
               interior=interior, negligible_R2=negligible),
          open("results/base_case.json", "w"), indent=2)
