"""Positive controls PC1-PC3, pre-registered. These must pass before any result
about optimal bridge duration is reportable."""
import json, numpy as np
from dataclasses import replace
from bridge_model import Params, draw_latents, metrics, median_os_at_T0, calibrate_k, CAL_TARGET

N = 400_000
rng = np.random.default_rng(20260907)
T_GRID = np.arange(0, 25, 1.0)
out = {}

base = Params()
L = draw_latents(base, N, rng)

# ---------------- PC1: external calibration anchor ----------------
pre = median_os_at_T0(base, L)
k_cal = calibrate_k(base, L)
base_cal = replace(base, k=k_cal)
post = median_os_at_T0(base_cal, L)
pc1_pass = abs(post - CAL_TARGET) < 0.25
out["PC1"] = dict(target=CAL_TARGET, median_before_calibration=pre,
                  k_uncalibrated=base.k, k_calibrated=k_cal,
                  median_after_calibration=post, passed=bool(pc1_pass))
print(f"PC1 median OS at T=0, target {CAL_TARGET} mo")
print(f"    before calibration (k={base.k}):  {pre:.2f} mo")
print(f"    calibrated k = {k_cal:.3f}  ->  {post:.2f} mo   {'PASS' if pc1_pass else 'FAIL'}")

# ---------------- PC2: no occult disease -> T* must be 0 ----------------
p2 = replace(base_cal, p_occ=0.0)
L2 = draw_latents(p2, N, rng)
r2 = [metrics(T, p2, L2)["rmst"] for T in T_GRID]
t2 = float(T_GRID[int(np.argmax(r2))])
pc2_pass = t2 == 0.0
out["PC2"] = dict(T_star=t2, rmst_at_0=r2[0], rmst_at_24=r2[-1], passed=bool(pc2_pass))
print(f"PC2 p_occ=0: T* = {t2:.0f} mo (must be 0)   RMST 0mo={r2[0]:.1f} 24mo={r2[-1]:.1f}   "
      f"{'PASS' if pc2_pass else 'FAIL'}")

# ---------------- PC3: no device hazard -> T* must run to horizon ----------------
p3 = replace(base_cal, dev_rate=1e-9)
L3 = draw_latents(p3, N, rng)
r3 = [metrics(T, p3, L3)["rmst"] for T in T_GRID]
t3 = float(T_GRID[int(np.argmax(r3))])
pc3_pass = t3 == float(T_GRID[-1])
out["PC3"] = dict(T_star=t3, rmst_at_0=r3[0], rmst_at_24=r3[-1], passed=bool(pc3_pass))
print(f"PC3 dev_rate=0: T* = {t3:.0f} mo (must be {T_GRID[-1]:.0f})   RMST 0mo={r3[0]:.1f} "
      f"24mo={r3[-1]:.1f}   {'PASS' if pc3_pass else 'FAIL'}")

out["all_passed"] = bool(pc1_pass and pc2_pass and pc3_pass)
out["k_calibrated"] = k_cal
json.dump(out, open("results/controls.json", "w"), indent=2)
print("\nALL CONTROLS PASSED" if out["all_passed"] else "\nCONTROL FAILURE - result not reportable")
