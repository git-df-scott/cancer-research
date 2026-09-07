"""E4/E5 + decision rules R1 and R3.

E4: fraction of plausible parameter space with an interior optimum. For every draw,
    k is re-calibrated to the same external 9-month anchor, so the explored space is
    exactly the space consistent with the one thing actually observed. Draws where the
    anchor is unreachable are discarded, not forced.
E5: one-at-a-time sensitivity of T*.
R3: does the result invert under Weibull unmasking with shape > 1 (growth-like,
    non-memoryless) rather than exponential?
"""
import json, numpy as np
from dataclasses import replace, asdict
from bridge_model import Params, draw_latents, metrics, calibrate_k

T_GRID = np.arange(0, 25, 1.0)
N_DRAW = 60_000

def evaluate(p, rng, n=N_DRAW):
    L = draw_latents(p, n, rng)
    k = calibrate_k(p, L)
    if k is None:
        return None
    pc = replace(p, k=k)
    m = [metrics(T, pc, L) for T in T_GRID]
    rmst = np.array([x["rmst"] for x in m])
    a60 = np.array([x["alive60"] for x in m])
    fut = np.array([x["futile_tx"] for x in m])
    i = int(np.argmax(rmst))
    return dict(k=k, T_star=float(T_GRID[i]), rmst_max=float(rmst[i]),
                rmst_gap=float(rmst.max() - rmst.min()),
                T_star_alive60=float(T_GRID[int(np.argmax(a60))]),
                futile_at_0=float(fut[0]), futile_at_6=float(fut[6]),
                futile_at_12=float(fut[12]),
                interior=bool(0.0 < T_GRID[i] < T_GRID[-1]))

# ---------------- E4: plausible parameter space ----------------
rng = np.random.default_rng(2015)
draws, discarded = [], 0
for _ in range(250):
    p = Params(p_occ=rng.uniform(0.35, 0.85), med_unmask=rng.uniform(3.0, 12.0),
               dev_rate=rng.uniform(0.010, 0.050), med_met=rng.uniform(4.0, 9.0))
    r = evaluate(p, rng)
    if r is None:
        discarded += 1
        continue
    r["draw"] = {kk: float(vv) for kk, vv in asdict(p).items()}
    draws.append(r)

Ts = np.array([d["T_star"] for d in draws])
gaps = np.array([d["rmst_gap"] for d in draws])
interior = np.array([d["interior"] for d in draws])
a60T = np.array([d["T_star_alive60"] for d in draws])

print(f"E4  {len(draws)} usable draws ({discarded} discarded: 9-month anchor unreachable)")
print(f"    interior optimum in {interior.mean()*100:.1f}% of draws   "
      f"(R1 threshold 50%)  -> {'SURVIVES' if interior.mean()>=0.5 else 'FAILS'}")
print(f"    T* (RMST-60)     median {np.median(Ts):.0f} mo, IQR "
      f"{np.percentile(Ts,25):.0f}-{np.percentile(Ts,75):.0f}, range {Ts.min():.0f}-{Ts.max():.0f}")
print(f"    RMST best-worst gap: median {np.median(gaps):.2f} mo, "
      f"{(gaps<1.0).mean()*100:.0f}% of draws below the 1.0-mo R2 negligibility threshold")
print(f"    T* by P(alive at 60mo): {(a60T==0).mean()*100:.0f}% of draws put it at T=0")
print(f"    futile transplant rate: {np.mean([d['futile_at_0'] for d in draws]):.3f} at T=0  ->  "
      f"{np.mean([d['futile_at_6'] for d in draws]):.3f} at 6mo  ->  "
      f"{np.mean([d['futile_at_12'] for d in draws]):.3f} at 12mo")

# ---------------- E5: one-at-a-time sensitivity ----------------
print("\nE5  one-at-a-time sensitivity of T* (RMST-60)")
base = Params()
for name, vals in [("p_occ", [0.35,0.5,0.65,0.8]), ("med_unmask", [3,6,9,12]),
                   ("dev_rate", [0.010,0.020,0.030,0.045]), ("med_met", [4,6,9])]:
    row = []
    for v in vals:
        r = evaluate(replace(base, **{name: v}), np.random.default_rng(7))
        row.append(f"{v:>6}:T*={r['T_star']:>2.0f}(gap {r['rmst_gap']:.2f})" if r else f"{v:>6}: n/a")
    print(f"    {name:<12} " + "  ".join(row))

# ---------------- R3: Weibull unmasking ----------------
print("\nR3  non-memoryless unmasking (Weibull shape > 1)")
r3 = {}
for shape in [1.0, 1.5, 2.0, 3.0]:
    r = evaluate(replace(base, unmask_shape=shape), np.random.default_rng(31))
    r3[shape] = r
    print(f"    shape {shape:>3.1f}: T* = {r['T_star']:>2.0f} mo, gap {r['rmst_gap']:.2f} mo, "
          f"k={r['k']:.2f}, T*(alive60)={r['T_star_alive60']:.0f}")
inv = r3[1.0]["interior"] != r3[2.0]["interior"]
print(f"    result inverts under shape=2.0? {inv}  -> "
      f"{'MODEL ARTEFACT (R3)' if inv else 'R3 passed'}")

json.dump(dict(n_draws=len(draws), discarded=discarded,
               interior_fraction=float(interior.mean()),
               T_star_median=float(np.median(Ts)),
               T_star_iqr=[float(np.percentile(Ts,25)), float(np.percentile(Ts,75))],
               gap_median=float(np.median(gaps)),
               frac_gap_below_1mo=float((gaps<1.0).mean()),
               frac_alive60_prefers_T0=float((a60T==0).mean()),
               weibull={str(k): v for k, v in r3.items()},
               draws=draws), open("results/robustness.json","w"), indent=2)
