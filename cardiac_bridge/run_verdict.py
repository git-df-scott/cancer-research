"""Apply the pre-registered decision rules R1-R3 and compute the organ-stewardship
trade-off that the base case surfaced but the pre-registration did not anticipate."""
import json, numpy as np
from scipy import stats

rb = json.load(open("results/robustness.json"))
bc = json.load(open("results/base_case.json"))

n = rb["n_draws"]; f = rb["interior_fraction"]
lo, hi = stats.binomtest(int(round(f*n)), n).proportion_ci(0.95)
print(f"R1  interior optimum in {f*100:.1f}% of {n} draws")
print(f"    95% CI [{lo*100:.1f}%, {hi*100:.1f}%]   threshold 50%")
print(f"    -> CI {'STRADDLES' if lo < 0.5 < hi else 'excludes'} 50%: "
      f"{'not meaningfully cleared' if lo < 0.5 else 'cleared'}")
print(f"R2  median best-worst RMST gap {rb['gap_median']:.2f} mo; "
      f"{rb['frac_gap_below_1mo']*100:.0f}% of draws below 1.0 mo -> not negligible")
w = rb["weibull"]
print(f"R3  T* by Weibull shape: " + ", ".join(f"{s}->{w[s]['T_star']:.0f}mo" for s in w) + " -> passed")
print(f"E2  {rb['frac_alive60_prefers_T0']*100:.0f}% of draws put the cure-maximising T at 0\n")

rows = {r["T"]: r for r in bc["rows"]}
print("Organ stewardship, per 100 patients reaching cardiectomy + TAH:")
print(f"{'bridge':>7} {'organs used':>12} {'futile organs':>14} {'futility rate':>14} {'5yr survivors':>14}")
base_alive = rows[0.0]["alive60"] * 100
for T in [0.0, 3.0, 6.0, 9.0, 12.0]:
    r = rows[T]
    used = r["tx_rate"] * 100
    fut = r["futile_tx"] * 100
    alive = r["alive60"] * 100
    print(f"{T:>7.0f} {used:>12.1f} {fut:>14.1f} {fut/used*100:>13.1f}% {alive:>14.1f}")
r6 = rows[6.0]
d_org = (rows[0.0]["tx_rate"] - r6["tx_rate"]) * 100
d_fut = (rows[0.0]["futile_tx"] - r6["futile_tx"]) * 100
d_alive = (rows[0.0]["alive60"] - r6["alive60"]) * 100
print(f"\n  6-month bridge vs immediate transplant, per 100 patients:")
print(f"    donor organs not consumed : {d_org:.1f}")
print(f"    futile transplants averted: {d_fut:.1f}")
print(f"    5-year survivors lost     : {d_alive:.1f}")
print(f"    organs saved per survivor lost: {d_org/d_alive:.1f}")
json.dump(dict(R1_ci=[lo, hi], R1_cleared=bool(lo >= 0.5), organs_saved_per_100=d_org,
               futile_averted_per_100=d_fut, survivors_lost_per_100=d_alive,
               organs_per_survivor_lost=d_org/d_alive),
          open("results/verdict.json", "w"), indent=2)
