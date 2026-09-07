"""Sections G and H on the bug-fixed model. Every parameter is evaluation-time, so all
conditions share ONE latent set (common random numbers) -- the earlier version redrew
between conditions and carried resampling noise into the comparisons."""
import json, numpy as np
from dataclasses import replace
from phase3_model import P3, draw_latents3, metrics3, calibrate_k3

N = 300_000
T_GRID = np.arange(0, 25, 1.0)
MCID = 0.05
L = draw_latents3(P3(), N, np.random.default_rng(1618))
P = replace(P3(), k=calibrate_k3(P3(), L))

def cal(p):
    kk = calibrate_k3(p, L)
    return replace(p, k=kk) if kk else p

def benefit(p):
    a = np.array([metrics3(float(T), cal(p), L)["alive60"] for T in T_GRID])
    i = int(np.argmax(a))
    return float(a[i] - a[0]), float(T_GRID[i])

# ---------------- G7: biology vs observation, on fixed code ----------------
print("G7  is the cytostasis effect biology, or the observation process?")
print("    tx_stasis = clearance window (BIOLOGY).  s_detect = detection clock (OBSERVATION).\n")
g7 = {}
gen0 = replace(P, p_R=0.30, p_dur=0.70)
print(f"{'condition':<42}{'benefit (pp)':>14}{'T*':>5}{'carried@T*':>13}")
for label, kw in [("neither (s=1, detect=1)",      dict(tx_stasis=1.0, s_detect=1.0)),
                  ("BIOLOGY only (s=3, detect=1)", dict(tx_stasis=3.0, s_detect=1.0)),
                  ("OBSERVATION only (s=1,det=3)", dict(tx_stasis=1.0, s_detect=3.0)),
                  ("both (s=3, detect=3)",         dict(tx_stasis=3.0, s_detect=3.0))]:
    p = replace(gen0, **kw); g, T = benefit(p)
    c = metrics3(T, cal(p), L)["carried"]
    g7[label] = dict(gain=g, T=T, carried=c)
    print(f"{label:<42}{g*100:>14.2f}{T:>5.0f}{c:>13.4f}")
b0 = g7["neither (s=1, detect=1)"]["gain"]
bio = g7["BIOLOGY only (s=3, detect=1)"]["gain"] - b0
obs = g7["OBSERVATION only (s=1,det=3)"]["gain"] - b0
tot = g7["both (s=3, detect=3)"]["gain"] - b0
print(f"\n    biology {bio*100:+.2f} pp | observation {obs*100:+.2f} pp | combined {tot*100:+.2f} pp")
print(f"    -> observation share of the effect: {obs/tot if abs(tot)>1e-9 else float('nan'):.0%}")

# ---------------- G1-G6: stacked attacks at the generous corner ----------------
print("\nG   attacks at the most generous evidence-admissible corner")
print("    (p_R=0.30, p_dur=0.70 -> pi=0.21; m_clear=3, a=1, s=3, unmask=9 mo, dev=0.010)\n")
gen = replace(P, p_R=0.30, p_dur=0.70, m_clear=3.0, a_clear=1.0,
              tx_stasis=3.0, med_unmask=9.0, dev_rate=0.010)
attacks = {}
for label, kw in [("no attack (best case)",                   {}),
                  ("G3 non-exponential clearance (a=2.5)",    dict(a_clear=2.5)),
                  ("G4 acquired escape (p_dur 0.70->0.35)",   dict(p_dur=0.35)),
                  ("G5 mid-quality centre (dev 0.025)",       dict(dev_rate=0.025)),
                  ("G5 poor centre (dev 0.045)",              dict(dev_rate=0.045)),
                  ("G6 IGR-anchored unmasking (3 mo)",        dict(med_unmask=3.0)),
                  ("no cytostasis (s=1)",                     dict(tx_stasis=1.0)),
                  ("ALL attacks together",                    dict(a_clear=2.5, p_dur=0.35,
                                                                   dev_rate=0.025, med_unmask=3.0))]:
    g, T = benefit(replace(gen, **kw))
    attacks[label] = dict(gain=g, T=T, meets=bool(g >= MCID))
    print(f"    {label:<42}{g*100:>7.2f} pp  T*={T:>2.0f}  {'MEETS MCID' if g>=MCID else 'below'}")

# ---------------- H: required vs observed ----------------
print("\nH   required durable clearance fraction vs what the evidence supports\n")
fav = replace(gen, p_dur=1.0)
rows = []
for pi in [0.02,0.05,0.07,0.10,0.15,0.21,0.25,0.30,0.40]:
    g, T = benefit(replace(fav, p_R=pi))
    rows.append(dict(pi=pi, gain=g, T=T))
    print(f"      pi={pi:.2f} -> {g*100:>6.2f} pp  T*={T:>2.0f} {'MEETS MCID' if g>=MCID else ''}")
req = next((r["pi"] for r in rows if r["gain"] >= MCID), None)
print(f"\n    REQUIRED pi, single most favourable admissible setting : >= {req}")
print(f"    ADMISSIBLE ceiling  (p_R 0.30 x p_dur 0.70)            : 0.21")
print(f"    CARDIAC-SPECIFIC evidence: DART cardiac 0/1; visceral 1/7 -> pi ~ 0.00-0.07")
print(f"      (1/7 responders x 2/4 durable = 0.07; adjuvant chemo in cardiac sarcoma ~ 0)")
print("\n    (secondary, NOT the pre-registered rule) required pi by MCID:")
for m in [0.02, 0.03, 0.05, 0.08]:
    r = next((x["pi"] for x in rows if x["gain"] >= m), None)
    print(f"      MCID {m*100:.0f} pp -> pi >= {r if r is not None else '>0.40'}")

json.dump(dict(g7=g7, obs_share=(obs/tot if abs(tot)>1e-9 else None), attacks=attacks,
               h_rows=rows, required_pi=req), open("results/phase3_attacks2.json","w"), indent=2)
