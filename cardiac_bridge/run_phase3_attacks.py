"""Sections G (attacks on the synergy) and H (immunotherapy window: required vs observed).
Pre-registered. MCID = 5.0 pp."""
import json, numpy as np
from dataclasses import replace
from phase3_model import P3, draw_latents3, metrics3, calibrate_k3

N = 250_000
T_GRID = np.arange(0, 25, 1.0)
MCID = 0.05
rng = np.random.default_rng(2024)
base = P3()
L = draw_latents3(base, N, rng)
P = replace(base, k=calibrate_k3(base, L))

def benefit(p, Lx=None):
    Lx = L if Lx is None else Lx
    a = np.array([metrics3(float(T), p, Lx)["alive60"] for T in T_GRID])
    i = int(np.argmax(a))
    return float(a[i] - a[0]), float(T_GRID[i])

# ---------------- G7: is the "synergy" biology or the observation process? -------------
print("G7  Phase 2 synergy: true biology, or an artefact of cytostasis delaying detection?")
print("    tx_stasis = clearance window (BIOLOGY). s_detect = detection clock (OBSERVATION).\n")
print(f"{'condition':<44}{'benefit (pp)':>14}{'T*':>6}")
g7 = {}
for label, kw in [("neither (s=1, detect=1)",                dict(tx_stasis=1.0, s_detect=1.0)),
                  ("BIOLOGY only (s=3, detect=1)",           dict(tx_stasis=3.0, s_detect=1.0)),
                  ("OBSERVATION only (s=1, detect=3)",       dict(tx_stasis=1.0, s_detect=3.0)),
                  ("both, as a real drug would (s=3,det=3)", dict(tx_stasis=3.0, s_detect=3.0))]:
    g, T = benefit(replace(P, p_R=0.30, p_dur=0.70, **kw))
    g7[label] = dict(gain=g, T=T)
    print(f"{label:<44}{g*100:>14.2f}{T:>6.0f}")
b0 = g7["neither (s=1, detect=1)"]["gain"]
bio = g7["BIOLOGY only (s=3, detect=1)"]["gain"] - b0
obs = g7["OBSERVATION only (s=1, detect=3)"]["gain"] - b0
tot = g7["both, as a real drug would (s=3,det=3)"]["gain"] - b0
print(f"\n    biology component    {bio*100:+.2f} pp")
print(f"    observation component{obs*100:+.2f} pp")
print(f"    combined             {tot*100:+.2f} pp")
share = obs / tot if abs(tot) > 1e-9 else float('nan')
print(f"    -> observation process accounts for {share:.0%} of the cytostasis effect")

# ---------------- G3/G4/G5/G6: does any benefit survive the attacks? ------------------
print("\nG   attacks at the most generous evidence-admissible corner (pi = 0.30 x 0.70 = 0.21)")
gen = replace(P, p_R=0.30, p_dur=0.70, m_clear=3.0, a_clear=1.0,
              tx_stasis=3.0, med_unmask=9.0, dev_rate=0.010)
Lg = draw_latents3(gen, N, np.random.default_rng(606))
kg = calibrate_k3(gen, Lg); gen = replace(gen, k=kg)
attacks = {}
for label, kw in [("generous corner, no attack",              {}),
                  ("G3 non-exponential clearance (a=2.5)",    dict(a_clear=2.5)),
                  ("G4 acquired escape (p_dur 0.70 -> 0.35)", dict(p_dur=0.35)),
                  ("G5 mid-quality device centre (0.025)",    dict(dev_rate=0.025)),
                  ("G5 poor device centre (0.045)",           dict(dev_rate=0.045)),
                  ("G6 faster unmasking (IGR: 3 mo)",         dict(med_unmask=3.0)),
                  ("all attacks together",                    dict(a_clear=2.5, p_dur=0.35,
                                                                   dev_rate=0.025, med_unmask=3.0))]:
    p = replace(gen, **kw)
    Lx = Lg if not ({"dev_rate","med_unmask"} & set(kw)) else draw_latents3(p, N, np.random.default_rng(606))
    if {"dev_rate","med_unmask"} & set(kw):
        kk = calibrate_k3(p, Lx)
        p = replace(p, k=kk) if kk else p
    g, T = benefit(p, Lx)
    attacks[label] = dict(gain=g, T=T, meets_mcid=bool(g >= MCID))
    print(f"    {label:<42}{g*100:>7.2f} pp  T*={T:>2.0f}  {'MEETS MCID' if g>=MCID else 'below MCID'}")

# ---------------- H: required efficacy vs observed efficacy ---------------------------
print("\nH   IMMUNOTHERAPY WINDOW: required durable clearance vs what the evidence supports")
print("    Sweeping pi at the single most favourable admissible setting for the hypothesis.")
fav = replace(gen, p_dur=1.0)
rows = []
for pi in [0.05, 0.10, 0.15, 0.21, 0.25, 0.30, 0.40, 0.50]:
    g, T = benefit(replace(fav, p_R=pi), Lg)
    rows.append((pi, g, T))
    print(f"      pi = {pi:.2f}  ->  {g*100:>6.2f} pp   T*={T:>2.0f}   "
          f"{'MEETS MCID' if g >= MCID else ''}")
req = next((pi for pi, g, _ in rows if g >= MCID), None)
print(f"\n    REQUIRED pi (best case for the hypothesis): >= {req}")
print(f"    OBSERVED/admissible pi ceiling             : 0.30 x 0.70 = 0.21")
print(f"    Evidence centre of mass for CARDIAC AS     : 0/1 DART cardiac; ~0 adjuvant chemo")
print(f"    -> regions {'DO NOT OVERLAP' if req is None or req > 0.21 else 'overlap'}")

# MCID sensitivity, reported but NOT the pre-registered decision rule
print("\n    (secondary, not the pre-registered rule) required pi at other MCIDs:")
for m in [0.02, 0.03, 0.05, 0.08]:
    r = next((pi for pi, g, _ in rows if g >= m), None)
    print(f"      MCID {m*100:.0f} pp -> required pi {r if r is not None else '>0.50'}")

json.dump(dict(g7=g7, obs_share=share, attacks=attacks,
               h_rows=[dict(pi=p_, gain=g_, T=t_) for p_, g_, t_ in rows], required_pi=req),
          open("results/phase3_attacks.json","w"), indent=2)
