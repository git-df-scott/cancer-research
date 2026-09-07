"""Analysis of experiment L1 against pre-registered P1-P5. Paired by seed within architecture."""
import json, numpy as np
from scipy.stats import wilcoxon, spearmanr
R = json.load(open('results/expL1.json'))
G = {}
for r in R: G.setdefault((r['arch'], r['tfi']), {})[r['seed']] = r
ARCHS = ['follicle', 'multi', 'dispersed']; TFIS = [0, 2, 4, 7, 14]

print('=== Descriptives (median over 12 seeds) ===')
print(f"{'arch':11s}{'TFI':>4}{'nB d28':>8}{'nB d42':>8}{'clr':>5}{'engFrac28':>11}{'meanE28':>9}{'contact/T':>11}{'pen':>6}")
for a in ARCHS:
    for k in TFIS:
        rs = list(G[(a, k)].values())
        med = lambda f: np.median([x[f] for x in rs])
        pen = [x['penetration'] for x in rs if x['penetration'] is not None]
        print(f"{a:11s}{k:>4}{med('nB28'):>8.0f}{med('nB42'):>8.0f}{sum(x['cleared'] for x in rs):>5}"
              f"{med('engFrac28'):>11.3f}{med('meanE28'):>9.4f}{med('contact_per_T'):>11.0f}"
              f"{(np.median(pen) if pen else float('nan')):>6.2f}")

print('\n=== P1 engaged fraction, continuous arm: follicle >=3x lower than dispersed ===')
ef = {a: np.median([x['engFrac28'] for x in G[(a, 0)].values()]) for a in ARCHS}
print('  ', {a: round(v, 3) for a, v in ef.items()}, f"| ratio dispersed/follicle = {ef['dispersed']/ef['follicle']:.1f}",
      '| P1', 'HELD' if ef['dispersed'] / ef['follicle'] >= 3 and ef['follicle'] < ef['multi'] < ef['dispersed'] else 'FAILED')

print('\n=== P2 exhaustion at d28, continuous arm: follicle >=2x lower than dispersed ===')
ee = {a: np.median([x['meanE28'] for x in G[(a, 0)].values()]) for a in ARCHS}
print('  ', {a: round(v, 4) for a, v in ee.items()}, f"| ratio = {ee['dispersed']/max(ee['follicle'],1e-9):.1f}",
      '| P2', 'HELD' if ee['dispersed'] / max(ee['follicle'], 1e-9) >= 2 else 'FAILED')

print('\n=== P3/P4 TFI benefit vs continuous, paired on day-42 burden (negative = TFI better) ===')
best = {}
for a in ARCHS:
    print(f'  -- {a} --')
    rows = []
    for k in TFIS[1:]:
        c, t = G[(a, 0)], G[(a, k)]
        seeds = sorted(set(c) & set(t))
        d = np.array([t[s]['nB42'] - c[s]['nB42'] for s in seeds], float)
        p = wilcoxon(d).pvalue if np.any(d != 0) else 1.0
        win = (np.median(d) < 0) and (p < 0.05)
        rows.append((k, np.median(d), p, win, int((d < 0).sum()), len(d)))
        print(f"     TFI{k:>3}: median diff {np.median(d):>8.0f}  better in {int((d<0).sum())}/{len(d)}  p={p:.3f}  {'BEATS CONTINUOUS' if win else ''}")
    anywin = any(r[3] for r in rows)
    bestrow = min(rows, key=lambda r: r[1])
    best[a] = bestrow[1]
    print(f"     -> any TFI beats continuous: {anywin}; best median diff {bestrow[1]:.0f} (TFI{bestrow[0]})")
print(f"\n  P3: dispersed has a winner AND follicle does not ->",
      'HELD' if any(wilcoxon(np.array([G[('dispersed',k)][s]['nB42']-G[('dispersed',0)][s]['nB42'] for s in sorted(G[('dispersed',0)])],float)).pvalue<0.05
                    and np.median([G[('dispersed',k)][s]['nB42']-G[('dispersed',0)][s]['nB42'] for s in sorted(G[('dispersed',0)])])<0 for k in TFIS[1:])
      and not any(wilcoxon(np.array([G[('follicle',k)][s]['nB42']-G[('follicle',0)][s]['nB42'] for s in sorted(G[('follicle',0)])],float)).pvalue<0.05
                  and np.median([G[('follicle',k)][s]['nB42']-G[('follicle',0)][s]['nB42'] for s in sorted(G[('follicle',0)])])<0 for k in TFIS[1:]) else 'FAILED')
print(f"  P4 ordering dispersed < multi < follicle (more negative = more benefit): "
      f"{ {a: round(best[a]) for a in ARCHS} } -> "
      f"{'HELD' if best['dispersed'] < best['multi'] < best['follicle'] else 'FAILED'}")

print('\n=== P5 TFI benefit vs continuous-arm engaged fraction, across 12 arch-TFI cells ===')
xs, ys = [], []
for a in ARCHS:
    for k in TFIS[1:]:
        c, t = G[(a, 0)], G[(a, k)]
        seeds = sorted(set(c) & set(t))
        xs.append(ef[a]); ys.append(np.median([t[s]['nB42'] - c[s]['nB42'] for s in seeds]))
rho = spearmanr(xs, ys)
print(f"  Spearman(engaged fraction, day-42 diff) = {rho.correlation:.2f} (p={rho.pvalue:.3f})")
print(f"  more engagement should mean MORE benefit i.e. more negative diff -> expect rho < -0.6:",
      'HELD' if rho.correlation < -0.6 else 'FAILED')
