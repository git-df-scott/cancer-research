"""
PHASE 1 analysis of experiment L1. Reports EVERY pre-registered endpoint and prediction,
regardless of outcome, with effect sizes and seed-level uncertainty.

The central quantity is the ARCHITECTURE x SCHEDULE INTERACTION: does the value of a
treatment-free interval differ between architectures? Everything else is supporting detail.
"""
import json, sys, itertools
import numpy as np
from scipy.stats import wilcoxon, mannwhitneyu, spearmanr

SRC = sys.argv[1] if len(sys.argv) > 1 else 'results/expL1.json'
R = json.load(open(SRC))
G = {}
for r in R:
    G.setdefault((r['arch'], r['tfi']), {})[r['seed']] = r
ARCHS = [a for a in ['follicle', 'multi', 'dispersed'] if any(k[0] == a for k in G)]
TFIS = sorted({k[1] for k in G})
SEEDS = sorted(set.intersection(*[set(G[k]) for k in G]))
print(f'source: {SRC}   architectures {ARCHS}   TFIs {TFIS}   seeds n={len(SEEDS)}\n')


def boot_ci(x, f=np.median, n=4000, seed=0):
    x = np.asarray(x, float)
    if len(x) == 0:
        return (np.nan, np.nan)
    rng = np.random.default_rng(seed)
    s = [f(rng.choice(x, len(x), replace=True)) for _ in range(n)]
    return tuple(np.percentile(s, [2.5, 97.5]))


def cliffs_delta(a, b):
    """Non-parametric effect size in [-1,1]. Robust, no distributional assumption."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    if not len(a) or not len(b):
        return np.nan
    gt = sum((x > y) for x in a for y in b)
    lt = sum((x < y) for x in a for y in b)
    return (gt - lt) / (len(a) * len(b))


# ------------------------------------------------------------------ 1. descriptives
print('=' * 110)
print('PHASE 1.1  DESCRIPTIVES  (median [bootstrap 95% CI] over seeds)')
print('=' * 110)
hdr = f"{'arch':10s}{'TFI':>4}{'nB d28':>16}{'nB d42':>16}{'clr':>4}{'kills':>8}{'nT_end':>7}{'engFrac28':>10}{'meanE28':>9}{'contact/T':>10}{'drugdays':>9}{'pen':>6}"
print(hdr)
for a in ARCHS:
    for k in TFIS:
        rs = [G[(a, k)][s] for s in SEEDS if s in G[(a, k)]]
        b28 = [x['nB28'] for x in rs]; b42 = [x['nB42'] for x in rs]
        c28 = boot_ci(b28); c42 = boot_ci(b42)
        pen = [x['penetration'] for x in rs if x.get('penetration') is not None]
        dd = [x.get('total_drug', np.nan) for x in rs]
        print(f"{a:10s}{k:>4}"
              f"{np.median(b28):>7.0f}[{c28[0]:.0f},{c28[1]:.0f}]".ljust(0)[:0] or
              f"{np.median(b28):>7.0f} [{c28[0]:>4.0f},{c28[1]:>4.0f}]"
              f"{np.median(b42):>7.0f} [{c42[0]:>4.0f},{c42[1]:>4.0f}]"
              f"{sum(x['cleared'] for x in rs):>4}"
              f"{np.median([x['kills'] for x in rs]):>8.0f}"
              f"{np.median([x['nT_end'] for x in rs]):>7.0f}"
              f"{np.median([x['engFrac28'] for x in rs]):>10.3f}"
              f"{np.median([x['meanE28'] for x in rs]):>9.4f}"
              f"{np.median([x['contact_per_T'] for x in rs]):>10.0f}"
              f"{(np.median(dd) if not all(np.isnan(dd)) else float('nan')):>9.1f}"
              f"{(np.median(pen) if pen else float('nan')):>6.2f}")

# ------------------------------------------------------------------ 2. preregistered P1-P5
print('\n' + '=' * 110)
print('PHASE 1.2  PRE-REGISTERED PREDICTIONS  (verdicts reported regardless of outcome)')
print('=' * 110)
ef = {a: np.median([G[(a, 0)][s]['engFrac28'] for s in SEEDS if s in G[(a, 0)]]) for a in ARCHS}
ee = {a: np.median([G[(a, 0)][s]['meanE28'] for s in SEEDS if s in G[(a, 0)]]) for a in ARCHS}
ct = {a: np.median([G[(a, 0)][s]['contact_per_T'] for s in SEEDS if s in G[(a, 0)]]) for a in ARCHS}
print('\nP1  engaged fraction under continuous dosing, follicle >=3x lower than dispersed,')
print('    and multi intermediate.')
for a in ARCHS:
    ci = boot_ci([G[(a, 0)][s]['engFrac28'] for s in SEEDS if s in G[(a, 0)]])
    print(f'      {a:10s} engFrac28 = {ef[a]:.3f}  [{ci[0]:.3f}, {ci[1]:.3f}]   contact-min/T = {ct[a]:.0f}')
if 'follicle' in ef and 'dispersed' in ef:
    ratio = ef['dispersed'] / ef['follicle']
    ordered = ('multi' not in ef) or (ef['follicle'] < ef['multi'] < ef['dispersed'])
    print(f'      ratio dispersed/follicle = {ratio:.1f}   ordering holds: {ordered}')
    print(f'      P1 -> {"HELD" if ratio >= 3 and ordered else "FAILED"}')

print('\nP2  mean exhaustion at day 28 under continuous dosing, follicle >=2x lower than dispersed.')
for a in ARCHS:
    ci = boot_ci([G[(a, 0)][s]['meanE28'] for s in SEEDS if s in G[(a, 0)]])
    print(f'      {a:10s} meanE28 = {ee[a]:.4f}  [{ci[0]:.4f}, {ci[1]:.4f}]')
if 'follicle' in ee and 'dispersed' in ee:
    r2 = ee['dispersed'] / max(ee['follicle'], 1e-9)
    print(f'      ratio = {r2:.1f}   P2 -> {"HELD" if r2 >= 2 else "FAILED"}')

print('\nP3  dispersed has >=1 TFI beating continuous (paired, day-42, Wilcoxon p<0.05);')
print('    follicle has none.')
winners = {}
detail = {}
for a in ARCHS:
    rows = []
    for k in [t for t in TFIS if t != 0]:
        c, t = G[(a, 0)], G[(a, k)]
        ss = [s for s in SEEDS if s in c and s in t]
        d = np.array([t[s]['nB42'] - c[s]['nB42'] for s in ss], float)
        p = wilcoxon(d).pvalue if np.any(d != 0) else 1.0
        ci = boot_ci(d)
        win = (np.median(d) < 0) and (p < 0.05)
        rows.append(dict(k=k, med=float(np.median(d)), ci=ci, p=float(p), win=win,
                         nbetter=int((d < 0).sum()), n=len(d), d=d))
    detail[a] = rows
    winners[a] = [r['k'] for r in rows if r['win']]
    print(f'    -- {a} --')
    for r in rows:
        print(f"       TFI{r['k']:>3}: median dnB42 = {r['med']:>8.0f}  [{r['ci'][0]:>7.0f},{r['ci'][1]:>7.0f}]"
              f"  better {r['nbetter']:>2}/{r['n']}  p={r['p']:.4f}  {'<-- beats continuous' if r['win'] else ''}")
    print(f"       winners: {winners[a] if winners[a] else 'none'}")
p3 = bool(winners.get('dispersed')) and not bool(winners.get('follicle'))
print(f'    P3 -> {"HELD" if p3 else "FAILED"}')

print('\nP4  benefit of best TFI ordered dispersed > multi > follicle (more negative = more benefit).')
best = {a: min(r['med'] for r in detail[a]) for a in ARCHS}
for a in ARCHS:
    print(f'      {a:10s} best median dnB42 = {best[a]:>8.0f}')
p4 = ('multi' in best and best['dispersed'] < best['multi'] < best['follicle']) if len(ARCHS) == 3 else None
print(f'    P4 -> {"HELD" if p4 else "FAILED"}')

print('\nP5  TFI benefit correlates with continuous-arm engaged fraction, Spearman rho < -0.6.')
xs, ys = [], []
for a in ARCHS:
    for r in detail[a]:
        xs.append(ef[a]); ys.append(r['med'])
rho = spearmanr(xs, ys)
print(f'      rho = {rho.correlation:.2f} (p={rho.pvalue:.3f}, n={len(xs)} cells)')
print(f'    P5 -> {"HELD" if rho.correlation < -0.6 else "FAILED"}')

# ------------------------------------------------------------------ 3. the interaction
print('\n' + '=' * 110)
print('PHASE 1.3  ARCHITECTURE x SCHEDULE INTERACTION  (the decisive quantity)')
print('=' * 110)
print('Difference-in-differences: (TFI - continuous) in dispersed  MINUS  (TFI - continuous) in follicle.')
print('Burden is normalised to each run\'s own n0 so architectures are comparable on one scale.\n')
print(f"{'TFI':>4}{'disp effect':>26}{'foll effect':>26}{'interaction':>26}{'MWU p':>9}{'cliff d':>9}")
inter_any = False
for k in [t for t in TFIS if t != 0]:
    out = {}
    for a in ('dispersed', 'follicle'):
        if (a, k) not in G:
            continue
        c, t = G[(a, 0)], G[(a, k)]
        ss = [s for s in SEEDS if s in c and s in t]
        out[a] = np.array([(t[s]['nB42'] - c[s]['nB42']) / t[s]['n0'] for s in ss], float)
    if len(out) < 2:
        continue
    dd = out['dispersed']; ff = out['follicle']
    ci_d, ci_f = boot_ci(dd), boot_ci(ff)
    p = mannwhitneyu(dd, ff).pvalue
    delta = cliffs_delta(dd, ff)
    inter = np.median(dd) - np.median(ff)
    if p < 0.05:
        inter_any = True
    print(f"{k:>4}{np.median(dd):>10.4f} [{ci_d[0]:>6.3f},{ci_d[1]:>6.3f}]"
          f"{np.median(ff):>10.4f} [{ci_f[0]:>6.3f},{ci_f[1]:>6.3f}]"
          f"{inter:>26.4f}{p:>9.4f}{delta:>9.2f}")
print(f"\nAny architecture x schedule interaction significant at p<0.05: {inter_any}")

# ------------------------------------------------------------------ 4. schedule ranking
print('\n' + '=' * 110)
print('PHASE 1.4  SCHEDULE RANKING BY ARCHITECTURE (best -> worst on median day-42 burden)')
print('=' * 110)
for a in ARCHS:
    rank = sorted(TFIS, key=lambda k: np.median([G[(a, k)][s]['nB42'] for s in SEEDS if s in G[(a, k)]]))
    vals = {k: np.median([G[(a, k)][s]['nB42'] for s in SEEDS if s in G[(a, k)]]) for k in TFIS}
    print(f"  {a:10s} " + '  >  '.join(f"TFI{k}({vals[k]:.0f})" for k in rank))
print('\n  Obertopp/Basanta (B-ALL, well-mixed) ranking was: short TFIs best, continuous worst.')

# ------------------------------------------------------------------ 5. trajectories
print('\n' + '=' * 110)
print('PHASE 1.5  BURDEN TRAJECTORY (median nB by day, continuous vs TFI7)')
print('=' * 110)
for a in ARCHS:
    for k in (0, 7):
        if (a, k) not in G:
            continue
        rs = [G[(a, k)][s] for s in SEEDS if s in G[(a, k)]]
        days = [7, 14, 21, 28, 35, 42]
        line = []
        for d in days:
            v = []
            for x in rs:
                tr = {t[0]: t[1] for t in x['traj']}
                v.append(tr.get(d, 0))
            line.append(f'd{d}:{np.median(v):.0f}')
        print(f"  {a:10s} TFI{k:<3} " + '  '.join(line))
