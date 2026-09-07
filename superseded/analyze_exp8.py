"""Analysis of experiment 8 (follicular architecture) against pre-registered F1-F5."""
import json, numpy as np
from scipy.stats import wilcoxon, mannwhitneyu, spearmanr
H = 6000; VF = 0.0472
runs = json.load(open('results/exp8.json'))
G = {}
for r in runs: G.setdefault((r['M'], r['apoptosis'], r['kills_quiescent'], r['arm']), {})[r['seed']] = r

def paired(M, a, kq):
    c, t = G[(M, a, kq, 'continuous')], G[(M, a, kq, 'threshold')]
    seeds = sorted(set(c) & set(t))
    ben = np.array([(t[s]['ttr'] or H) - (c[s]['ttr'] or H) for s in seeds])
    lag = np.array([(t[s]['t0'] or np.nan) - (c[s]['t0'] or np.nan) for s in seeds if t[s]['t0'] and c[s]['t0']])
    vr = np.array([t[s]['v'] / c[s]['v'] for s in seeds if t[s]['v'] and c[s]['v']])
    ttr_c = np.array([c[s]['ttr'] or H for s in seeds])
    p = wilcoxon(ben).pvalue if (ben != 0).any() else 1.0
    return dict(ben=ben, lag=lag, vr=vr, ttr_c=ttr_c, p=p, better=int((ben > 0).sum()), n=len(seeds),
                depth=np.median([t[s]['depth'] for s in seeds]), hyp=np.median([t[s]['hypoxic0'] for s in seeds]),
                quies=np.mean([t[s]['quiescent'] for s in seeds]), D=np.median([dutyc(t[s]) for s in seeds]),
                vt=np.median([t[s]['v'] for s in seeds if t[s]['v']]) / VF, vc=np.median([c[s]['v'] for s in seeds if c[s]['v']]) / VF,
                erad_c=sum(c[s]['eradicated'] for s in seeds))

def dutyc(r):
    tr = np.array(r['traj']); m = (tr[:, 2] >= 100) & (tr[:, 2] <= 2500)
    return float((tr[m, 5] > 0).mean()) if m.sum() else np.nan

S = {}
print(f"{'M':>3s}{'a':>7s}{'drug':>9s}{'n':>3s}{'depth':>6s}{'hyp0':>6s}{'qui':>5s}{'TTRc':>6s}{'benefit':>8s}{'IQR':>13s}{'better':>7s}{'p':>9s}{'lag':>6s}{'v_t/vf':>7s}{'v_c/vf':>7s}{'D':>5s}{'claim':>6s}")
for a in (0.0005, 0.01):
    for kq in (False, True):
        for M in (1, 4, 9, 16):
            s = paired(M, a, kq); S[(M, a, kq)] = s
            q = np.percentile(s['ben'], [25, 75])
            claim = np.median(s['ben']) > 0 and s['p'] < 0.05
            print(f"{M:3d}{a:7.4f}{'nonspec' if kq else 'cycle':>9s}{s['n']:3d}{s['depth']:6.1f}{s['hyp']:6.2f}{s['quies']:5.1f}{np.median(s['ttr_c']):6.0f}{np.median(s['ben']):8.0f}  [{q[0]:5.0f},{q[1]:5.0f}]{s['better']:4d}/{s['n']:<3d}{s['p']:9.1e}{np.median(s['lag']) if len(s['lag']) else np.nan:6.0f}{s['vt']:7.2f}{s['vc']:7.2f}{s['D']:5.2f}{'yes' if claim else 'NO':>6s}")

print('\n--- F1 jail: lag gain vs M, a=0.0005, cycle-specific')
lags = [np.median(S[(M, 0.0005, False)]['lag']) for M in (1, 4, 9, 16)]
print('  lag medians M=1,4,9,16:', [f'{x:.0f}' for x in lags], '| monotone:', all(np.diff(lags) <= 0), '| >100 at M=1:', lags[0] > 100, '| <40 at M=16:', lags[-1] < 40)
hyp = []; lg = []
for (M, a, kq), s in S.items():
    t = G[(M, a, kq, 'threshold')]; c = G[(M, a, kq, 'continuous')]
    for sd in t:
        if sd in c and t[sd]['t0'] and c[sd]['t0']: hyp.append(t[sd]['hypoxic0']); lg.append(t[sd]['t0'] - c[sd]['t0'])
rho = spearmanr(hyp, lg)
print(f'  Spearman(hypoxic0, lag gain) over {len(hyp)} pairs: {rho.correlation:.2f} (p={rho.pvalue:.1e}) | >0.6:', rho.correlation > 0.6)

print('\n--- F2 confinement ∝ depth, a=0.0005, kills-quiescent drug')
b1 = np.median(S[(1, 0.0005, True)]['ben']); d1 = S[(1, 0.0005, True)]['depth']
for M in (4, 9, 16):
    s = S[(M, 0.0005, True)]; rb = np.median(s['ben']) / b1; rd = s['depth'] / d1
    print(f'  M={M:2d}: benefit ratio {rb:.2f}, depth ratio {rd:.2f}, ratio/ratio {rb/rd:.2f} | within 1.5x:', 1/1.5 <= rb/rd <= 1.5)
for M in (1, 4, 9, 16):
    s = S[(M, 0.0005, True)]; print(f'  M={M:2d}: benefit per site of depth = {np.median(s["ben"])/s["depth"]:.1f} steps/site (law ~10.5)')
dep = []; bb = []
for M in (1, 4, 9, 16):
    t = G[(M, 0.0005, True, 'threshold')]; c = G[(M, 0.0005, True, 'continuous')]
    for sd in t: dep.append(t[sd]['depth']); bb.append((t[sd]['ttr'] or H) - (c[sd]['ttr'] or H))
print(f'  Spearman(depth, benefit), 96 pairs: {spearmanr(dep, bb).correlation:.2f}')

print('\n--- F3 drug-class effect by M, a=0.0005')
for M in (1, 4, 9, 16):
    bs = S[(M, 0.0005, False)]['ben']; bn = S[(M, 0.0005, True)]['ben']
    diff = np.median(bs) - np.median(bn); p = mannwhitneyu(bs, bn).pvalue
    print(f'  M={M:2d}: cycle {np.median(bs):5.0f}  nonspec {np.median(bn):5.0f}  lost {diff:5.0f}  MWU p={p:.1e}')
print('  M=1 lost >150:', np.median(S[(1,0.0005,False)]['ben'])-np.median(S[(1,0.0005,True)]['ben']) > 150, '| M=16 lost <50:', abs(np.median(S[(16,0.0005,False)]['ben'])-np.median(S[(16,0.0005,True)]['ben'])) < 50)

print('\n--- F4 aggressive regime a=0.01, M=16')
for kq in (False, True):
    s = S[(16, 0.01, kq)]; claim = np.median(s['ben']) > 0 and s['p'] < 0.05
    print(f"  {'nonspec' if kq else 'cycle':8s}: benefit {np.median(s['ben']):5.0f} ({s['better']}/{s['n']}, p={s['p']:.2e}), TTR_cont {np.median(s['ttr_c']):.0f}, claim rule {'passes' if claim else 'FAILS'} | F4 holds:", (not claim) or np.median(s['ben']) < 60)

print('\n--- F5 cure (continuous drug, mutation on)')
try:
    cure = json.load(open('results/exp8_cure.json'))
    for kq in (False, True):
        for M in (1, 16):
            rs = [r for r in cure if r['M'] == M and r['kills_quiescent'] == kq]
            er = [r for r in rs if r['eradicated']]
            te = np.median([r['t_erad'] for r in er]) if er else np.nan
            print(f"  {'nonspec' if kq else 'cycle':8s} M={M:2d}: eradicated {len(er)}/{len(rs)}  median t_erad {te:.0f}  r_start>0: {sum(r['r_start']>0 for r in rs)}  median TTR of non-eradicated {np.median([r['ttr'] or H for r in rs if not r['eradicated']]) if len(er)<len(rs) else np.nan:.0f}")
    e = lambda M, kq: sum(r['eradicated'] for r in cure if r['M'] == M and r['kills_quiescent'] == kq)
    t = lambda M, kq: np.median([r['t_erad'] for r in cure if r['M'] == M and r['kills_quiescent'] == kq and r['eradicated']] or [np.nan])
    print(f'  cycle: M16-M1 = {e(16,False)-e(1,False)} (>=6?), t ratio {t(16,False)/t(1,False):.2f} (<0.7?) | nonspec M1 - cycle M1 = {e(1,True)-e(1,False)} (>=4?) | M effect nonspec = {e(16,True)-e(1,True)}')
except FileNotFoundError:
    print('  (cure file not yet written)')
