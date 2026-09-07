"""Figures for experiment 8: results/discovery/fig7_exp8.png"""
import json, os, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
H = 6000; VF = 0.0472
runs = json.load(open('results/exp8.json'))
G = {}
for r in runs: G.setdefault((r['M'], r['apoptosis'], r['kills_quiescent'], r['arm']), {})[r['seed']] = r
MS = [1, 4, 9, 16]
os.makedirs('results/discovery', exist_ok=True)

def stats(M, a, kq):
    c, t = G[(M, a, kq, 'continuous')], G[(M, a, kq, 'threshold')]
    seeds = sorted(set(c) & set(t))
    ben = np.array([(t[s]['ttr'] or H) - (c[s]['ttr'] or H) for s in seeds])
    lag = np.array([t[s]['t0'] - c[s]['t0'] for s in seeds if t[s]['t0'] and c[s]['t0']])
    conf = np.array([np.sqrt(2500 / np.pi) * (1 / t[s]['v'] - 1 / c[s]['v']) for s in seeds if t[s]['v'] and c[s]['v']])
    depth = np.array([t[s]['depth'] for s in seeds]); hyp = np.array([t[s]['hypoxic0'] for s in seeds])
    return ben, lag, conf, depth, hyp

fig, ax = plt.subplots(2, 2, figsize=(11, 8.5))
# A: benefit vs M, four conditions
a0 = ax[0, 0]
for a, kq, col, ls, lab in [(0.0005, False, 'C0', '-', 'slow turnover, cycle-specific'), (0.0005, True, 'C0', '--', 'slow turnover, kills quiescent'),
                            (0.01, False, 'C3', '-', 'fast turnover, cycle-specific'), (0.01, True, 'C3', '--', 'fast turnover, kills quiescent')]:
    med = []; lo = []; hi = []
    for M in MS:
        b = stats(M, a, kq)[0]; med.append(np.median(b)); q = np.percentile(b, [25, 75]); lo.append(q[0]); hi.append(q[1])
    a0.errorbar(MS, med, yerr=[np.array(med) - lo, np.array(hi) - med], color=col, ls=ls, marker='o', capsize=3, label=lab)
a0.axhline(0, color='k', lw=0.5); a0.set_xscale('log'); a0.set_xticks(MS); a0.set_xticklabels(MS)
a0.set_xlabel('number of follicles M (same total size)'); a0.set_ylabel('TTR benefit of size-triggered dosing (steps)')
a0.set_title('A  Benefit collapses with follicular architecture'); a0.legend(fontsize=8)
# B: decomposition at slow turnover, cycle-specific
a1 = ax[0, 1]; w = 0.35
for i, M in enumerate(MS):
    ben, lag, conf, depth, hyp = stats(M, 0.0005, False)
    a1.bar(i - w / 2, np.median(conf), w, color='C0', label='confinement gain' if i == 0 else None)
    a1.bar(i + w / 2, np.median(lag), w, color='C1', label='release-lag (jail) gain' if i == 0 else None)
    a1.text(i, max(np.median(conf), np.median(lag)) + 15, f'hyp {np.median(hyp):.2f}\ndepth {np.median(depth):.0f}', ha='center', fontsize=8)
a1.set_xticks(range(4)); a1.set_xticklabels([f'M={M}' for M in MS]); a1.set_ylabel('steps (front-fit decomposition, medians)')
a1.set_title('B  Jail vanishes with hypoxia; confinement tracks depth'); a1.legend(fontsize=8)
# C: benefit vs depth, kills-quiescent drug (no jail), all runs
a2 = ax[1, 0]
for M, mk in zip(MS, 'osD^'):
    ben, lag, conf, depth, hyp = stats(M, 0.0005, True)
    a2.scatter(depth, ben, marker=mk, s=25, alpha=0.7, label=f'M={M}')
d = np.linspace(0, 45, 10); a2.plot(d, 10.5 * d, 'k--', lw=1, label='law: 10.5 steps / site')
a2.set_xlabel('implant depth (sites to nearest empty site)'); a2.set_ylabel('TTR benefit (steps)'); a2.set_title('C  No jail: benefit ∝ jacket depth'); a2.legend(fontsize=8)
# D: cure
a3 = ax[1, 1]
try:
    cure = json.load(open('results/exp8_cure.json'))
    xs = []; ys = []; labs = []
    for kq in (False, True):
        for M in (1, 16):
            rs = [r for r in cure if r['M'] == M and r['kills_quiescent'] == kq]
            ys.append(sum(r['eradicated'] for r in rs)); labs.append(f"M={M}\n{'kills quiescent' if kq else 'cycle-specific'}")
    a3.bar(range(4), ys, color=['C0', 'C0', 'C2', 'C2']); a3.set_xticks(range(4)); a3.set_xticklabels(labs, fontsize=8)
    a3.set_ylabel('tumours eradicated of 24 (continuous drug, mutation on)'); a3.set_title('D  Cure under continuous dosing')
    for i, y in enumerate(ys): a3.text(i, y + 0.3, str(y), ha='center')
except FileNotFoundError:
    a3.text(0.5, 0.5, 'cure runs pending', ha='center', transform=a3.transAxes)
plt.tight_layout(); plt.savefig('results/discovery/fig7_exp8.png', dpi=130)
print('wrote results/discovery/fig7_exp8.png')
