"""
Phase F (architecture comparison) and Phase G (trafficking attack), run as one experiment.

Pre-registered in PHASE2_PREREG.md sections 6 and 7.1, committed before any result here existed.
This runs ONLY if the Phase E positive control passed; the pass set is read from disk and the
script refuses to run otherwise (PHASE2_PREREG.md section 9, stopping rule 1).

Phase G is not deferred to the end. TRAFFICKING_PREREG.md already established, before any
alternative was run, that absolute volume exclusion is probably the wrong rule: the 11 um/min
T-cell speed that calibrates the model was itself measured in densely cellular lymph-node cortex,
and at swap_prob = 0.5 the follicle engaged fraction rises from 0.048 to 0.576 - essentially the
dispersed value - with T cells reaching 93% of the way to the follicle centre without any killing.
So the trafficking rule is run as a first-class axis, not as an afterthought.

DECISION RULE, unchanged from TRAFFICKING_PREREG.md T2: if the architecture x schedule interaction
shrinks by more than half at swap_prob = 0.5, the effect is classified MODEL ARTEFACT (category B)
and is not presented as a follicular-lymphoma prediction.

usage:  python exp_arch2.py F     # full schedule sweep, absolute exclusion (swap_prob = 0)
        python exp_arch2.py G     # trafficking sweep at swap_prob 0.5 and 1.0
        python exp_arch2.py A     # analysis only, from stored results
"""
import json, os, sys, time
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lymphoid import Lymphoid
import exp_schedule as S
import schedules as SCH
from exp_poscontrol import load_combos, N_T_OVERRIDE as _NT_OVERRIDE
from exp_poscontrol4 import run as _run4, combos as _combos4

OUT = S.OUT
ARCHS = ['dispersed', 'multi', 'follicle']
SEEDS = list(range(8))
# Phase G uses a reduced schedule set for cost: continuous, the shortest repeating interval, the
# Monday-Friday regimen, and the clinical 7-day interval. Chosen before any Phase F result was
# examined, on the grounds that these span duty cycle and interval length.
G_SCHEDULES = ['A_cont', 'B_6on1off', 'B_MO_FR', 'A_tfi7']


def run(args):
    arch, sched_name, rep_label, params, influx, seed, swap, n_t = args
    p = dict(params); p['swap_prob'] = swap
    r = _run4((arch, sched_name, rep_label, p, influx, seed, n_t))
    r['swap_prob'] = swap
    return r


def passing_combos():
    """Gate on the NEWEST positive-control verdict that exists. Phase 2 wrote expE_*, Phase 3
    expE3_*, Phase 4 expE4_*; all are preserved. Returns (label, params, influx, n_t) so the
    initial T-cell number travels with the regime instead of sitting in a module global."""
    for fn, phase in (('expE4_poscontrol_verdict.json', 4),
                      ('expE3_poscontrol_verdict.json', 3),
                      ('expE_poscontrol_verdict.json', 2)):
        if os.path.exists(f'{OUT}/{fn}'):
            break
    else:
        raise SystemExit('no positive-control verdict found')
    print(f'gating on {fn} (phase {phase})')
    v = json.load(open(f'{OUT}/{fn}'))
    if not v['passes']:
        print('POSITIVE CONTROL FAILED - MODEL NOT VALIDATED FOR SCHEDULING QUESTION')
        print('PHASE2_PREREG.md section 9, stopping rule 1: the architecture comparison is NOT run.')
        sys.exit(1)
    keys = {(p[0], float(p[1])) for p in v['passes']}
    if phase == 4:
        allc = [(l, p, i, n) for l, p, i, n in _combos4()]
    else:
        allc = [(l, p, i, _NT_OVERRIDE if _NT_OVERRIDE is not None else S.N_T)
                for l, p, i in load_combos()]
    combos = [c for c in allc if (c[0], c[2]) in keys]
    print(f'Positive control passed in {len(keys)} combination(s); running those only.')
    for label, _p, infl, n_t in combos:
        print(f'    [{label}] influx {infl:.1e}  N_T {n_t}')
    return combos


def launch(jobs, out_file, tag):
    print(f'{tag}: {len(jobs)} runs', flush=True)
    t0 = time.time(); res = []
    with Pool(4) as pool:
        for i, r in enumerate(pool.imap_unordered(run, jobs, chunksize=1), 1):
            res.append(r)
            if i % 20 == 0 or i == len(jobs):
                el = time.time() - t0
                print(f'  {i}/{len(jobs)}  {el/60:.1f} min, eta {(el/i*(len(jobs)-i))/60:.0f} min',
                      flush=True)
                json.dump(res, open(f'{OUT}/{out_file}.partial', 'w'))
    json.dump(res, open(f'{OUT}/{out_file}', 'w'))
    return res


# --------------------------------------------------------------------------- analysis
def benefit(res, arch, sched, rep, infl, swap):
    """Per-seed benefit of a schedule over continuous, as a fraction of the initial burden.
    Positive means the schedule left LESS disease than continuous dosing."""
    sub = {r['seed']: r for r in res
           if r['arch'] == arch and r['rep'] == rep and r['influx'] == infl
           and r.get('swap_prob', 0.0) == swap}
    out = {}
    for s in SEEDS:
        a = [r for r in res if r['arch'] == arch and r['rep'] == rep and r['influx'] == infl
             and r.get('swap_prob', 0.0) == swap and r['seed'] == s and r['sched'] == sched]
        c = [r for r in res if r['arch'] == arch and r['rep'] == rep and r['influx'] == infl
             and r.get('swap_prob', 0.0) == swap and r['seed'] == s and r['sched'] == SCH.CONTINUOUS]
        if a and c:
            out[s] = (c[0]['nB42'] - a[0]['nB42']) / float(c[0]['n0'])
    return out


def cliffs_delta(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if not len(x) or not len(y):
        return float('nan')
    gt = sum((xi > y).sum() for xi in x); lt = sum((xi < y).sum() for xi in x)
    return (gt - lt) / float(len(x) * len(y))


def holm(pvals):
    order = np.argsort(pvals); adj = np.empty(len(pvals)); run_max = 0.0
    for rank, i in enumerate(order):
        v = min(1.0, (len(pvals) - rank) * pvals[i]); run_max = max(run_max, v); adj[i] = run_max
    return adj


def analyse(res, swaps):
    from scipy import stats
    combos = sorted({(r['rep'], r['influx']) for r in res})
    report = []
    for rep, infl in combos:
        for swap in swaps:
            names = sorted({r['sched'] for r in res if r.get('swap_prob', 0.0) == swap})
            if not names:
                continue
            print(f'\n{"="*82}')
            print(f'[{rep}] influx {infl:.1e}   swap_prob = {swap}')
            print(f'{"="*82}')
            print(f'  {"schedule":12s} {"duty":>5s} | '
                  + ' '.join(f'{a[:9]:>10s}' for a in ARCHS)
                  + f' | {"disp-foll":>10s} {"p":>8s} {"p_holm":>8s} {"cliff":>6s}')
            rows, praw = [], []
            for name in [n for n, _o, _f in SCH.ALL if n in names]:
                on, off = SCH.by_name(name)
                b = {a: benefit(res, a, name, rep, infl, swap) for a in ARCHS}
                med = {a: (np.median(list(b[a].values())) if b[a] else float('nan')) for a in ARCHS}
                if name == SCH.CONTINUOUS:
                    print(f'  {name:12s} {SCH.duty_cycle(on,off):5.2f} | '
                          + ' '.join(f'{med[a]:10.4f}' for a in ARCHS) + ' | (reference)')
                    continue
                d, f = list(b['dispersed'].values()), list(b['follicle'].values())
                p = stats.mannwhitneyu(d, f, alternative='two-sided').pvalue if d and f else 1.0
                rows.append((name, on, off, med, med['dispersed'] - med['follicle'],
                             p, cliffs_delta(d, f)))
                praw.append(p)
            padj = holm(np.array(praw)) if praw else []
            for (name, on, off, med, inter, p, cd), pa in zip(rows, padj):
                print(f'  {name:12s} {SCH.duty_cycle(on,off):5.2f} | '
                      + ' '.join(f'{med[a]:10.4f}' for a in ARCHS)
                      + f' | {inter:10.4f} {p:8.4f} {pa:8.4f} {cd:6.2f}')
                report.append(dict(rep=rep, influx=infl, swap=swap, sched=name,
                                   duty=SCH.duty_cycle(on, off),
                                   **{f'benefit_{a}': (None if med[a] != med[a] else float(med[a]))
                                      for a in ARCHS},
                                   interaction=float(inter), p=float(p), p_holm=float(pa),
                                   cliffs_delta=float(cd)))
    return report


def trafficking_verdict(report):
    print(f'\n{"="*82}\nPHASE G DECISION RULE (TRAFFICKING_PREREG.md T2)\n{"="*82}')
    base = {(r['rep'], r['influx'], r['sched']): r for r in report if r['swap'] == 0.0}
    any05 = [r for r in report if r['swap'] == 0.5]
    if not any05:
        print('  swap_prob = 0.5 not yet run.')
        return None
    print(f'  {"combination":34s} {"|inter| swap0":>14s} {"|inter| swap0.5":>16s} {"ratio":>7s}  verdict')
    shrunk, kept = 0, 0
    for r in any05:
        k = (r['rep'], r['influx'], r['sched'])
        if k not in base:
            continue
        i0, i5 = abs(base[k]['interaction']), abs(r['interaction'])
        ratio = i5 / i0 if i0 > 0 else float('nan')
        v = 'shrinks >half -> ARTEFACT' if (i0 > 0 and ratio < 0.5) else 'survives'
        shrunk += int(i0 > 0 and ratio < 0.5); kept += int(not (i0 > 0 and ratio < 0.5))
        print(f'  [{r["rep"]}|{r["influx"]:.0e}|{r["sched"]:11s}] {i0:14.4f} {i5:16.4f} '
              f'{ratio:7.2f}  {v}')
    print(f'\n  {shrunk} of {shrunk+kept} comparisons shrink by more than half at swap_prob = 0.5.')
    if shrunk > kept:
        print('  VERDICT: the architecture x schedule interaction is MODEL-DEPENDENT.')
        print('  Classification B (model artefact). It is NOT presented as an FL prediction.')
    else:
        print('  VERDICT: the interaction survives the trafficking attack at swap_prob = 0.5.')
    return dict(shrunk=shrunk, kept=kept)


if __name__ == '__main__':
    mode = (sys.argv[1] if len(sys.argv) > 1 else 'A').upper()
    if mode == 'F':
        combos = passing_combos()
        jobs = [(a, n, lab, p, infl, s, 0.0, n_t)
                for lab, p, infl, n_t in combos for a in ARCHS
                for n, _o, _f in SCH.ALL for s in SEEDS]
        res = launch(jobs, 'expF_arch.json', 'Phase F, architecture comparison, swap_prob = 0')
        rep = analyse(res, [0.0])
        json.dump(rep, open(f'{OUT}/expF_arch_report.json', 'w'), indent=1)
    elif mode == 'G':
        combos = passing_combos()
        jobs = [(a, n, lab, p, infl, s, sw, n_t)
                for lab, p, infl, n_t in combos for a in ARCHS
                for n in G_SCHEDULES for s in SEEDS for sw in (0.5, 1.0)]
        res = launch(jobs, 'expG_traffic.json', 'Phase G, trafficking attack')
        both = json.load(open(f'{OUT}/expF_arch.json')) + res
        rep = analyse(both, [0.0, 0.5, 1.0])
        json.dump(rep, open(f'{OUT}/expG_traffic_report.json', 'w'), indent=1)
        trafficking_verdict(rep)
    else:
        res = []
        for fn in ('expF_arch.json', 'expG_traffic.json'):
            if os.path.exists(f'{OUT}/{fn}'):
                res += json.load(open(f'{OUT}/{fn}'))
        swaps = sorted({r.get('swap_prob', 0.0) for r in res})
        rep = analyse(res, swaps)
        if 0.5 in swaps:
            trafficking_verdict(rep)
