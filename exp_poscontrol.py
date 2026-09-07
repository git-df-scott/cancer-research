"""
Phase E: the positive control. Dispersed / leukaemia-like architecture ONLY.

Pre-registered in PHASE2_PREREG.md section 3, before this was run.

WHAT MUST BE REPRODUCED
-----------------------
That an appropriate short, REPEATING treatment-free interval outperforms continuous exposure in
well-mixed leukaemia-like disease. This is the setting where the effect is known, from both the
reference model (Obertopp/Basanta, bioRxiv 2025, PPR1121269) and, at matched cumulative exposure,
from Philipp's own wet-lab data.

PASS CRITERION, FIXED IN ADVANCE
--------------------------------
PC-1 (required): at least one non-continuous schedule gives lower malignant burden at the primary
endpoint than continuous dosing, by a one-sided Wilcoxon signed-rank test over paired seeds at
p < 0.05, AND with a median paired reduction of at least 10% of the continuous arm's burden. The
10% floor exists so a significant-but-negligible difference cannot be counted as a pass.

PC-2 (reported, not required): the winning schedule should be a short repeating interval, off
fraction <= 3 days per cycle. A pass driven only by a long interval is flagged as anomalous.

Every (exhaustion representative x influx regime) combination is run and reported regardless of
outcome. If PC-1 fails in ALL of them, the run STOPS and reports
POSITIVE CONTROL FAILED - MODEL NOT VALIDATED FOR SCHEDULING QUESTION.

Nothing in this file is tuned. The parameters come from calib_philipp.py, which was scored only
against Philipp's measurements, and the admissible influx regimes come from exp_influx_check.py,
which was scored only against Philipp's patient data.
"""
import json, os, sys, time
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lymphoid import Lymphoid
import exp_schedule as S
import schedules as SCH

OUT = S.OUT
DAYS = 42
SEEDS = list(range(8))
ARCH = 'dispersed'
# Phase 3: initial T-cell number selected by the schedule-blind sweep in exp_ntsweep.py, per
# PHASE3_PREREG.md section 3. N_T = 1400 is E:T 1:4, matching Philipp's assay, and was the only
# value admissible in any combination. Phase 2 used exp_schedule.N_T = 200 (E:T 1:28).
# Set N_T_OVERRIDE = None to reproduce the Phase 2 positive control exactly.
N_T_OVERRIDE = 1400
MIN_EFFECT = 0.10      # median paired reduction must be at least this fraction of continuous


def build(arch, seed, params, influx):
    m = Lymphoid(L=S.L, seed=seed, dt=S.DT, p_kill=S.PARAMS['p_kill'],
                 p_div=S.PARAMS['p_div'], t_div=S.PARAMS['t_div'], t_influx=influx,
                 exhaust_model='twostate', **params)
    if arch == 'follicle':
        m.seed_follicles([(S.L // 2, S.L // 2)], S.R_BIG)
    elif arch == 'multi':
        r, d, c = S.R_BIG // 2, 30, S.L // 2
        m.seed_follicles([(c - d, c - d), (c - d, c + d), (c + d, c - d), (c + d, c + d)], r)
    else:
        m.seed_dispersed(S.N_B, occupancy=0.95)
    m.seed_tcells(S.N_T if N_T_OVERRIDE is None else N_T_OVERRIDE)
    return m


def run(args):
    arch, sched_name, rep_label, params, influx, seed = args
    on, off = SCH.by_name(sched_name)
    m = build(arch, seed, params, influx)
    n0 = m.nB
    fn = SCH.schedule_fn(on, off)
    per_day = int(1440 / S.DT)
    snap, traj = {}, []
    day_cleared = None
    integ = 0.0
    for step in range(int(DAYS * 1440 / S.DT)):
        m.step(fn(m))
        integ += m.history[-1]['nB'] * S.DT / 1440.0
        if step % per_day == per_day - 1:
            day = (step + 1) // per_day
            h = m.history[-1]
            nT = int(m.T.sum())
            traj.append((day, h['nB'], nT, round(h['meanE'], 4), round(h['engFrac'], 4),
                         round(float(m.Ed[m.T].mean()) if nT else 0.0, 4)))
            if day in (7, 14, 21, 28, 35, 42):
                snap[day] = h['nB']
        if m.nB == 0:
            day_cleared = m.t / 1440.0
            break
    if day_cleared is not None:
        integ += 0.0
        for d in (7, 14, 21, 28, 35, 42):
            snap.setdefault(d, 0)
    hist = m.history
    h28 = hist[:min(len(hist), int(28 * 1440 / S.DT))]
    nT = int(m.T.sum())
    return dict(arch=arch, sched=sched_name, rep=rep_label, influx=influx, seed=seed, n0=n0,
                nB28=snap.get(28, 0), nB42=snap.get(42, 0), snap=snap,
                cleared=bool(m.nB == 0), day_cleared=day_cleared,
                integ_burden=float(integ), kills=int(m.kills), nT_end=nT,
                duty=SCH.duty_cycle(on, off),
                drug_days=float(sum(r['drug'] for r in hist) * S.DT / 1440.0),
                meanE_end=float(m.E[m.T].mean()) if nT else 0.0,
                Ed_end=float(m.Ed[m.T].mean()) if nT else 0.0,
                Er_end=float(m.Er[m.T].mean()) if nT else 0.0,
                f_end=float(1.0 - m.E[m.T].mean()) if nT else 0.0,
                contact_per_T=float(m.C[m.T].mean() / 1440.0) if nT else 0.0,
                meanE28=float(np.mean([r['meanE'] for r in h28])),
                engFrac28=float(np.mean([r['engFrac'] for r in h28])),
                traj=traj)


def load_combos():
    """Every (exhaustion representative x influx) pair admissible FOR THAT REPRESENTATIVE, not
    only the intersection. Using the intersection alone would leave zero influx as the sole
    regime, and zero influx has its own unrealism: with no replenishment the T-cell pool falls
    from 200 to about 25 by day 28 through background death. The union keeps the one regime that
    is admissible on the external functional criterion AND holds the pool roughly stable
    (rec_low at 2.5e-5, pool 160 at day 28), which is the most physiologically sensible member.
    Admissibility was decided in exp_influx_check.py against Philipp's patient data alone."""
    cal = json.load(open(f'{OUT}/calib_philipp.json'))
    verd = json.load(open(f'{OUT}/expC_influx_verdict.json'))
    reps = {k: v for k, v in cal['representatives'].items() if k != 'best_fit'}
    combos = []
    for key, ok in sorted(verd['verdict'].items()):
        if not ok:
            continue
        label, infl = key.split('|')
        if label not in reps:
            continue
        ch = reps[label]
        p = dict(k_exh=ch['k_per_min'], frac_durable=ch['rho'],
                 recover_tau_r=ch['tau_r_min'], c50_exh=ch['c50_min'], n_exh=4.0)
        combos.append((label, p, float(infl)))
    return combos



def logrank(t1, e1, t2, e2):
    """Two-sample log-rank test. t = time, e = 1 if the event (clearance) was observed, 0 if the
    run was censored at the horizon. Implemented here because the pre-registered floor-effect
    contingency in PHASE2_PREREG.md section 5 specifies log-rank, and neither scipy nor lifelines
    provides it in this environment."""
    from scipy import stats as _st
    t1, e1, t2, e2 = map(np.asarray, (t1, e1, t2, e2))
    times = np.unique(np.concatenate([t1[e1 == 1], t2[e2 == 1]]))
    if not times.size:
        return 1.0
    O1 = E1 = V = 0.0
    for t in times:
        n1, n2 = int((t1 >= t).sum()), int((t2 >= t).sum())
        n = n1 + n2
        d1, d2 = int(((t1 == t) & (e1 == 1)).sum()), int(((t2 == t) & (e2 == 1)).sum())
        d = d1 + d2
        if n < 2 or d == 0:
            continue
        O1 += d1
        E1 += d * n1 / n
        V += d * (n1 / n) * (1 - n1 / n) * (n - d) / (n - 1)
    if V <= 0:
        return 1.0
    z2 = (O1 - E1) ** 2 / V
    return float(_st.chi2.sf(z2, 1))


def clearance(r):
    """Time to clearance in days, with runs that never clear censored at the 42-day horizon."""
    return (r['day_cleared'], 1) if r['cleared'] else (float(DAYS), 0)


def analyse(res, arch, combos, out_name):
    from scipy import stats
    print(f'\n{"="*78}\nPOSITIVE CONTROL, architecture = {arch}\n{"="*78}')
    passes, table = [], []
    for label, _p, infl in combos:
        sub = [r for r in res if r['rep'] == label and r['influx'] == infl and r['arch'] == arch]
        if not sub:
            continue
        cont = {r['seed']: r for r in sub if r['sched'] == SCH.CONTINUOUS}
        cleared_frac = np.mean([cont[s]['cleared'] for s in cont])
        endpoint = 'nB42' if cleared_frac <= 0.5 else 'time_to_clearance'
        print(f'\n--- [{label}]  influx {infl:.1e}   '
              f'continuous arm clears {cleared_frac*100:.0f}% of seeds -> '
              f'primary endpoint = {endpoint} ---')
        base = np.array([cont[s]['nB42'] for s in sorted(cont)], float)
        print(f'    continuous: median nB42 = {np.median(base):.0f}  '
              f'integrated burden = {np.median([cont[s]["integ_burden"] for s in sorted(cont)]):.0f}  '
              f'end function f = {np.median([cont[s]["f_end"] for s in sorted(cont)]):.3f}')
        print(f'    {"schedule":12s} {"duty":>5s} {"med nB42":>9s} {"vs cont":>9s} '
              f'{"p(1-sided)":>10s} {"integ":>9s} {"f_end":>6s} {"Ed_end":>6s} verdict')
        for name, on, off in SCH.ALL:
            arm = {r['seed']: r for r in sub if r['sched'] == name}
            if len(arm) < 3:
                continue
            v = np.array([arm[s]['nB42'] for s in sorted(arm)], float)
            b = np.array([cont[s]['nB42'] for s in sorted(arm)], float)
            d = v - b
            if name == SCH.CONTINUOUS:
                p = np.nan; red = 0.0
            elif endpoint == 'nB42':
                p = (stats.wilcoxon(d, alternative='less').pvalue
                     if np.any(d != 0) else 1.0)
                red = -np.median(d) / max(np.median(b), 1e-9)
            else:
                # Pre-registered floor-effect contingency (PHASE2_PREREG.md section 5): the
                # continuous arm clears more than half its seeds, so day-42 burden has no dynamic
                # range. Switch to time to clearance, censored at the horizon, compared by
                # log-rank, with the same 10% minimum effect applied to the median.
                ta, ea = zip(*[clearance(arm[s]) for s in sorted(arm)])
                tb, eb = zip(*[clearance(cont[s]) for s in sorted(arm)])
                p = logrank(ta, ea, tb, eb)
                mb = np.median(tb)
                red = (mb - np.median(ta)) / max(mb, 1e-9)
            ok = (name != SCH.CONTINUOUS and p < 0.05 and red >= MIN_EFFECT)
            if ok:
                passes.append((label, infl, name, float(p), float(red)))
            print(f'    {name:12s} {SCH.duty_cycle(on,off):5.2f} {np.median(v):9.0f} '
                  f'{-np.median(d):+9.0f} {p:10.4f} '
                  f'{np.median([arm[s]["integ_burden"] for s in sorted(arm)]):9.0f} '
                  f'{np.median([arm[s]["f_end"] for s in sorted(arm)]):6.3f} '
                  f'{np.median([arm[s]["Ed_end"] for s in sorted(arm)]):6.3f} '
                  f'{"PASS" if ok else ""}')
            table.append(dict(rep=label, influx=infl, sched=name, duty=SCH.duty_cycle(on, off),
                              endpoint=endpoint,
                              med_nB42=float(np.median(v)), med_diff=float(np.median(d)),
                              p=float(p) if p == p else None, reduction=float(red),
                              med_integ=float(np.median([arm[s]['integ_burden'] for s in sorted(arm)])),
                              f_end=float(np.median([arm[s]['f_end'] for s in sorted(arm)])),
                              passes=bool(ok)))
    print(f'\n{"="*78}')
    if passes:
        print('POSITIVE CONTROL PASSED in the following combinations:')
        for label, infl, name, p, red in passes:
            on, off = SCH.by_name(name)
            short = off <= 3
            print(f'    [{label}] influx {infl:.1e}  {name} '
                  f'({on} on / {off} off)  p={p:.4f}  reduction={red*100:.0f}%'
                  f'{"" if short else "   <-- PC-2 FLAG: not a short interval"}')
    else:
        print('POSITIVE CONTROL FAILED - MODEL NOT VALIDATED FOR SCHEDULING QUESTION')
        print('Per PHASE2_PREREG.md section 9 rule 1, the architecture comparison is NOT run.')
    print(f'{"="*78}')
    json.dump(dict(table=table, passes=passes, min_effect=MIN_EFFECT),
              open(f'{OUT}/{out_name}', 'w'), indent=1)
    return passes


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1].upper() == 'A':
        # Analysis only, from stored raw results. The simulation process imports this module at
        # launch, so a later fix to analyse() is not picked up by a run already in flight; this
        # recomputes the verdict from results/expE3_poscontrol.json with the current code.
        res = json.load(open(f'{OUT}/expE3_poscontrol.json'))
        analyse(res, ARCH, load_combos(), 'expE3_poscontrol_verdict.json')
        sys.exit(0)
    combos = load_combos()
    jobs = [(ARCH, name, label, p, infl, s)
            for label, p, infl in combos
            for name, _on, _off in SCH.ALL
            for s in SEEDS]
    print(f'Phase E positive control: {len(jobs)} runs '
          f'({len(combos)} parameter combinations x {len(SCH.ALL)} schedules x {len(SEEDS)} seeds)',
          flush=True)
    t0 = time.time()
    res = []
    with Pool(4) as pool:
        for i, r in enumerate(pool.imap_unordered(run, jobs, chunksize=1), 1):
            res.append(r)
            if i % 20 == 0 or i == len(jobs):
                el = time.time() - t0
                print(f'  {i}/{len(jobs)}  {el/60:.1f} min, eta {(el/i*(len(jobs)-i))/60:.0f} min',
                      flush=True)
                json.dump(res, open(f'{OUT}/expE3_poscontrol_partial.json', 'w'))
    json.dump(res, open(f'{OUT}/expE3_poscontrol.json', 'w'))
    analyse(res, ARCH, combos, 'expE3_poscontrol_verdict.json')
