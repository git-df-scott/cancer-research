"""
Phase 4 positive control. Dispersed / leukaemia-like architecture ONLY.

Pre-registered in PHASE4_PREREG.md section 5, committed before this ran. The criterion is
UNCHANGED from Phases 2 and 3, and the analysis is imported from exp_poscontrol rather than
rewritten, so the test cannot drift between phases:

  PC-1 (required): at least one non-continuous schedule reduces day-42 burden versus continuous
  at one-sided Wilcoxon p < 0.05 AND by a median of at least 10%.

Regimes come from exp_window4.py, which simulated the CONTINUOUS arm only and therefore could not
express a preference about schedule ranking. Both post-timing variants are carried: the sweep
cannot distinguish them, because tau_sys only acts while the drug is OFF and the sweep never turns
it off. They can differ here, which is precisely why both are run.

THIS IS THE THIRD AND FINAL MODEL. PHASE4_PREREG.md section 5 commits that if PC-1 fails in every
admissible regime, the project is written up as a negative and no fourth mechanism is attempted.
"""
import json, os, sys, time
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lymphoid import Lymphoid
import exp_schedule as S
import schedules as SCH
from exp_poscontrol import analyse, DAYS, SEEDS

OUT = S.OUT


def run(args):
    arch, sched_name, tag, params, influx, seed, n_t = args
    on, off = SCH.by_name(sched_name)
    # systemic is enabled only when the caller supplies calibrated systemic parameters, so this
    # runner also reproduces the Phase 2 and Phase 3 configurations exactly when they are absent.
    # It must NOT default to on with k_sys = 0: the compartment would still accumulate exposure
    # and hand it to arrivals, which changes behaviour whenever the Hill lag c50 > 0.
    m = Lymphoid(L=S.L, seed=seed, dt=S.DT, p_kill=S.PARAMS['p_kill'],
                 p_div=S.PARAMS['p_div'], t_div=S.PARAMS['t_div'], t_influx=influx,
                 exhaust_model='twostate', systemic=('k_sys' in params), **params)
    m.seed_dispersed(S.N_B, occupancy=0.95)
    m.seed_tcells(n_t)
    n0 = m.nB
    fn = SCH.schedule_fn(on, off)
    per_day = int(1440 / S.DT)
    snap, traj, integ, day_cleared = {}, [], 0.0, None
    for step in range(int(DAYS * 1440 / S.DT)):
        m.step(fn(m))
        integ += m.history[-1]['nB'] * S.DT / 1440.0
        if step % per_day == per_day - 1:
            day = (step + 1) // per_day
            h = m.history[-1]; nT = int(m.T.sum())
            traj.append((day, h['nB'], nT, round(h['meanE'], 4), round(h['engFrac'], 4),
                         round(h['E_sys'], 4)))
            if day in (7, 14, 21, 28, 35, 42):
                snap[day] = h['nB']
        if m.nB == 0:
            day_cleared = m.t / 1440.0
            for d in (7, 14, 21, 28, 35, 42):
                snap.setdefault(d, 0)
            break
    nT = int(m.T.sum())
    h28 = m.history[:min(len(m.history), int(28 * 1440 / S.DT))]
    return dict(arch=arch, sched=sched_name, rep=tag, influx=influx, seed=seed, n0=n0, n_t=n_t,
                nB28=snap.get(28, 0), nB42=snap.get(42, 0), snap=snap,
                cleared=bool(m.nB == 0), day_cleared=day_cleared,
                integ_burden=float(integ), kills=int(m.kills), nT_end=nT,
                duty=SCH.duty_cycle(on, off),
                drug_days=float(sum(r['drug'] for r in m.history) * S.DT / 1440.0),
                meanE_end=float(m.E[m.T].mean()) if nT else 0.0,
                Ed_end=float(m.Ed[m.T].mean()) if nT else 0.0,
                Er_end=float(m.Er[m.T].mean()) if nT else 0.0,
                f_end=float(1.0 - m.E[m.T].mean()) if nT else 0.0,
                E_sys_end=float(m.history[-1]['E_sys']),
                contact_per_T=float(m.C[m.T].mean() / 1440.0) if nT else 0.0,
                meanE28=float(np.mean([r['meanE'] for r in h28])),
                engFrac28=float(np.mean([r['engFrac'] for r in h28])), traj=traj)


def combos():
    cal = json.load(open(f'{OUT}/calib_philipp.json'))
    ch = cal['representatives']['rec_low']
    per_cell = dict(k_exh=ch['k_per_min'], frac_durable=ch['rho'],
                    recover_tau_r=ch['tau_r_min'], c50_exh=ch['c50_min'], n_exh=4.0)
    sysc = json.load(open(f'{OUT}/calib_systemic.json'))
    verdict = json.load(open(f'{OUT}/expW4_verdict.json'))
    if not verdict['window_open']:
        print('THE WINDOW IS CLOSED. Per PHASE4_PREREG.md section 4 the positive control is NOT run.')
        sys.exit(1)
    out = []
    for a in verdict['admissible']:
        if not a['both']:
            continue
        v = sysc[a['tag']]
        p = dict(per_cell, k_sys=v['k_sys_per_min'], tau_sys=v['tau_sys_min'],
                 rho_sys=v['rho_sys'])
        label = f"{a['tag']}_NT{a['n_t']}"
        out.append((label, p, float(a['influx']), int(a['n_t'])))
    return out


if __name__ == '__main__':
    cs = combos()
    print(f'Phase 4 positive control: {len(cs)} admissible regimes')
    for label, _p, infl, n_t in cs:
        print(f'    {label}  influx {infl:.1e}  N_T {n_t}')
    jobs = [('dispersed', name, label, p, infl, s, n_t)
            for label, p, infl, n_t in cs
            for name, _o, _f in SCH.ALL for s in SEEDS]
    print(f'{len(jobs)} runs\n', flush=True)
    t0 = time.time(); res = []
    with Pool(4) as pool:
        for i, r in enumerate(pool.imap_unordered(run, jobs, chunksize=1), 1):
            res.append(r)
            if i % 20 == 0 or i == len(jobs):
                el = time.time() - t0
                print(f'  {i}/{len(jobs)}  {el/60:.1f} min, eta '
                      f'{(el/i*(len(jobs)-i))/60:.0f} min', flush=True)
                json.dump(res, open(f'{OUT}/expE4_poscontrol_partial.json', 'w'))
    json.dump(res, open(f'{OUT}/expE4_poscontrol.json', 'w'))
    analyse(res, 'dispersed', [(l, p, i) for l, p, i, _n in cs],
            'expE4_poscontrol_verdict.json')
