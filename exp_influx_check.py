"""
Phase C: external population-level check on T-cell influx.

THE PROBLEM THIS TESTS
----------------------
Experiment L1's diagnosis was that fresh T-cell influx dilutes the exhaustion pool: the T-cell
population grows from 200 to ~860 over a run, so most effectors present at day 28 are recent
arrivals carrying no accrued exhaustion, and the population mean is held low regardless of
architecture or schedule. If that dilution is an artefact of the influx rate, the whole model is
structurally unable to show exhaustion-driven scheduling effects, and no amount of recalibrating
the exhaustion clock will fix it.

THE EXTERNAL ANCHOR
-------------------
Philipp et al. Figure 1B, which is patient data and is independent of the in vitro curve used to
calibrate the clock. In relapsed/refractory B-cell precursor ALL patients on continuous
blinatumomab infusion, ex vivo specific lysis by PERIPHERAL T cells fell

    day 0 = 73.1%   ->   day 14 = 17.4%   (P = .0176)

i.e. to 0.238 of baseline, and recovered to 48.5% after cessation (not significantly different
from baseline). That is a POPULATION-level collapse, measured in a compartment where T-cell
redistribution, recruitment and blinatumomab-driven expansion are all occurring. So in the real
system, influx does NOT prevent the effector pool from losing most of its function.

PRE-REGISTERED ADMISSIBILITY CRITERION (fixed before running)
-------------------------------------------------------------
An influx regime is ADMISSIBLE only if, under continuous dosing in the dispersed (leukaemia-like)
architecture that matches the patient setting, the mean functional capacity of the whole T-cell
pool has fallen to <= 0.40 of its day-0 value by day 14.

The threshold is 0.40 rather than the measured 0.238 deliberately, to be generous to the model:
the lattice compartment is not peripheral blood, and the model's engaged fraction is below 1 even
in the dispersed arm, so some shortfall is expected for reasons unrelated to influx. A regime that
cannot even reach 0.40 is not marginally wrong, it is qualitatively wrong.

Regimes failing this are EXCLUDED before the architecture experiment. This criterion is about
reproducing an external measurement. It says nothing about, and is not permitted to be adjusted
by, any treatment-schedule ranking.
"""
import json, os, sys, time
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lymphoid import Lymphoid
import exp_schedule as S

OUT = S.OUT
INFLUX = [0.0, 2.5e-5, 1.0e-4, 4.0e-4, 1.6e-3]      # 1.0e-4 is the value used in experiment L1
SEEDS = list(range(4))
DAYS = 28
ADMISSIBLE_MAX_F14 = 0.40


def run(args):
    label, params, infl, seed = args
    m = Lymphoid(L=S.L, seed=seed, dt=S.DT, p_kill=S.PARAMS['p_kill'],
                 p_div=S.PARAMS['p_div'], t_div=S.PARAMS['t_div'], t_influx=infl,
                 exhaust_model='twostate', **params)
    m.seed_dispersed(S.N_B, occupancy=0.95)
    m.seed_tcells(S.N_T)
    per_day = int(1440 / S.DT)
    snap = {}
    for step in range(int(DAYS * 1440 / S.DT)):
        m.step(1.0)
        if step % per_day == per_day - 1:
            day = (step + 1) // per_day
            if day in (7, 14, 28):
                nT = int(m.T.sum())
                snap[day] = dict(f=float(1.0 - m.E[m.T].mean()) if nT else 0.0,
                                 nT=nT, nB=m.nB,
                                 Ed=float(m.Ed[m.T].mean()) if nT else 0.0,
                                 engFrac=float(m.history[-1]['engFrac']),
                                 frac_naive=float((m.C[m.T] < 60).mean()) if nT else 0.0)
        if m.nB == 0:
            break
    return dict(label=label, influx=infl, seed=seed, snap=snap, n0T=S.N_T)


if __name__ == '__main__':
    cal = json.load(open(f'{OUT}/calib_philipp.json'))
    reps = cal['representatives']
    jobs = []
    for label, ch in reps.items():
        if label == 'best_fit':
            continue                      # best_fit duplicates rec_low's region of the class
        p = dict(k_exh=ch['k_per_min'], frac_durable=ch['rho'],
                 recover_tau_r=ch['tau_r_min'], c50_exh=ch['c50_min'], n_exh=4.0)
        for infl in INFLUX:
            for s in SEEDS:
                jobs.append((label, p, infl, s))
    t0 = time.time()
    print(f'Phase C: {len(jobs)} runs, {DAYS} simulated days, dispersed architecture, '
          f'continuous dosing.\n')
    res = []
    with Pool(4) as pool:
        for i, r in enumerate(pool.imap_unordered(run, jobs, chunksize=1), 1):
            res.append(r)
            if i % 10 == 0 or i == len(jobs):
                el = time.time() - t0
                print(f'  {i}/{len(jobs)}  {el/60:.1f} min, eta {(el/i*(len(jobs)-i))/60:.0f} min',
                      flush=True)
                json.dump(res, open(f'{OUT}/expC_influx_partial.json', 'w'))
    json.dump(res, open(f'{OUT}/expC_influx.json', 'w'))

    print(f'\nExternal anchor: patient peripheral T-cell function fell to 0.238 of baseline by')
    print(f'day 14 of continuous infusion (Philipp Fig 1B). Admissible if model <= '
          f'{ADMISSIBLE_MAX_F14:.2f}.\n')
    verdict = {}
    for label in sorted({r['label'] for r in res}):
        print(f'  [{label}]')
        print('    influx     f(d7)   f(d14)   f(d28)   T pool d28   naive frac d28   verdict')
        for infl in INFLUX:
            sub = [r for r in res if r['label'] == label and r['influx'] == infl]
            g = lambda d, k: float(np.mean([r['snap'][d][k] for r in sub if d in r['snap']]))
            f14 = g(14, 'f')
            ok = f14 <= ADMISSIBLE_MAX_F14
            verdict[(label, infl)] = ok
            print(f'    {infl:.1e}  {g(7,"f"):6.3f}  {f14:6.3f}  {g(28,"f"):6.3f}   '
                  f'{g(28,"nT"):9.0f}   {g(28,"frac_naive"):13.3f}   '
                  f'{"ADMISSIBLE" if ok else "EXCLUDED"}')
        print()
    adm = sorted({infl for (lab, infl), ok in verdict.items() if ok}
                 & set.intersection(*[{i for (l, i), o in verdict.items() if l == lab and o}
                                      for lab in {k[0] for k in verdict}]))
    print(f'Admissible under EVERY exhaustion representative: {adm}')
    json.dump(dict(verdict={f'{k[0]}|{k[1]:.1e}': v for k, v in verdict.items()},
                   admissible_all=adm, threshold=ADMISSIBLE_MAX_F14),
              open(f'{OUT}/expC_influx_verdict.json', 'w'), indent=1)
