"""
Phase 4 window sweep. Pre-registered in PHASE4_PREREG.md section 4, committed before this ran.

THE DECISIVE QUESTION, with both answers pre-committed:
  Does any regime satisfy BOTH tumour control AND binding T-cell exhaustion? In the previous model
  this was false in 24 of 24 cells and the window was empty, which is why three experiments failed.
  If it is still empty here, the systemic compartment does not repair it, the empty-window finding
  is robust, and the project stops.

CONTINUOUS DOSING ONLY. No treatment-free interval or schedule comparison is computed, so the
selection is structurally incapable of expressing a preference about schedule ranking.

Admissible requires BOTH:
  1. external check, unchanged from Phase C: mean pool function <= 0.40 of day-0 by day 14
  2. dynamic range: median nB42 > 0; cleared <= 25%; median nB42 < n0; peak < 70% of capacity
"""
import json, os, sys, time
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lymphoid import Lymphoid
import exp_schedule as S

OUT = S.OUT
NT_GRID = [200, 400, 800, 1400]
INFLUX = [0.0, 2.5e-5, 1.0e-4, 4.0e-4]
SEEDS = list(range(4))
DAYS, CAP, CEILING = 42, S.L * S.L, 0.70


def run(args):
    tag, per_cell, sysp, n_t, infl, seed = args
    m = Lymphoid(L=S.L, seed=seed, dt=S.DT, p_kill=S.PARAMS['p_kill'],
                 p_div=S.PARAMS['p_div'], t_div=S.PARAMS['t_div'], t_influx=infl,
                 exhaust_model='twostate', systemic=True, **per_cell, **sysp)
    m.seed_dispersed(S.N_B, occupancy=0.95)
    m.seed_tcells(n_t)
    n0 = m.nB
    per_day = int(1440 / S.DT)
    peak, snap = n0, {}
    for step in range(int(DAYS * 1440 / S.DT)):
        m.step(1.0)                                    # CONTINUOUS ONLY
        peak = max(peak, m.history[-1]['nB'])
        if step % per_day == per_day - 1:
            day = (step + 1) // per_day
            nT = int(m.T.sum())
            if day in (7, 14, 28, 42):
                snap[day] = dict(nB=m.history[-1]['nB'], nT=nT,
                                 f=float(1.0 - m.E[m.T].mean()) if nT else 0.0,
                                 E_sys=float(m.history[-1]['E_sys']))
        if m.nB == 0:
            for d in (7, 14, 28, 42):
                snap.setdefault(d, dict(nB=0, nT=int(m.T.sum()), f=0.0, E_sys=0.0))
            break
    nT = int(m.T.sum())
    return dict(tag=tag, n_t=n_t, influx=infl, seed=seed, n0=n0, peak=int(peak),
                nB42=(m.history[-1]['nB'] if m.nB else 0), cleared=bool(m.nB == 0),
                nT_end=nT, f_end=float(1.0 - m.E[m.T].mean()) if nT else 0.0,
                E_sys_end=float(m.history[-1]['E_sys']), snap=snap, kills=int(m.kills))


if __name__ == '__main__':
    cal = json.load(open(f'{OUT}/calib_philipp.json'))
    ch = cal['representatives']['rec_low']
    per_cell = dict(k_exh=ch['k_per_min'], frac_durable=ch['rho'],
                    recover_tau_r=ch['tau_r_min'], c50_exh=ch['c50_min'], n_exh=4.0)
    sysc = json.load(open(f'{OUT}/calib_systemic.json'))
    variants = {k: dict(k_sys=v['k_sys_per_min'], tau_sys=v['tau_sys_min'],
                        rho_sys=v['rho_sys']) for k, v in sysc.items()}
    jobs = [(tag, per_cell, sp, n_t, infl, s)
            for tag, sp in variants.items() for n_t in NT_GRID for infl in INFLUX for s in SEEDS]
    print(f'Phase 4 window sweep: {len(jobs)} runs, continuous arm only', flush=True)
    t0 = time.time(); res = []
    with Pool(4) as pool:
        for i, r in enumerate(pool.imap_unordered(run, jobs, chunksize=1), 1):
            res.append(r)
            if i % 20 == 0 or i == len(jobs):
                el = time.time() - t0
                print(f'  {i}/{len(jobs)}  {el/60:.1f} min, eta {(el/i*(len(jobs)-i))/60:.0f} min',
                      flush=True)
                json.dump(res, open(f'{OUT}/expW4.json.partial', 'w'))
    json.dump(res, open(f'{OUT}/expW4.json', 'w'))

    print(f'\n{"="*96}\nTHE DECISIVE QUESTION: does any regime give BOTH control AND binding '
          f'exhaustion?\n{"="*96}')
    print(f'  {"variant":10s} {"N_T":>5s} {"influx":>8s} | {"nB42/n0":>8s} {"f(d14)":>7s} '
          f'{"f_end":>7s} {"E_sys":>7s} {"nT_end":>7s} {"peak%":>6s} | ext? range? BOTH?')
    ok_any, adm = False, []
    for tag in sorted(variants):
        for n_t in NT_GRID:
            for infl in INFLUX:
                sub = [r for r in res if r['tag'] == tag and r['n_t'] == n_t and r['influx'] == infl]
                if not sub:
                    continue
                n0 = np.median([r['n0'] for r in sub])
                med = np.median([r['nB42'] for r in sub])
                f14 = np.median([r['snap'][14]['f'] for r in sub])
                fend = np.median([r['f_end'] for r in sub])
                clr = float(np.mean([r['cleared'] for r in sub]))
                pk = np.median([r['peak'] for r in sub]) / CAP
                ext = f14 <= 0.40
                rng = (med > 0) and (clr <= 0.25) and (med < n0) and (pk < CEILING)
                both = (med < n0) and (f14 < 0.5)
                ok_any |= (ext and rng and both)
                if ext and rng:
                    adm.append(dict(tag=tag, n_t=n_t, influx=infl, nB42=float(med),
                                    ratio=float(med/n0), f14=float(f14), f_end=float(fend),
                                    both=bool(both)))
                print(f'  {tag:10s} {n_t:5d} {infl:8.1e} | {med/n0:8.2f} {f14:7.3f} {fend:7.3f} '
                      f'{np.median([r["E_sys_end"] for r in sub]):7.3f} '
                      f'{np.median([r["nT_end"] for r in sub]):7.0f} {100*pk:5.1f}% | '
                      f'{"yes" if ext else "no ":>4s} {"yes" if rng else "no ":>6s} '
                      f'{"** YES **" if (ext and rng and both) else ""}')
        print()
    print('=' * 96)
    if ok_any:
        print('THE WINDOW IS OPEN. Regimes exist with tumour control AND binding exhaustion.')
        print('Admissible regimes carried to the Phase 4 positive control:')
        for a in [x for x in adm if x['both']]:
            print(f"    {a['tag']} N_T={a['n_t']} influx={a['influx']:.1e}  "
                  f"burden {a['ratio']:.2f}x n0, pool function at d14 {a['f14']:.3f}")
    else:
        print('THE WINDOW IS STILL EMPTY.')
        print('Per PHASE4_PREREG.md section 4: STOP. The systemic compartment does not repair it,')
        print('the empty-window finding is robust to this repair, and classification A stands.')
        print('Do NOT widen the grid, adjust k_sys, or attempt a fourth mechanism.')
    print('=' * 96)
    json.dump(dict(window_open=bool(ok_any), admissible=adm),
              open(f'{OUT}/expW4_verdict.json', 'w'), indent=1)
