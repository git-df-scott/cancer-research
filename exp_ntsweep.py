"""
Phase 3 sweep: locate the initial T-cell number at which the model has dynamic range.

Pre-registered in PHASE3_PREREG.md section 3, committed before this was run.

THIS IS CALIBRATION, NOT A HYPOTHESIS TEST. It simulates the CONTINUOUS arm ONLY. No
treatment-free interval, no schedule comparison and no architecture comparison is computed here,
so the selection is structurally incapable of expressing a preference about schedule ranking.
That is the safeguard that makes changing N_T a dynamic-range fix rather than tuning.

WINDOW CRITERION, fixed before running. N_T is admissible if, in the continuous arm:
  1. median nB42 > 0                      the tumour is not eradicated  (no floor)
  2. cleared fraction <= 0.25
  3. median nB42 < n0                     the drug demonstrably does something
  4. max burden at any time < 70% of lattice capacity                   (no ceiling)
Among qualifying values take the one whose median nB42 is closest to n0/2. If none qualifies,
STOP and report that; do not widen the grid in search of a pass.
"""
import json, os, sys, time
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lymphoid import Lymphoid
import exp_schedule as S
import schedules as SCH
from exp_poscontrol import load_combos

OUT = S.OUT
NT_GRID = [200, 400, 800, 1400, 2000, 2800]
SEEDS = list(range(6))
DAYS = 42
CAP = S.L * S.L
CEILING = 0.70


def run(args):
    rep, params, influx, n_t, seed = args
    m = Lymphoid(L=S.L, seed=seed, dt=S.DT, p_kill=S.PARAMS['p_kill'],
                 p_div=S.PARAMS['p_div'], t_div=S.PARAMS['t_div'], t_influx=influx,
                 exhaust_model='twostate', **params)
    m.seed_dispersed(S.N_B, occupancy=0.95)
    m.seed_tcells(n_t)
    n0 = m.nB
    nT0 = int(m.T.sum())
    per_day = int(1440 / S.DT)
    traj, peak = [], n0
    for step in range(int(DAYS * 1440 / S.DT)):
        m.step(1.0)                                   # CONTINUOUS ONLY
        peak = max(peak, m.history[-1]['nB'])
        if step % per_day == per_day - 1:
            h = m.history[-1]
            traj.append((step // per_day + 1, h['nB'], int(m.T.sum()), round(h['meanE'], 4)))
        if m.nB == 0:
            break
    nT = int(m.T.sum())
    return dict(rep=rep, influx=influx, n_t=n_t, seed=seed, n0=n0, nT0=nT0,
                nB42=(m.history[-1]['nB'] if m.nB else 0), peak=int(peak),
                cleared=bool(m.nB == 0), day_cleared=(m.t / 1440.0 if m.nB == 0 else None),
                kills=int(m.kills), nT_end=nT,
                f_end=float(1.0 - m.E[m.T].mean()) if nT else 0.0, traj=traj)


if __name__ == '__main__':
    combos = load_combos()
    jobs = [(rep, p, infl, n_t, s)
            for rep, p, infl in combos for n_t in NT_GRID for s in SEEDS]
    print(f'Phase 3 N_T sweep: {len(jobs)} runs, continuous arm only', flush=True)
    t0 = time.time(); res = []
    with Pool(4) as pool:
        for i, r in enumerate(pool.imap_unordered(run, jobs, chunksize=1), 1):
            res.append(r)
            if i % 20 == 0 or i == len(jobs):
                el = time.time() - t0
                print(f'  {i}/{len(jobs)}  {el/60:.1f} min, eta {(el/i*(len(jobs)-i))/60:.0f} min',
                      flush=True)
                json.dump(res, open(f'{OUT}/expNT_sweep.json.partial', 'w'))
    json.dump(res, open(f'{OUT}/expNT_sweep.json', 'w'))

    print(f'\n{"="*86}\nWINDOW CRITERION (PHASE3_PREREG.md section 3), continuous arm only\n{"="*86}')
    print(f'  lattice capacity {CAP}, ceiling {CEILING*100:.0f}% = {CAP*CEILING:.0f} cells\n')
    qualify = {}
    for rep, infl in sorted({(r['rep'], r['influx']) for r in res}):
        print(f'--- [{rep}] influx {infl:.1e} ---')
        print(f'  {"N_T":>5s} {"E:T":>7s} {"med nB42":>9s} {"% n0":>6s} {"cleared":>8s} '
              f'{"peak %cap":>10s} {"f_end":>6s}  verdict')
        for n_t in NT_GRID:
            sub = [r for r in res if r['rep'] == rep and r['influx'] == infl and r['n_t'] == n_t]
            if not sub:
                continue
            n0 = np.median([r['n0'] for r in sub])
            med = np.median([r['nB42'] for r in sub])
            clr = float(np.mean([r['cleared'] for r in sub]))
            pk = np.median([r['peak'] for r in sub]) / CAP
            ok = (med > 0) and (clr <= 0.25) and (med < n0) and (pk < CEILING)
            why = [] if ok else (
                (['eradicated'] if med <= 0 else []) +
                (['clears>25%'] if clr > 0.25 else []) +
                (['no effect'] if med >= n0 else []) +
                (['hits ceiling'] if pk >= CEILING else []))
            if ok:
                qualify[(rep, infl, n_t)] = abs(med - n0 / 2)
            print(f'  {n_t:5d} {"1:"+f"{S.N_B/n_t:.0f}":>7s} {med:9.0f} {100*med/n0:5.0f}% '
                  f'{100*clr:7.0f}% {100*pk:9.1f}% '
                  f'{np.median([r["f_end"] for r in sub]):6.3f}  '
                  f'{"ADMISSIBLE" if ok else ", ".join(why)}')
        print()
    if not qualify:
        print('NO VALUE OF N_T QUALIFIES. Per PHASE3_PREREG.md section 3, STOP and report.')
        print('The grid is NOT widened in search of a pass.')
    else:
        best_nt = {}
        for (rep, infl, n_t), dist in sorted(qualify.items(), key=lambda kv: kv[1]):
            best_nt.setdefault((rep, infl), n_t)
        print('SELECTED N_T (closest median nB42 to n0/2) per combination:')
        for (rep, infl), n_t in sorted(best_nt.items()):
            print(f'    [{rep}] influx {infl:.1e}  ->  N_T = {n_t}  (E:T 1:{S.N_B/n_t:.0f})')
        counts = {}
        for v in best_nt.values():
            counts[v] = counts.get(v, 0) + 1
        consensus = max(counts, key=counts.get)
        print(f'\n  Most frequently selected across combinations: N_T = {consensus}')
        json.dump(dict(selected={f'{k[0]}|{k[1]:.1e}': v for k, v in best_nt.items()},
                       consensus=int(consensus), grid=NT_GRID, ceiling=CEILING),
                  open(f'{OUT}/expNT_sweep_verdict.json', 'w'), indent=1)
        print(f'  written {OUT}/expNT_sweep_verdict.json')
