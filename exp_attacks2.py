"""
Phase H: causal tests. Exhaustion knockout, dose matching, exhaustion distribution, and the
geometry continuum. Pre-registered in PHASE2_PREREG.md sections 6 (H5) and 7.2-7.4.

These run only if an architecture x schedule interaction survives Phase G. Each is designed to
destroy the proposed causal chain, not to support it.

  K  exhaustion knockout      k_exh = 0. Architecture and delivered exposure unchanged. If the
                              interaction survives essentially unchanged without any exhaustion
                              at all, the chain architecture -> contact -> exhaustion -> schedule
                              is FALSE and the real cause must be found and reported.

  D  dose matching            kill hazard is linear in drug, so a schedule with duty cycle d is
                              matched by a continuous arm at level d. If an interval does not beat
                              its own dose-matched twin, there is no timing effect, only dose.
                              Philipp already ran this experimentally: her TFI arm at day 21 has
                              had the same 14 days of cumulative exposure as the continuous arm at
                              day 14, and is at 58.7% specific lysis against 34.9%.

  P  exhaustion distribution  exhaustion among engaged vs unengaged T cells and by radius. If the
                              population mean is low while the invasion front is fully exhausted
                              and continuously replaced from a fresh reservoir, that is a DIFFERENT
                              mechanism and the claim must be restated, not preserved.

  C  geometry continuum       packing varied continuously at fixed total malignant burden, rather
                              than three named categories. N follicles of radius sqrt(N_B/(N pi))
                              tile the same total area, so total malignant perimeter grows as
                              sqrt(N) while burden is held constant. Tests H5: whether schedule
                              behaviour is predicted by a dimensionless quantity - engaged
                              fraction, contact dwell time per unit treatment time, the ratio of
                              exhaustion to recovery timescale, or surface-to-volume - rather than
                              by a fitted geometry label.

usage:  python exp_attacks2.py K|D|P|C
"""
import json, os, sys, time
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lymphoid import Lymphoid, nbr_sum
import exp_schedule as S
import schedules as SCH
from exp_poscontrol import build
from exp_arch2 import ARCHS, SEEDS, passing_combos

OUT = S.OUT
DAYS = 42
SUBSET = ['A_cont', 'B_6on1off', 'B_MO_FR', 'A_tfi7']
# Geometry continuum: follicle counts tiling the same total malignant area.
N_FOLL = [1, 4, 9, 16, 36, 64]


def seed_continuum(m, n_foll):
    """n_foll follicles of equal radius, same total malignant area as the single-follicle case,
    laid out on a square grid inside the lattice."""
    total = np.pi * S.R_BIG ** 2
    r = int(round(np.sqrt(total / (n_foll * np.pi))))
    side = int(round(np.sqrt(n_foll)))
    step = S.L / float(side + 1)
    centres = [(int(round(step * (i + 1))), int(round(step * (j + 1))))
               for i in range(side) for j in range(side)][:n_foll]
    m.seed_follicles(centres, max(r, 1))
    return r, len(centres)


def measure(m, day):
    eng = (nbr_sum(m.B) > 0) & m.T
    une = m.T & ~eng
    c = S.L // 2
    yy, xx = np.mgrid[:S.L, :S.L]
    rr = np.hypot(yy - c, xx - c)
    bins = []
    for lo, hi in ((0, 14), (14, 28), (28, 42), (42, 60)):
        msk = m.T & (rr >= lo) & (rr < hi)
        bins.append(dict(lo=lo, hi=hi, n=int(msk.sum()),
                         E=(float(m.E[msk].mean()) if msk.any() else None),
                         Ed=(float(m.Ed[msk].mean()) if msk.any() else None)))
    Ev = m.E[m.T]
    return dict(day=day, n_eng=int(eng.sum()), n_une=int(une.sum()),
                E_eng=(float(m.E[eng].mean()) if eng.any() else None),
                E_une=(float(m.E[une].mean()) if une.any() else None),
                Ed_eng=(float(m.Ed[eng].mean()) if eng.any() else None),
                Ed_une=(float(m.Ed[une].mean()) if une.any() else None),
                E_p50=(float(np.percentile(Ev, 50)) if Ev.size else None),
                E_p90=(float(np.percentile(Ev, 90)) if Ev.size else None),
                E_max=(float(Ev.max()) if Ev.size else None), radial=bins)


def run(args):
    mode, arch, sched_name, rep, params, influx, seed, extra = args
    p = dict(params)
    if mode == 'K':
        p['k_exh'] = 0.0                      # no exhaustion at all
    m = Lymphoid(L=S.L, seed=seed, dt=S.DT, p_kill=S.PARAMS['p_kill'],
                 p_div=S.PARAMS['p_div'], t_div=S.PARAMS['t_div'], t_influx=influx,
                 exhaust_model='twostate', **p)
    n_foll = None
    if mode == 'C':
        r, n_foll = seed_continuum(m, extra)
    elif arch == 'follicle':
        m.seed_follicles([(S.L // 2, S.L // 2)], S.R_BIG)
    elif arch == 'multi':
        r, d, c = S.R_BIG // 2, 30, S.L // 2
        m.seed_follicles([(c - d, c - d), (c - d, c + d), (c + d, c - d), (c + d, c + d)], r)
    else:
        m.seed_dispersed(S.N_B, occupancy=0.95)
    m.seed_tcells(S.N_T)
    n0 = m.nB
    perim = int(((nbr_sum(m.B) < 8) & m.B).sum())     # malignant cells with a non-malignant neighbour

    on, off = SCH.by_name(sched_name)
    if mode == 'D':
        lvl = SCH.duty_cycle(on, off)
        fn = lambda s: lvl                    # dose-matched continuous twin
    else:
        fn = SCH.schedule_fn(on, off)

    per_day = int(1440 / S.DT)
    snap, prof, traj = {}, [], []
    integ = 0.0
    for step in range(int(DAYS * 1440 / S.DT)):
        m.step(fn(m))
        integ += m.history[-1]['nB'] * S.DT / 1440.0
        if step % per_day == per_day - 1:
            day = (step + 1) // per_day
            h = m.history[-1]; nT = int(m.T.sum())
            traj.append((day, h['nB'], nT, round(h['meanE'], 4), round(h['engFrac'], 4)))
            if day in (7, 14, 21, 28, 35, 42):
                snap[day] = h['nB']
            if mode == 'P' and day in (7, 14, 28, 42):
                prof.append(measure(m, day))
        if m.nB == 0:
            for d in (7, 14, 21, 28, 35, 42):
                snap.setdefault(d, 0)
            break
    nT = int(m.T.sum())
    h28 = m.history[:min(len(m.history), int(28 * 1440 / S.DT))]
    return dict(mode=mode, arch=arch, sched=sched_name, rep=rep, influx=influx, seed=seed,
                n_foll=n_foll, n0=n0, perimeter0=perim,
                nB28=snap.get(28, 0), nB42=snap.get(42, 0), snap=snap,
                cleared=bool(m.nB == 0), integ_burden=float(integ), kills=int(m.kills),
                nT_end=nT, duty=SCH.duty_cycle(on, off),
                drug_days=float(sum(r['drug'] for r in m.history) * S.DT / 1440.0),
                meanE_end=float(m.E[m.T].mean()) if nT else 0.0,
                Ed_end=float(m.Ed[m.T].mean()) if nT else 0.0,
                f_end=float(1.0 - m.E[m.T].mean()) if nT else 0.0,
                contact_per_T=float(m.C[m.T].mean() / 1440.0) if nT else 0.0,
                engFrac28=float(np.mean([r['engFrac'] for r in h28])),
                meanE28=float(np.mean([r['meanE'] for r in h28])),
                profile=prof, traj=traj)


if __name__ == '__main__':
    mode = sys.argv[1].upper()
    combos = passing_combos()
    if mode == 'K':
        jobs = [('K', a, n, lab, p, infl, s, None)
                for lab, p, infl in combos for a in ARCHS for n in SUBSET for s in SEEDS]
    elif mode == 'D':
        jobs = [('D', a, n, lab, p, infl, s, None)
                for lab, p, infl in combos for a in ARCHS
                for n in SUBSET if n != SCH.CONTINUOUS for s in SEEDS]
    elif mode == 'P':
        jobs = [('P', a, n, lab, p, infl, s, None)
                for lab, p, infl in combos for a in ARCHS
                for n in ('A_cont', 'B_MO_FR') for s in SEEDS[:5]]
    elif mode == 'C':
        jobs = [('C', 'continuum', n, lab, p, infl, s, nf)
                for lab, p, infl in combos for nf in N_FOLL
                for n in SUBSET for s in SEEDS[:6]]
    else:
        raise SystemExit('mode must be K, D, P or C')
    print(f'Phase H mode {mode}: {len(jobs)} runs', flush=True)
    t0 = time.time(); res = []
    with Pool(4) as pool:
        for i, r in enumerate(pool.imap_unordered(run, jobs, chunksize=1), 1):
            res.append(r)
            if i % 20 == 0 or i == len(jobs):
                el = time.time() - t0
                print(f'  {i}/{len(jobs)}  {el/60:.1f} min, eta '
                      f'{(el/i*(len(jobs)-i))/60:.0f} min', flush=True)
                json.dump(res, open(f'{OUT}/expH_{mode}.json.partial', 'w'))
    json.dump(res, open(f'{OUT}/expH_{mode}.json', 'w'))
    print(f'written {OUT}/expH_{mode}.json')
