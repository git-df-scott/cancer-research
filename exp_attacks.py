"""
Experiment L2: the attack battery. Phases 3, 4 and 5 of the mandate, plus the Phase 1
exhaustion-distribution measurement that L1 does not record.

Every arm here is a control designed to DESTROY the L1 hypothesis, not to support it. Arms and
decision rules were written before any L1 result was examined (see docs/calibration/TRAFFICKING_PREREG.md for the
Phase 5 commitment, which was made before running anything).

  P3  exhaustion knockout      exhaust_tonic = 0. If the architecture x schedule interaction
                               survives without exhaustion, the proposed causal chain
                               (architecture -> engagement -> exhaustion -> scheduling) is FALSE.
  P4  dose-matched control     a continuous arm delivering the same integrated engager exposure as
                               each TFI arm. Kill hazard is linear in `drug`, so TFI-k is matched by
                               a continuous level (28-k)/28. If a TFI arm does not beat its own
                               dose-matched twin, there is no timing effect at all, only dose.
  P5  trafficking knockout     swap_prob 0.5 and 1.0. T cells squeeze past malignant B cells instead
                               of requiring a vacancy. Pre-registered: if the interaction shrinks by
                               more than half at swap_prob = 0.5, the effect is model-dependent
                               (category B) and must NOT be reported as an FL prediction.
  PROF exhaustion distribution  exhaustion among ENGAGED vs UNENGAGED T cells, and by radius, under
                               continuous dosing. Tests whether a low MEAN hides a fully exhausted
                               front that is continuously replaced from a fresh reservoir - which
                               would be a different mechanism and would require restating the claim.
"""
import json, os, sys, time
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lymphoid import Lymphoid, nbr_sum
import exp_schedule as S

OUT = S.OUT
TFIS_A = [0, 2, 7]


def build(arch, seed, **kw):
    p = dict(S.PARAMS); p.update(kw)
    T = Lymphoid(L=S.L, seed=seed, dt=S.DT, **p)
    if arch == 'follicle':
        T.seed_follicles([(S.L // 2, S.L // 2)], S.R_BIG)
    elif arch == 'multi':
        r, d, c = S.R_BIG // 2, 30, S.L // 2
        T.seed_follicles([(c - d, c - d), (c - d, c + d), (c + d, c - d), (c + d, c + d)], r)
    else:
        T.seed_dispersed(S.N_B, occupancy=0.95)
    T.seed_tcells(S.N_T)
    return T


def run(args):
    arch, mode, k, seed = args
    kw = {}
    if mode == 'noexh':
        kw['exhaust_tonic'] = 0.0
    if mode == 'swap05':
        kw['swap_prob'] = 0.5
    if mode == 'swap10':
        kw['swap_prob'] = 1.0
    T = build(arch, seed, **kw)
    n0 = T.nB

    if mode == 'dosematch':
        lvl = (S.CYCLE - k) / S.CYCLE
        sched = lambda s: lvl
    else:
        sched = S.schedule_fn(k)

    nsteps = int(S.DAYS * 1440 / S.DT)
    per_day = int(1440 / S.DT)
    snap, prof, traj = {}, [], []
    for step in range(nsteps):
        T.step(sched(T))
        if step % per_day == per_day - 1:
            day = (step + 1) // per_day
            h = T.history[-1]
            traj.append((day, h['nB'], h['nT'], round(h['meanE'], 4), round(h['engFrac'], 4)))
            if day in (7, 14, 28, 42):
                snap[day] = h['nB']
            if mode == 'profile' and day in (7, 14, 28, 42):
                eng = (nbr_sum(T.B) > 0) & T.T
                une = T.T & ~eng
                c = S.L // 2
                yy, xx = np.mgrid[:S.L, :S.L]
                rr = np.hypot(yy - c, xx - c)
                bins = []
                for lo, hi in ((0, 14), (14, 28), (28, 42), (42, 60)):
                    m = T.T & (rr >= lo) & (rr < hi)
                    bins.append(dict(lo=lo, hi=hi, n=int(m.sum()),
                                     meanE=(float(T.E[m].mean()) if m.any() else None)))
                Ev = T.E[T.T]
                prof.append(dict(day=day, n_eng=int(eng.sum()), n_une=int(une.sum()),
                                 E_eng=(float(T.E[eng].mean()) if eng.any() else None),
                                 E_une=(float(T.E[une].mean()) if une.any() else None),
                                 E_p50=(float(np.percentile(Ev, 50)) if Ev.size else None),
                                 E_p90=(float(np.percentile(Ev, 90)) if Ev.size else None),
                                 E_max=(float(Ev.max()) if Ev.size else None),
                                 radial=bins))
        if T.nB == 0:
            break
    h = T.history
    h28 = h[:min(len(h), int(28 * 1440 / S.DT))]
    return dict(arch=arch, mode=mode, k=k, seed=seed, n0=n0,
                nB7=snap.get(7, 0), nB14=snap.get(14, 0), nB28=snap.get(28, 0), nB42=snap.get(42, 0),
                cleared=bool(T.nB == 0), kills=T.kills, nT_end=int(T.T.sum()),
                meanE_end=float(h[-1]['meanE']),
                meanE28=float(np.mean([r['meanE'] for r in h28])),
                engFrac28=float(np.mean([r['engFrac'] for r in h28])),
                contact_per_T=float(T.cum_contact / max(int(T.T.sum()), 1)),
                total_drug=float(sum(r['drug'] for r in h) * S.DT / 1440.0),
                profile=prof, traj=traj)


if __name__ == '__main__':
    t0 = time.time()
    A = ('follicle', 'dispersed')
    jobs = []
    for a in A:
        for s in S.SEEDS:
            for k in TFIS_A:
                jobs.append((a, 'noexh', k, s))       # P3
                jobs.append((a, 'swap05', k, s))      # P5
            for k in (2, 7):
                jobs.append((a, 'dosematch', k, s))   # P4
            for k in (0, 7):
                jobs.append((a, 'swap10', k, s))      # P5 extreme
        for s in S.SEEDS[:6]:
            jobs.append((a, 'profile', 0, s))         # exhaustion distribution
    print(f'exp L2 attacks: {len(jobs)} runs', flush=True)
    res = []
    with Pool(6) as p:
        for i, r in enumerate(p.imap_unordered(run, jobs, chunksize=1), 1):
            res.append(r)
            if i % 10 == 0 or i == len(jobs):
                el = time.time() - t0
                print(f'  {i}/{len(jobs)}  {el/60:.1f} min elapsed, eta {(el/i*(len(jobs)-i))/60:.0f} min', flush=True)
                json.dump(res, open(f'{OUT}/expL2_partial.json', 'w'))
    json.dump(res, open(f'{OUT}/expL2.json', 'w'))
    msg = f'expL2: {len(jobs)} runs in {time.time()-t0:.0f}s'
    open(f'{OUT}/expL2_log.txt', 'w').write(msg + '\n')
    print(msg)
