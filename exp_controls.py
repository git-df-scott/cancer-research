"""
Experiment L2 (PRE-REGISTERED before running): controls that could kill the L1 claim.

L1 asks whether tissue architecture changes the value of a treatment-free interval (TFI) for a
T-cell engager. Before any such result can be believed, three alternative explanations have to
be excluded. Each is a way the L1 headline could be an artefact rather than a finding.

C1. IS IT JUST DOSE?
    A TFI arm delivers less total drug-time than continuous. If continuous wins in the follicle
    arm, the trivial explanation is "more drug kills more", not "exhaustion did not bind".
    Control: dose-matched continuous arms. In this model the kill hazard is linear in `drug`, so
    an arm held at drug = (28-k)/28 continuously delivers the same integrated drug as the TFI-k
    arm. If TFI-k and its dose-matched continuous twin perform the SAME, the schedule is doing
    nothing beyond dose and there is no timing effect to report.
    (This is the same control the parent solid-tumour project used; it is what separated a
    genuine timing effect from a dose effect there.)

C2. IS EXHAUSTION ACTUALLY THE MECHANISM?
    L1's story is: dense architecture -> low engaged fraction -> slow exhaustion -> TFIs not worth
    much. If that story is right, then DISABLING exhaustion entirely should abolish the
    architecture-dependence of the TFI benefit. If the TFI ordering is unchanged with exhaustion
    switched off, the mechanism is not exhaustion, and L1's interpretation is wrong even if its
    numbers are right.
    Control: exhaust_tonic = 0, everything else identical.

C3. DOES THE MEAN HIDE THE DISTRIBUTION?
    L1 reports MEAN exhaustion across all T cells. In a follicle the T cells at the invasion front
    are in permanent contact while the reservoir behind them is idle, so the mean can be low while
    the cells that actually matter are fully exhausted. If so, a stalled front could coexist with a
    low mean, and "exhaustion does not bind in follicles" would be false.
    Control: record the exhaustion distribution split by engagement state (engaged vs not) and by
    distance from the follicle centre, under continuous dosing.

PRE-REGISTERED PREDICTIONS
--------------------------
D1 (C1) In `dispersed`, at least one TFI arm beats its OWN dose-matched continuous twin on day-42
   burden (median paired difference < 0, Wilcoxon p < 0.05). That is a genuine timing effect.
D2 (C1) In `follicle`, no TFI arm beats its dose-matched twin by that rule. Whatever difference
   exists between TFI and full continuous dosing in the follicle arm is explained by dose alone.
D3 (C2) With exhaustion disabled, the TFI-versus-continuous difference becomes statistically
   indistinguishable between `follicle` and `dispersed`: the architecture x schedule interaction
   that L1 reports disappears. Quantitatively, the gap between the best-TFI benefit in dispersed
   and in follicle shrinks by more than half relative to L1.
D4 (C3) Under continuous dosing in `follicle`, mean exhaustion among ENGAGED T cells is more than
   3x the mean among unengaged T cells, and the engaged-cell mean exceeds 0.15 by day 28. If this
   holds, the honest statement is not "exhaustion does not accumulate in follicles" but "exhaustion
   is confined to a small front sub-population that is continuously replaced from a fresh
   reservoir" - a different and more interesting mechanism, and the write-up must say so.

FALSIFICATION SUMMARY. If D1/D2 fail (TFI beats dose-matched control in the follicle too), the
architecture claim collapses into a dose claim. If D3 fails (removing exhaustion changes nothing),
the mechanism is misattributed. If D4 holds strongly, the mechanism is REPLACEMENT of exhausted
front cells rather than ABSENCE of exhaustion, and the claim must be restated. All three readings
are written down before running.
"""
import json, os, sys, time
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lymphoid import Lymphoid, nbr_sum
import exp_schedule as S

OUT = S.OUT
TFIS_C = [2, 7, 14]


def run(args):
    arch, mode, k, seed = args
    """mode: 'tfi' | 'dosematch' | 'tfi_noexh' | 'cont_noexh' | 'profile'"""
    params = dict(S.PARAMS)
    if mode in ('tfi_noexh', 'cont_noexh'):
        params['exhaust_tonic'] = 0.0
    T = Lymphoid(L=S.L, seed=seed, **params)
    if arch == 'follicle':
        T.seed_follicles([(S.L // 2, S.L // 2)], S.R_BIG)
    elif arch == 'multi':
        r, d, c = S.R_BIG // 2, 30, S.L // 2
        T.seed_follicles([(c - d, c - d), (c - d, c + d), (c + d, c - d), (c + d, c + d)], r)
    else:
        T.seed_dispersed(S.N_B, occupancy=0.95)
    T.seed_tcells(S.N_T)
    n0 = T.nB

    if mode in ('tfi', 'tfi_noexh'):
        sched = S.schedule_fn(k)
    elif mode in ('dosematch',):
        lvl = (S.CYCLE - k) / S.CYCLE           # same integrated drug as TFI-k
        sched = lambda s: lvl
    else:                                        # cont_noexh, profile
        sched = lambda s: 1.0

    prof = []
    snap = {}
    for step in range(S.DAYS * 1440):
        T.step(sched(T))
        if step % 1440 == 1439:
            day = (step + 1) // 1440
            if day in (7, 14, 28, 42):
                snap[day] = T.history[-1]['nB']
            if mode == 'profile' and day in (7, 14, 21, 28):
                eng = (nbr_sum(T.B) > 0) & T.T
                une = T.T & ~eng
                c = S.L // 2
                yy, xx = np.mgrid[:S.L, :S.L]
                rr = np.hypot(yy - c, xx - c)
                bins = []
                for lo, hi in ((0, 14), (14, 28), (28, 42), (42, 60)):
                    m = T.T & (rr >= lo) & (rr < hi)
                    bins.append(dict(lo=lo, hi=hi, n=int(m.sum()),
                                     meanE=float(T.E[m].mean()) if m.any() else None))
                prof.append(dict(day=day,
                                 n_eng=int(eng.sum()), n_une=int(une.sum()),
                                 E_eng=float(T.E[eng].mean()) if eng.any() else None,
                                 E_une=float(T.E[une].mean()) if une.any() else None,
                                 E_max=float(T.E[T.T].max()) if T.T.any() else None,
                                 radial=bins))
        if T.nB == 0:
            break
    h = T.history
    return dict(arch=arch, mode=mode, k=k, seed=seed, n0=n0,
                nB28=snap.get(28, 0), nB42=snap.get(42, 0), cleared=bool(T.nB == 0),
                kills=T.kills, nT_end=int(T.T.sum()),
                meanE_end=float(h[-1]['meanE']),
                engFrac28=float(np.mean([r['engFrac'] for r in h[:min(len(h), 28 * 1440)]])),
                total_drug=float(sum(r['drug'] for r in h) / 1440.0),   # drug-days delivered
                profile=prof)


if __name__ == '__main__':
    t0 = time.time()
    jobs = []
    # C1 dose-matched twins, and the TFI arms re-run alongside them for a clean paired contrast
    for a in ('follicle', 'dispersed'):
        for k in TFIS_C:
            for s in S.SEEDS:
                jobs.append((a, 'dosematch', k, s))
    # C2 exhaustion disabled
    for a in ('follicle', 'dispersed'):
        for s in S.SEEDS:
            jobs.append((a, 'cont_noexh', 0, s))
        for k in TFIS_C:
            for s in S.SEEDS:
                jobs.append((a, 'tfi_noexh', k, s))
    # C3 exhaustion distribution, continuous, fewer seeds (measurement not comparison)
    for a in ('follicle', 'dispersed'):
        for s in S.SEEDS[:6]:
            jobs.append((a, 'profile', 0, s))
    print(f'exp L2: {len(jobs)} runs', flush=True)
    with Pool(6) as p:
        res = p.map(run, jobs, chunksize=1)
    json.dump(res, open(f'{OUT}/expL2.json', 'w'))
    msg = f'expL2: {len(jobs)} runs in {time.time()-t0:.0f}s'
    open(f'{OUT}/expL2_log.txt', 'w').write(msg + '\n')
    print(msg)
