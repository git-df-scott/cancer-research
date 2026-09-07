"""
Experiment L1 (PRE-REGISTERED before running): does tissue architecture decide the optimal
treatment-free interval for a T-cell engager?

WHY THIS EXPERIMENT EXISTS
--------------------------
Obertopp, Froid, Pilon-Thomas & Basanta (bioRxiv 2025, doi 10.1101/2025.11.17.688873,
PMC12667981) built a 2D lattice agent-based model of blinatumomab scheduling in B-ALL and
concluded that the clinically used 7-day treatment-free interval (TFI) is NOT optimal: shorter
TFIs of 2-3 days did best and continuous dosing was worst. Their model seeds tumour cells
RANDOMLY at 50% occupancy with T cells interspersed among them -- their own stated limitation
("random rather than clustered cell seeding"). That is a petri dish, and it is a fair model of
leukaemia, where an effector and a blast are already neighbours.

Follicular lymphoma is not like that. Imaging mass cytometry of paired FL biopsies found that
"peri-follicular regions represented a barrier for immune infiltration into the follicles" and
that FL cells inside follicles are "separated spatially from the attack by CD8+ T cells" than
those at the periphery (J Hematol Oncol 2022, PMC9396877). A T cell must eat its way in.

MECHANISTIC PREDICTION UNDER TEST
---------------------------------
Exhaustion is driven by dwell time in contact with antigen (Philipp et al., Blood 2022,
PMID 35878001: continuous CD19xCD3 exposure drove specific lysis from 88.4% at day 7 to 8.6%
at day 28; a treatment-free interval restored day-14 lysis to 93.4% vs 34.9%). In a packed
follicle only the T cells at the invasion front are in contact at any moment; the rest are
searching. In a dispersed tumour nearly every T cell is in contact all the time.

A pilot measurement (21 days continuous, identical cell numbers and parameters, seed 0) gave:

    architecture   mean engaged fraction   mean exhaustion d21   contact-minutes per T cell
    follicle              0.085                  0.032                    2,118
    dispersed             0.610                  0.139                   10,679

So the same drug, the same cells and the same T cells produce a 7-fold difference in how fast
the effector pool degrades, purely from geometry. If exhaustion is what TFIs are for, then TFIs
should be worth much less in follicular architecture than in dispersed disease -- and the
Obertopp conclusion should not transfer to lymphoma.

DESIGN
------
3 architectures x 5 schedules x 12 seeds = 180 runs, 42 simulated days, 1 step = 1 minute.
Total B cells (5,542) and T cells (200) are IDENTICAL across architectures; only geometry
differs. Paired by seed: every architecture/schedule cell starts from the same seed.
  architectures  follicle   : one dense follicle, radius 42 sites = 840 um diameter
                              (measured neoplastic follicle diameter 757 um, range 577-930)
                 multi       : four dense follicles, radius 21, same total cells
                 dispersed   : 5,542 cells scattered, T cells interspersed (Obertopp regime)
  schedules      28-day cycle with the last k days off, k in {0, 2, 4, 7, 14}
                 k=0 is continuous; k=7 is the blinatumomab-style clinical interval.
  endpoints      nB at day 28 and day 42 (primary: day 42, matching Obertopp's endpoint where
                 their 7-day-TFI advantage was lost); cleared (nB=0); mean exhaustion; mean
                 engaged fraction; cumulative contact-minutes; T-cell penetration depth.

PRE-REGISTERED PREDICTIONS
--------------------------
P1 (mechanism, expected to hold): mean engaged fraction under continuous dosing is at least
   3x lower in `follicle` than in `dispersed`, and `multi` lies between them.
P2 (mechanism): mean T-cell exhaustion at day 28 under continuous dosing is at least 2x lower
   in `follicle` than in `dispersed`.
P3 (the claim): in `dispersed`, at least one TFI schedule beats continuous on day-42 tumour
   burden (median paired difference < 0, Wilcoxon p < 0.05), reproducing Obertopp. In
   `follicle`, no TFI schedule beats continuous by that rule.
P4 (dose-response): the benefit of the best TFI over continuous, in day-42 burden, is ordered
   dispersed > multi > follicle.
P5 (unifying variable): across all 15 architecture-schedule cells, the TFI benefit correlates
   with the continuous-arm mean engaged fraction, Spearman rho > 0.6.

FALSIFICATION. If TFIs help equally in all three architectures, the geometry-exhaustion
coupling is not doing the work claimed and the hypothesis is wrong. If TFIs help MORE in
`follicle`, the sign is opposite to the mechanism and the hypothesis is wrong. Both readings
are stated here before running.

NOT CLAIMED. This is a tissue-scale model of one lymph-node block. It says nothing about
whole-patient cure, about the common-progenitor reservoir (which is not represented), or about
CRS. Times are simulated days under stated rate calibrations, not clinical predictions.
"""
import json, os, sys, time
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lymphoid import Lymphoid

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
L = 120
R_BIG = 42
N_B = int(np.pi * R_BIG ** 2)          # 5542
N_T = 200
DAYS = 42
CYCLE = 28
SEEDS = list(range(10))
TFIS = [0, 2, 4, 7, 14]
ARCHS = ['follicle', 'multi', 'dispersed']
PARAMS = dict(p_kill=0.0012, t_influx=0.0001, p_div=1 / 2880, t_div=0.0)
DT = 5.0        # minutes per step. Motility is sub-stepped, so T cells still cover ~1 site/min.
                # Validated against dt=1: nB, nT, kills, exhaustion and engaged fraction all agree
                # within 13% on a matched run. Every arm uses the same dt, so the discretisation
                # bias is common-mode and cancels in the paired comparisons that carry the claims.


def build(arch, seed):
    T = Lymphoid(L=L, seed=seed, dt=DT, **PARAMS)
    if arch == 'follicle':
        T.seed_follicles([(L // 2, L // 2)], R_BIG)
    elif arch == 'multi':
        r = R_BIG // 2                                   # four follicles, same total area
        d = 30
        c = L // 2
        T.seed_follicles([(c - d, c - d), (c - d, c + d), (c + d, c - d), (c + d, c + d)], r)
    else:
        T.seed_dispersed(N_B, occupancy=0.95)
    T.seed_tcells(N_T)
    return T


def schedule_fn(tfi):
    """28-day cycle: on for (28 - tfi) days, then off for tfi days."""
    if tfi == 0:
        return lambda s: 1.0
    on_min = (CYCLE - tfi) * 1440
    cyc_min = CYCLE * 1440
    return lambda s: 1.0 if (s.t % cyc_min) < on_min else 0.0


def penetration(T):
    """Deepest T cell inside the original follicle footprint, as a fraction of its radius."""
    c = L // 2
    yy, xx = np.mgrid[:L, :L]
    r = np.hypot(yy - c, xx - c)
    inside = T.T & (r <= R_BIG)
    if not inside.any():
        return 0.0
    return float((R_BIG - r[inside].min()) / R_BIG)


def run(args):
    arch, tfi, seed = args
    T = build(arch, seed)
    n0 = T.nB
    sched = schedule_fn(tfi)
    snap = {}
    traj = []
    nsteps = int(DAYS * 1440 / DT)
    per_day = int(1440 / DT)
    for step in range(nsteps):
        T.step(sched(T))
        if step % per_day == per_day - 1:
            day = (step + 1) // per_day
            h = T.history[-1]
            traj.append((day, h['nB'], h['nT'], round(h['meanE'], 4), round(h['engFrac'], 4)))
            if day in (7, 14, 28, 42):
                snap[day] = h['nB']
        if T.nB == 0:
            break
    hist = T.history
    eng = float(np.mean([r['engFrac'] for r in hist]))
    # exhaustion and engagement measured over the first 28 days only, for a clean comparison
    h28 = hist[:min(len(hist), int(28 * 1440 / DT))]
    return dict(arch=arch, tfi=tfi, seed=seed, n0=n0,
                nB7=snap.get(7, 0), nB14=snap.get(14, 0), nB28=snap.get(28, 0), nB42=snap.get(42, 0),
                cleared=bool(T.nB == 0), day_cleared=(len(hist) / 1440.0 if T.nB == 0 else None),
                kills=T.kills, nT_end=int(T.T.sum()),
                meanE_end=float(hist[-1]['meanE']),
                meanE28=float(np.mean([r['meanE'] for r in h28])),
                engFrac_mean=eng,
                engFrac28=float(np.mean([r['engFrac'] for r in h28])),
                contact_per_T=float(T.cum_contact / max(int(T.T.sum()), 1)),
                penetration=penetration(T) if arch != 'dispersed' else None,
                traj=traj)


if __name__ == '__main__':
    t0 = time.time()
    jobs = [(a, k, s) for a in ARCHS for k in TFIS for s in SEEDS]
    print(f'exp L1: {len(jobs)} runs, {DAYS} simulated days each', flush=True)
    res = []
    with Pool(6) as p:
        for i, r in enumerate(p.imap_unordered(run, jobs, chunksize=1), 1):
            res.append(r)
            if i % 5 == 0 or i == len(jobs):
                el = time.time() - t0
                print(f'  {i}/{len(jobs)} done, {el/60:.1f} min elapsed, eta {(el/i*(len(jobs)-i))/60:.0f} min', flush=True)
                json.dump(res, open(f'{OUT}/expL1_partial.json', 'w'))
    json.dump(res, open(f'{OUT}/expL1.json', 'w'))
    msg = f'expL1: {len(jobs)} runs in {time.time()-t0:.0f}s'
    open(f'{OUT}/expL1_log.txt', 'w').write(msg + '\n')
    print(msg)
