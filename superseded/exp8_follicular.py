"""
Experiment 8 (pre-registered before running): does the spatial mechanism of adaptive therapy survive
the architecture of a lymphoma?

Motivation. Every result in this project (HANDOFF.md, DISCOVERY.md) concerns a resistant clone buried
in ONE compact mass whose interior is hypoxic. Its benefit from size-triggered dosing decomposes into
(i) confinement: the clone can only grow into drug-cleared space, so v_thr = 1.21 * dbar * v_free, and
(ii) the quiescence jail: while the core is hypoxic the clone cannot divide at all, and a cycle-
specific drug cannot touch it. Non-Hodgkin lymphoma has neither property by default. Follicular
lymphoma grows as many small follicles (~100-300 um) inside a vascular lymph node, each too small for
a hypoxic core; the front-line agents (anti-CD20 antibodies, T-cell engagers, CAR-T) kill regardless
of cell cycle; and aggressive subtypes (DLBCL, Burkitt) have very high turnover, whereas follicular
lymphoma has low turnover. Clinically, follicular lymphoma is the one NHL routinely managed with
watch-and-wait and intermittent treatment, i.e. de facto adaptive scheduling, yet the long-term
Ardeshna trial (Lancet Haematol 2025) found CONTINUOUS rituximab maintenance gave the longest time
to next treatment with no resistance penalty. The model makes a specific prediction about why.

Feasibility facts (follicular.py, seed 0, L = 160, N0 = 5,000): M follicles on a grid do not merge at
N0; hypoxic fraction at N0 = 0.77 (M=1), 0.54 (M=4), 0.09 (M=9), 0.006 (M=16); follicle radius 40,
20, 13, 10. At apoptosis 0.01 the compact mass hollows (as in exp5), so the clone is implanted at the
DEEPEST VIABLE POINT of the chosen nodule (maximum of the distance-to-empty transform), which is the
centre for a solid nodule and the mid-ring for a hollowed one. Its depth is recorded per run.

Design. M in {1, 4, 9, 16} (same total N0, L = 160 throughout, so geometry is the only change)
x drug class {cycle-specific, kills-quiescent} x apoptosis {0.0005 (a/p = 0.013, follicular-like),
0.01 (a/p = 0.25, aggressive-like)} x arm {continuous, thr0.5}, dose 0.8, compact 4-cell clone with
resist 0.9, no new mutations, 24 tumours regrown per (M, a). 768 runs. Endpoints as PILOT_PLAN.md:
TTR = first step with R >= 2,500, TTP = N >= 6,000, horizon 6,000, paired by tumour. Front fit
sqrt(R/pi) = v (t - t0) over 100 <= R <= 2,500 gives v and the release lag t0 per run.
Cure sub-experiment: mut_rate = 0.002, no implanted clone, continuous drug only, M in {1, 16} x drug
class, a = 0.0005, 24 tumours: 96 runs. Records eradication and time to eradication.

Pre-registered predictions (claim rule as always: median paired difference of the right sign AND
Wilcoxon p < 0.05; "benefit" = TTR_thr - TTR_cont per tumour):
  F1 Jail. At a = 0.0005 with the cycle-specific drug, the release-lag gain (t0_thr - t0_cont from
     the front fit, medians) falls monotonically with M: > 100 steps at M = 1, < 40 steps at M = 16.
     Across all threshold runs, Spearman(hypoxic fraction at N0, lag gain) > 0.6.
  F2 Confinement scales with jacket depth. With the kills-quiescent drug (no jail), a = 0.0005, the
     median benefit is proportional to implant depth: benefit(M) / benefit(1) lies within a factor
     of 1.5 of depth(M) / depth(1) for M = 4, 9, 16, with the law's slope about 10 steps per site of
     depth (1/(0.6 v_free) - 1/(0.85 v_free) = 10.5 at v_free = 0.047).
  F3 Drug class matters only where there is a jail. At a = 0.0005 the benefit lost by switching to
     the kills-quiescent drug is > 150 steps at M = 1 (exp2 found 189) and < 50 steps at M = 16.
  F4 Aggressive-lymphoma regime. At a = 0.01 and M = 16 the benefit fails the claim rule or its
     median is < 60 steps (< 10% of the continuous-arm TTR): the spatial benefit is gone for a
     high-turnover multi-follicular tumour, whatever the drug class.
  F5 Cure. Under continuous drug with mutation on, the compact mass shelters the tumour: eradication
     is more frequent at M = 16 than at M = 1 by >= 6 of 24 for the cycle-specific drug, and among
     eradicated tumours the median time to eradication at M = 16 is < 0.7 x that at M = 1. With the
     kills-quiescent drug eradication at M = 1 rises by >= 4 of 24 relative to the cycle-specific
     drug, and the M effect narrows.
  Reading if F1-F4 hold: size-triggered dosing has nothing to offer a follicular, cycle-agnostic,
  low-turnover lymphoma beyond a small confinement term of order 100 steps, and nothing at all to
  an aggressive one; continuous dosing is favoured, in line with Ardeshna 2025. If F2 fails with a
  benefit that does NOT shrink with depth, the confinement mechanism is not geometric and the
  project's central law is wrong in this regime; both readings are stated in advance.
"""
import json, os, time
import numpy as np
from multiprocessing import Pool
from scipy import ndimage
import pilot
from tumor import Tumor
from frontfit import fit
import follicular as F

OUT = os.path.join(os.path.dirname(__file__), 'results')
MS = [1, 4, 9, 16]; RATES = [0.0005, 0.01]; DRUGS = [False, True]; K = 4


def implant_deepest(T, M, K=K):
    lab, ncomp = ndimage.label(T.alive, structure=np.ones((3, 3)))
    sizes = ndimage.sum(T.alive, lab, range(1, ncomp + 1))
    big = [i + 1 for i, s in enumerate(sizes) if s >= 0.5 * pilot.N0 / M]
    if not big: big = [int(np.argmax(sizes)) + 1]
    cents = ndimage.center_of_mass(T.alive, lab, big)
    lb = big[int(np.argmin([(cy - F.L / 2) ** 2 + (cx - F.L / 2) ** 2 for cy, cx in cents]))]
    comp = lab == lb
    edt = ndimage.distance_transform_edt(T.alive)          # distance to nearest empty site
    edt[~comp] = -1
    iy, ix = np.unravel_index(np.argmax(edt), edt.shape)
    depth = float(edt[iy, ix])
    ys0, xs0 = np.nonzero(comp)
    d = (ys0 - iy) ** 2 + (xs0 - ix) ** 2
    pick = np.argsort(d)[:K]
    T.resist[ys0[pick], xs0[pick]] = 0.9
    T._relax_o2()
    return ys0[pick], xs0[pick], depth, float(np.sqrt(comp.sum() / np.pi)), int(ncomp)


def treat(T, arm, dose=pilot.DOSE):
    state = dict(on=True); E = 0.0; ttp = ttr = None; traj = []; t_erad = None
    for st in range(pilot.HORIZON):
        if arm == 'continuous': dse = dose
        else:
            n = T.n
            if state['on'] and n < 0.5 * pilot.N0: state['on'] = False
            elif not state['on'] and n >= pilot.N0: state['on'] = True
            dse = dose if state['on'] else 0.0
        T.step(dse); E += dse
        h = T.history[-1]; n, r = h['n'], h['n_res']
        if st % 5 == 0: traj.append((st, n, r, h['n_sens_q'], round(E, 1), dse))
        if ttp is None and n >= pilot.PROG_N: ttp = st + 1
        if ttr is None and r >= pilot.PROG_R: ttr = st + 1
        if n == 0: t_erad = st + 1; break
        if ttp is not None and ttr is not None: break
    return dict(ttp=ttp, ttr=ttr, E=round(E, 1), t_erad=t_erad, eradicated=T.n == 0, traj=traj,
                r_end=T.history[-1]['n_res'], n_end=T.n)


def run(args):
    seed, M, a, kq, arm = args
    params = dict(pilot.BASE, mut_rate=0.0, apoptosis=a, drug_kills_quiescent=kq)
    T = F.make_follicular(seed, M, params)
    hyp0 = float((T.alive & (T.o2 < T.hypoxia)).sum() / T.n)
    T = T.copy(800000 + 1000 * seed + 100 * MS.index(M) + 10 * RATES.index(a) + 2 * int(kq) + (arm == 'threshold'))
    ys, xs, depth, rf, ncomp = implant_deepest(T, M)
    quies = int((T.o2[ys, xs] < T.hypoxia).sum()); tpre = T.t
    out = treat(T, arm)
    f = fit(out['traj']) or {}
    return dict(seed=seed, M=M, apoptosis=a, kills_quiescent=kq, arm=arm, t_pre=tpre,
                hypoxic0=hyp0, depth=depth, r_f=rf, ncomp=ncomp, quiescent=quies,
                v=f.get('v'), t0=f.get('t0'), r2=f.get('r2'), **out)


def run_cure(args):
    seed, M, kq = args
    params = dict(pilot.BASE, mut_rate=0.002, apoptosis=0.0005, drug_kills_quiescent=kq)
    T = F.make_follicular(seed, M, params)
    hyp0 = float((T.alive & (T.o2 < T.hypoxia)).sum() / T.n)
    r0 = int((T.alive & (T.resist > 0.5)).sum()); mr0 = float(T.resist[T.alive].mean())
    T = T.copy(900000 + 1000 * seed + 100 * MS.index(M) + 2 * int(kq))
    out = treat(T, 'continuous')
    return dict(seed=seed, M=M, kills_quiescent=kq, hypoxic0=hyp0, r_start=r0, mean_resist0=mr0, **out)


if __name__ == '__main__':
    t0 = time.time()
    jobs = [(s, M, a, kq, arm) for M in MS for a in RATES for kq in DRUGS for s in pilot.SEEDS for arm in ('continuous', 'threshold')]
    with Pool(6) as p:
        res = p.map(run, jobs, chunksize=2)
    json.dump(res, open(f'{OUT}/exp8.json', 'w'))
    print(f'exp8 main: {len(jobs)} runs in {time.time()-t0:.0f}s', flush=True)
    t1 = time.time()
    jobs = [(s, M, kq) for M in (1, 16) for kq in DRUGS for s in pilot.SEEDS]
    with Pool(6) as p:
        res = p.map(run_cure, jobs, chunksize=2)
    json.dump(res, open(f'{OUT}/exp8_cure.json', 'w'))
    open(f'{OUT}/exp8_log.txt', 'w').write(f'exp8: 768 main runs in {t1-t0:.0f}s; 96 cure runs in {time.time()-t1:.0f}s\n')
    print(f'exp8 cure: {len(jobs)} runs in {time.time()-t1:.0f}s', flush=True)
