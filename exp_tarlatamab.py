"""
Experiment T (PRE-REGISTERED, written before any ABM run): does a half-life-extended engager
dosed Q2W deliver a treatment-free interval at all?

THE CLAIM UNDER TEST
--------------------
Every treatment-free-interval (TFI) result in this field descends from blinatumomab: terminal
half-life ~2 hours, continuous infusion, so stopping the pump takes concentration to zero within
hours and T cells genuinely disengage. Tarlatamab is a half-life-extended BiTE -- 5.8 days
(DeLLphi-300, PMID 39589690) to 11.2 days (population PK, PMID 40261494) -- given 10 mg every 14
days after a 1 mg step-up prime.

With an 11-day half-life on a 14-day interval, concentration never approaches zero. If exhaustion
is driven by cumulative engaged dwell time, the extension that made tarlatamab dosable may have
removed the disengagement window the TFI literature assumes exists -- while clinical resistance
to tarlatamab with DLL3 RETAINED, accompanied by Treg enrichment and CD8 TIM3/LAG3 exhaustion
markers, is now documented (PMC13082099; PMID 41903703). A 2026 review states the TFI approach
"has not yet been formally evaluated in tarlatamab-treated SCLC patients" (PMC13457367).

PRE-REGISTERED PREDICTIONS
--------------------------
T1 (occupancy floor). At steady state under 10 mg Q2W, target occupancy falls by less than 25%
   from its peak across the inter-dose interval.
   FAILS IF: the drop exceeds 25%.

T2 (dwell-time equivalence). Cumulative engaged T-cell-minutes over 84 days under the clinical
   Q2W schedule are within 10% of a matched-AUC continuous infusion.
   FAILS IF: the difference exceeds 10%.

T3 (exhaustion equivalence). Mean exhaustion at day 84 under clinical Q2W is not lower than under
   matched-AUC continuous infusion by more than 0.05 in absolute terms.
   FAILS IF: Q2W is lower by more than 0.05.

T4 (fractionation has no purchase). A fractionated schedule at matched total exposure (5 mg
   weekly) does not reduce day-84 exhaustion relative to clinical Q2W by more than 0.05.
   FAILS IF it does -- which would mean peak-shaving matters and the trough argument is wrong.

T5 (where it breaks). There exists an EC50 above which T1 fails. Report that boundary explicitly.

RESULT ALREADY OBTAINED FROM THE PK LAYER ALONE, BEFORE ANY ABM RUN
-------------------------------------------------------------------
T1 is EC50-dependent, and the dependence is sharp. Computed from `pk.py` at steady state:

    EC50        occupancy drop from peak, per inter-dose interval
                 t1/2 = 5.8 d      t1/2 = 11.2 d
      0.1 nM        2.9%              2.1%
      1.0 nM       22.3%             17.3%      <- in vitro anchor (~1 nM, PMC13082099)
     10.0 nM       69.1%             62.1%
    100.0 nM       87.4%             83.8%

So T1 HOLDS at the in-vitro-anchored EC50, but only just -- 17-22% against a 25% threshold -- and
FAILS decisively at 10 nM and above. The hypothesis therefore lives or dies on a quantity that is
NOT established in vivo: the in vitro co-culture EC50 of ~1 nM is a different measurement in a
different system, and using it as an in vivo potency is an assumption, not a citation.

This is recorded here, before the ABM runs, so it cannot later be presented as a robust prediction.
Any final classification must be at most D (robust computational prediction) CONDITIONAL on EC50,
and the conditionality must appear in the claim itself, not in a footnote.

WHAT THE ABM ADDS BEYOND THIS TABLE
-----------------------------------
The occupancy table is a statement about drug exposure. It says nothing about whether T cells are
in contact when the drug is present, which is the quantity exhaustion actually integrates. In
spatially structured disease a T cell can be saturated with drug and still not be engaged. T2/T3
are therefore not implied by T1, and the ABM is what separates them.
"""
import json, os, sys, itertools
from multiprocessing import Pool
import numpy as np

from lymphoid import Lymphoid
import pk as PK

OUT = 'results/expT'
os.makedirs(OUT, exist_ok=True)

MW_KDA = 105.0                    # half-life-extended BiTE, ~105 kDa
DAYS = 84
DT = 5.0
SEEDS = list(range(8))
EC50_nM_GRID = [0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0]
THALF_GRID = [PK.THALF_D300, PK.THALF_POPPK]


def nM_to_mgL(x_nM):
    return x_nM * MW_KDA * 1e-3


def regimen(kind, t_half):
    """Build the dose history for a named schedule. All deliver comparable total mg where stated."""
    m = PK.TwoCompartmentPK(t_half_terminal=t_half)
    if kind == 'q2w':
        m.add_regimen(until_day=DAYS)
    elif kind == 'weekly_half':
        m.add_dose(0.0, 1.0)
        t = 7.0
        while t <= DAYS:
            m.add_dose(t, 5.0)
            t += 7.0
    elif kind == 'weekly_matched':
        # Fractionation at MATCHED windowed exposure. The plain weekly_half arm delivers ~11% less
        # AUC over the horizon, so any advantage it showed could be a dose effect rather than a
        # timing one. T4 asks whether peak-shaving helps at equal exposure, so the comparator has
        # to be equal-exposure. Linear PK means a single scale factor suffices.
        m.add_dose(0.0, 1.0)
        tt = 7.0
        while tt < DAYS:
            m.add_dose(tt, 5.0)
            tt += 7.0
        ref = PK.TwoCompartmentPK(t_half_terminal=t_half).add_regimen(until_day=DAYS)
        m.scale_doses(ref.auc(0.0, DAYS) / m.auc(0.0, DAYS))
    elif kind == 'q4w_double':
        m.add_dose(0.0, 1.0)
        t = 7.0
        while t <= DAYS:
            m.add_dose(t, 20.0)
            t += 28.0
    elif kind == 'infusion':
        # Matched on AUC over the SIMULATED WINDOW, not on total milligrams. See
        # pk.add_infusion_auc_matched for why dose-matching biased this comparator by ~9%.
        ref = PK.TwoCompartmentPK(t_half_terminal=t_half).add_regimen(until_day=DAYS)
        m.add_infusion_auc_matched(ref.auc(0.0, DAYS), 0.0, DAYS)
    elif kind == 'q2w_holiday':
        # genuine drug holiday: two cycles on, one cycle off, repeating
        t, i = 0.0, 0
        m.add_dose(0.0, 1.0)
        t = 7.0
        while t <= DAYS:
            if (i % 3) != 2:
                m.add_dose(t, 10.0)
            t += 14.0
            i += 1
    else:
        raise ValueError(kind)
    return m


def occupancy_stats(m, ec50_mgL, day_lo=56, day_hi=70):
    t = np.linspace(day_lo, day_hi, 5001)
    occ = m.occupancy(t, ec50_mgL)
    return dict(occ_min=float(occ.min()), occ_max=float(occ.max()),
                occ_mean=float(occ.mean()),
                drop_pct=float(100 * (1 - occ.min() / occ.max())) if occ.max() else 0.0)


def run_one(job):
    kind, t_half, ec50_nM, seed = job
    tag = f'{kind}__th{t_half}__ec{ec50_nM}__s{seed}'
    path = f'{OUT}/{tag}.json'
    if os.path.exists(path):
        return path
    ec50 = nM_to_mgL(ec50_nM)
    m = regimen(kind, t_half)
    sched = m.as_schedule(ec50)

    sim = Lymphoid(L=120, seed=seed, dt=DT, p_kill=0.0012, t_influx=1e-4, t_div=0.0)
    sim.seed_dispersed(5542, occupancy=0.95)
    sim.seed_tcells(200)

    steps_per_day = int(1440 / DT)
    snap = {}
    for day in range(DAYS):
        sim.run(steps_per_day, schedule=sched, record_every=10**9, stop_when_clear=False)
        nT = int(sim.T.sum())
        snap[day + 1] = dict(nB=sim.nB, nT=nT,
                             meanE=float(sim.E[sim.T].mean()) if nT else 0.0,
                             occ=float(sched(sim)))
    rec = dict(kind=kind, t_half=t_half, ec50_nM=ec50_nM, seed=seed,
               total_mg=m.total_mg(), auc_window=m.auc(0.0, DAYS),
               dwell_Tmin=float(sim.cum_engaged_min),
               nB84=sim.nB, meanE84=snap[DAYS]['meanE'], kills=int(sim.kills),
               occupancy=occupancy_stats(m, ec50), snap=snap)
    with open(path, 'w') as f:
        json.dump(rec, f)
    return path


def pk_only_table():
    """T1 and T5 need no ABM. Compute and write them first; they are the honest headline."""
    rows = []
    for t_half, ec_nM in itertools.product(THALF_GRID, EC50_nM_GRID):
        m = regimen('q2w', t_half)
        s = occupancy_stats(m, nM_to_mgL(ec_nM))
        s.update(t_half=t_half, ec50_nM=ec_nM, T1_holds=bool(s['drop_pct'] < 25.0))
        rows.append(s)
        print(f"t1/2={t_half:5.1f}d EC50={ec_nM:6.1f}nM  occ {s['occ_min']:.3f}->{s['occ_max']:.3f}"
              f"  drop {s['drop_pct']:5.1f}%  T1 {'HOLDS' if s['T1_holds'] else 'FAILS'}", flush=True)
    with open(f'{OUT}/pk_only.json', 'w') as f:
        json.dump(rows, f, indent=1)
    return rows


def jobs():
    kinds = ['q2w', 'infusion', 'weekly_matched', 'weekly_half', 'q4w_double', 'q2w_holiday']
    return [(k, PK.THALF_POPPK, ec, s)
            for k in kinds for ec in (1.0, 10.0) for s in SEEDS]


if __name__ == '__main__':
    what = sys.argv[1] if len(sys.argv) > 1 else 'pk'
    if what == 'pk':
        pk_only_table()
    else:
        js = jobs()
        nproc = int(sys.argv[2]) if len(sys.argv) > 2 else os.cpu_count()
        print(f'{len(js)} runs on {nproc} procs', flush=True)
        with Pool(nproc) as p:
            for i, path in enumerate(p.imap_unordered(run_one, js), 1):
                print(f'{i}/{len(js)} {os.path.basename(path)}', flush=True)
