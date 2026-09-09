"""
Experiment R (PRE-REGISTERED, written before any run): is the model broken, or was the
experiment misconfigured?

THE PROBLEM THIS EXISTS TO SETTLE
---------------------------------
`docs/findings/FINDINGS.md` classifies experiment L1 as UNINFORMATIVE because its positive control failed: the
dispersed arm was meant to reproduce the reference result (short treatment-free intervals beat
continuous dosing in well-mixed, leukaemia-like disease) and produced the opposite ranking.
`docs/findings/FINDINGS.md` and `docs/plan/HANDOFF.md` both attribute that failure to exhaustion being miscalibrated --
the model reaches 12% functional loss at day 28 where Philipp et al. measured ~90% -- and make
recalibration the blocking first task.

That diagnosis was never tested. Reading the reference's published methods (PMC12667981) against
`exp_schedule.py` shows the dispersed arm was not a replication of the reference in the first
place. It differs in four ways:

  1. SCHEDULE SHAPE. Reference: "recurring treatment-free intervals (2-7 days, respectively),
     alternating with periods of TCE dosing" -- a duty cycle -- plus a uniform 14-day rest applied
     to every arm including continuous. L1: one interruption at the end of a 28-day cycle, and the
     continuous arm dosed for all 42 days. See `schedules.py`. Measured exposure over 42 days:
     reference CONT 0.667, L1 tfi=0 1.000. L1's continuous arm received 1.5x the reference's.
  2. E:T RATIO. Reference 1:4 (2,400 T : 9,600 tumour). L1 1:27.7 (200 T : 5,542 B).
  3. T-CELL PROLIFERATION. Reference: active T cells divide after killing (5 h refractory).
     L1: `t_div=0.0`, disabled entirely.
  4. T-CELL INFLUX. Reference: none, fixed population. L1: `t_influx=1e-4`, 200 -> ~860 per run.

Items 1 and 2 are each independently capable of producing the observed ranking:

  - Under the L1 encoding, TFI_k IS "continuous minus k days of drug". It cannot beat continuous
    unless exhaustion repays the lost exposure, so TFI0 > TFI2 > TFI4 > TFI7 > TFI14 is close to
    guaranteed by construction -- which is exactly what L1 reported, in all three architectures.
  - At 1:27.7 the system is target-rich and effector-poor. Every T cell has a target in contact
    essentially always, killing is effector-limited, and 25% less drug-time costs ~25% of the
    killing directly. Exhaustion cannot become the binding constraint at that ratio under ANY
    calibration.

So the exhaustion gap may be real and still not be the cause. This experiment separates them.

PRE-REGISTERED PREDICTIONS
--------------------------
R1 (replication). Under faithful reference conditions -- duty-cycle schedules, 1:4 E:T, 50%
   occupancy, T-cell proliferation on, no influx -- at least one duty-cycle arm achieves lower
   median tumour burden than CONT at day 28.
   FAILS IF: CONT is <= every duty-cycle arm at day 28.

R2 (exhaustion magnitude). Under CONT at 1:4 E:T with no influx, mean exhaustion exceeds 0.50 by
   day 16. The reference reports >80%; 0.50 is the weaker claim this model must clear for
   exhaustion to be capable of binding at all.
   FAILS IF: mean exhaustion at day 16 under CONT is < 0.50.

R3 (attribution). In the one-at-a-time ablation from the L1 configuration toward the reference,
   the single largest movement in the CONT-minus-best-TFI contrast is produced by the schedule
   encoding or the E:T ratio -- not by influx.
   FAILS IF: influx alone accounts for the largest movement.

R4 (assumption robustness). The R1 verdict is unchanged across on_days in {2, 3, 5, 7}.
   FAILS IF: R1 holds for some on_days and fails for others. Then the duty cycle is doing the work
   and the result is reported as assumption-dependent, not as a replication.

DECISION RULE, FIXED IN ADVANCE
-------------------------------
  R1 holds and R2 holds  -> control passes. The model was sound; L1 was misconfigured.
                            docs/findings/FINDINGS.md's root-cause attribution is wrong and gets corrected.
                            Re-run the architecture comparison with correct schedules.
  R1 holds, R2 fails     -> schedule encoding was the bug; exhaustion is still undercooked.
                            Recalibrate, but the scope is narrowed and L1's diagnosis was still
                            wrong about the cause.
  R1 fails               -> genuine model failure under faithful conditions. Proceed to full
                            exhaustion recalibration as docs/plan/HANDOFF.md originally specified.

In every branch the ablation table is the deliverable, because it attributes the failure
quantitatively instead of asserting a cause. That is the thing L1 did not do.

WHAT THIS EXPERIMENT DOES NOT DO
--------------------------------
It does not tune any parameter toward any ranking. Every configuration below is either taken from
the reference's stated methods or preserved verbatim from L1. The one quantity the reference does
not state -- ON-days between OFF-intervals -- is swept, not chosen (R4).

DEVIATIONS FROM THE REFERENCE, RECORDED NOT HIDDEN
--------------------------------------------------
  - Grid. Reference is 150x160 = 24,000 sites. `Lymphoid` is square, so L=155 (24,025 sites) is
    used. 0.1% difference in area, common-mode across arms.
  - T-cell division. The reference gates division on a 5-hour refractory period since the last
    division. `Lymphoid` has no per-cell division clock; it uses a per-minute hazard gated on
    drug presence and target adjacency. `t_div` is set so mean division rate is comparable. This
    is an approximation and is listed here as such.
  - Exhaustion mechanism. The reference accumulates PD-1 and applies a hard threshold, after
    which the cell cannot kill and dies with ~40% probability within 10 days. `Lymphoid` uses a
    continuous exhaustion scalar scaling the kill hazard. These are different mechanisms. R2 is
    the check on whether the difference matters, and it is a pre-registered endpoint, not a
    post-hoc excuse.
"""
import json, os, sys, itertools
from multiprocessing import Pool
import numpy as np

from lymphoid import Lymphoid
import schedules as S

OUT = 'results/expR'
os.makedirs(OUT, exist_ok=True)

# ----------------------------------------------------------------- configurations
L_REF = 155                     # 24,025 sites ~ the reference's 150x160 = 24,000
N_TOTAL = 12000                 # 50% occupancy
N_B_REF = 9600                  # 80% of occupied
N_T_REF = 2400                  # 20% of occupied -> E:T 1:4
DAYS = 42
DOSE_DAYS = 28
DT = 5.0
SEEDS = list(range(10))

# L1's configuration, preserved verbatim from exp_schedule.py
L1 = dict(L=120, n_B=5542, n_T=200, occupancy=0.95,
          p_kill=0.0012, t_influx=1e-4, t_div=0.0,
          sched_kind='legacy', dose_days=42)

# The reference's configuration, from its published methods
REF = dict(L=L_REF, n_B=N_B_REF, n_T=N_T_REF, occupancy=0.5,
           p_kill=0.0012, t_influx=0.0, t_div=1.0 / 2880,
           sched_kind='duty', dose_days=DOSE_DAYS)


def build(cfg, seed):
    m = Lymphoid(L=cfg['L'], seed=seed, dt=DT,
                 p_kill=cfg['p_kill'], t_influx=cfg['t_influx'],
                 t_div=cfg['t_div'], p_div=1 / 2880)
    m.seed_dispersed(cfg['n_B'], occupancy=cfg['occupancy'])
    m.seed_tcells(cfg['n_T'])
    return m


def make_sched(cfg, arm):
    """arm is ('cont',) or ('duty', on_days, off_days) or ('legacy', tfi)."""
    if arm[0] == 'cont':
        if cfg['sched_kind'] == 'legacy':
            return S.legacy_l1(0)
        return S.continuous(dose_days=cfg['dose_days'])
    if arm[0] == 'duty':
        return S.duty_cycle(arm[1], arm[2], dose_days=cfg['dose_days'])
    return S.legacy_l1(arm[1])


def run_one(job):
    tag, cfg, arm, seed = job
    path = f'{OUT}/{tag}__{"_".join(str(a) for a in arm)}__s{seed}.json'
    if os.path.exists(path):
        return path
    m = build(cfg, seed)
    sched = make_sched(cfg, arm)
    n0 = m.nB
    steps_per_day = int(1440 / DT)
    snap = {}
    for day in range(DAYS):
        m.run(steps_per_day, schedule=sched, record_every=10**9, stop_when_clear=False)
        nT = int(m.T.sum())
        snap[day + 1] = dict(nB=m.nB, nT=nT,
                             meanE=float(m.E[m.T].mean()) if nT else 0.0,
                             engFrac=float(m._n_engaged / nT) if nT else 0.0)
    rec = dict(tag=tag, arm=list(arm), seed=seed, n0=n0,
               nB16=snap[16]['nB'], nB28=snap[28]['nB'], nB42=snap[42]['nB'],
               meanE16=snap[16]['meanE'], meanE28=snap[28]['meanE'],
               kills=int(m.kills), nT_end=int(m.T.sum()),
               exposure=S.exposure_fraction(sched, DAYS, DT),
               snap=snap)
    with open(path, 'w') as f:
        json.dump(rec, f)
    return path


def jobs_replication(on_days_sweep=(2, 3, 5, 7), off_days_sweep=(2, 3, 7)):
    """R1/R2/R4: faithful reference conditions, duty cycle swept."""
    out = [('ref', REF, ('cont',), s) for s in SEEDS]
    for on_d, off_d in itertools.product(on_days_sweep, off_days_sweep):
        out += [('ref', REF, ('duty', on_d, off_d), s) for s in SEEDS]
    return out


def jobs_ablation():
    """R3: one-at-a-time toggles from L1 toward REF. Each row changes exactly one factor."""
    out = []
    variants = {
        'L1': L1,
        'L1+sched':  {**L1, 'sched_kind': 'duty', 'dose_days': DOSE_DAYS},
        'L1+ET':     {**L1, 'L': L_REF, 'n_B': N_B_REF, 'n_T': N_T_REF, 'occupancy': 0.5},
        'L1+prolif': {**L1, 't_div': 1.0 / 2880},
        'L1+noinflux': {**L1, 't_influx': 0.0},
    }
    for tag, cfg in variants.items():
        if cfg['sched_kind'] == 'legacy':
            arms = [('legacy', 0), ('legacy', 2), ('legacy', 7)]
        else:
            arms = [('cont',), ('duty', 5, 2), ('duty', 5, 7)]
        for arm in arms:
            out += [(tag, cfg, arm, s) for s in SEEDS]
    return out


if __name__ == '__main__':
    # HARD GATE. R and T previously constructed Lymphoid with shipped defaults, so running
    # them would have used uncalibrated parameters while appearing to run the calibrated
    # experiment (Codex review, finding 3). They now refuse unless a VALIDATED calibration
    # artifact exists. A good fit is not enough: validation requires a provenance-clean
    # external prediction, which provenance.py shows is currently unavailable.
    from calibration import require_validated
    _calib = require_validated('experiment R')

    what = sys.argv[1] if len(sys.argv) > 1 else 'all'
    nproc = int(sys.argv[2]) if len(sys.argv) > 2 else os.cpu_count()
    jobs = []
    if what in ('all', 'rep'):
        jobs += jobs_replication()
    if what in ('all', 'abl'):
        jobs += jobs_ablation()
    print(f'{len(jobs)} runs on {nproc} procs', flush=True)
    with Pool(nproc) as p:
        for i, path in enumerate(p.imap_unordered(run_one, jobs), 1):
            print(f'{i}/{len(jobs)} {os.path.basename(path)}', flush=True)
