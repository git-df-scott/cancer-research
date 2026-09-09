"""
Phase E: fit latent-dynamics models against the CORRECTED observation model.

WHY THE "LINEAR IS IMPOSSIBLE" RESULT IS NOW ITSELF SUPERSEDED
--------------------------------------------------------------
Earlier work argued that a linear exhaustion mechanism cannot reproduce Philipp's curve, by
inverting `lysis = naive x (1 - E)` and observing that the data demand E(14)/E(7) between 5.6 and
unbounded while linear accrual forces exactly 2.00. That argument was reported as surviving the
probe bug because it depends on no simulation.

It does not survive Phase D. It assumed the observation model
    Y = 1 - E
which is precisely the assumption Phase D rejects. The measured observable is 72 h target
depletion at E:T 1:1 normalised to a no-engager control, and that is a SATURATING function of
per-cell killing capacity, not a linear one. Measured directly (see below), the readout is flat at
100% for E from 0 to 0.7 at the shipped kill rate, then collapses -- a compressive nonlinearity
that can turn a linear latent decay into a sigmoid observable all by itself.

So the arithmetic was sound about the model it described and wrong about this system. Both the
linear and the threshold mechanisms must be re-fitted against the corrected observation model, and
neither can be ruled out a priori. That is what this module does.

p_kill IS A FREE PARAMETER HERE, AND HAS TO BE
-----------------------------------------------
`lymphoid.py` sets p_kill from Halle et al. (2-16 targets/CTL/day, in vivo). The calibration assay
is in vitro at E:T 1:1 over 72 h. Using an in vivo CTL killing rate as an in vitro assay parameter
is a mapping error of exactly the kind that produced the previous round of corrections, so it is
not carried over unexamined.

It also sets the readout's operating point, measured over a synthetic exhaustion sweep:

    p_kill    kills/T/day   E=0.0   E=0.3   E=0.6   E=0.9
    4.2e-03      6.05       100.0   100.0   100.0    59.7     <- shipped value: saturated, blind
    1.0e-03      1.44       100.0    96.9    85.3    12.8
    3.0e-04      0.43        80.9    66.2    44.7     5.6
    1.0e-04      0.14        48.4    31.9    19.7     2.2

Philipp's range is 88.4% down to 8.6%, so the assay only has discriminating power near
p_kill ~ 5e-4. At the shipped value the observable is nearly independent of exhaustion across the
entire biologically relevant range, which would make any fit meaningless.

This is an identifiability statement, not a licence to tune: p_kill and the exhaustion parameters
are constrained jointly by these four numbers, and reporting a fitted exhaustion rate without its
accompanying p_kill would be reporting half a result.

FIXED BEFORE FITTING
--------------------
Targets, weighting and acceptance are declared here and must not be changed after seeing results.

  targets     Philipp Fig 2E and Fig 3E, four values, n=6 donors, mean +/- SEM
  weighting   uniform, unweighted SSE in percentage-point units
  held out    d14_tfi is HELD OUT (Phase F). Fitting uses the three continuous-arm points only.
              d14_tfi is a recovery observation; nothing in the fit informs recover_tau, so it is
              a genuine out-of-sample prediction rather than a re-description.
  tolerance   predeclared: |predicted - 93.4| <= 15 percentage points. Philipp report mean +/- SEM
              over n=6 donors; the SEM is not given numerically in the text available, so 15 points
              is a stated practical margin, not a derived one, and is recorded as such.
"""
import json, os, sys, itertools
from multiprocessing import Pool
import numpy as np

import philipp_assay as PA
from lymphoid import Lymphoid
from exhaustion import ThresholdLymphoid

OUT = os.environ.get('FIT_OUT', 'results/calib2')
os.makedirs(OUT, exist_ok=True)

FIT_TARGETS = ('d7_cont', 'd14_cont', 'd28_cont')   # d14_tfi is held out
HELD_OUT = 'd14_tfi'
HELD_OUT_TOL = 15.0


def fit_loss(meas):
    """SSE over the FIT targets only. The held-out point never enters this."""
    return sum((meas[k] - PA.PHILIPP[k]) ** 2 for k in FIT_TARGETS)


def cell(job):
    model, p_kill, tonic, extra = job
    tag = f"{model}_pk{p_kill:.2e}_t{tonic:.2e}" + \
          ("" if not extra else "_" + "_".join(f"{k}{v}" for k, v in sorted(extra.items())))
    path = f'{OUT}/{tag}.json'
    if os.path.exists(path):
        return path
    cls = {'M0': Lymphoid, 'M1': ThresholdLymphoid}[model]
    params = dict(p_kill=p_kill, exhaust_tonic=tonic, exhaust_per_kill=0.0, **extra)
    per_seed = [PA.evaluate(model_cls=cls, seed=s, assay_seed=900 + s, **params) for s in (0, 1)]
    meas = {k: float(np.mean([m[k] for m in per_seed])) for k in PA.PHILIPP}
    spread = {k: float(np.std([m[k] for m in per_seed])) for k in PA.PHILIPP}
    row = dict(model=model, p_kill=p_kill, exhaust_tonic=tonic, extra=extra,
               measured=meas, seed_sd=spread, target=PA.PHILIPP,
               fit_loss=fit_loss(meas), full_loss=PA.loss(meas),
               err={k: round(meas[k] - PA.PHILIPP[k], 1) for k in PA.PHILIPP},
               held_out=HELD_OUT,
               held_out_err=round(meas[HELD_OUT] - PA.PHILIPP[HELD_OUT], 1),
               held_out_pass=bool(abs(meas[HELD_OUT] - PA.PHILIPP[HELD_OUT]) <= HELD_OUT_TOL))
    with open(path, 'w') as f:
        json.dump(row, f, indent=1)
    return path


def grid(model, mode='base'):
    # Base grid. Its optimum landed on the boundary in BOTH dimensions (highest p_kill, lowest
    # tonic), so it did not bracket the optimum and its best cell must not be quoted as a fit.
    pks = [3e-4, 5e-4, 7e-4, 1e-3]
    tonics = [2.5e-5, 5e-5, 1e-4, 2e-4]
    if mode == 'extend':
        pks = [1e-3, 2e-3, 4e-3]
        tonics = [6e-6, 1.2e-5, 2.5e-5]
    if mode == 'clean':
        # Phase 1 clean M1 reselection. See M1_RESELECTION_PREREG.md, committed before running.
        # theta and hill ranges are defined from the reachable domain and structural limits, NOT
        # from the previous optimum, which came from the superseded destructive-probe scan.
        pks = [5e-4, 1e-3, 2e-3]
        tonics = [2.5e-5, 5e-5, 1e-4]
    if model == 'M0':
        return [('M0', pk, t, {}) for pk, t in itertools.product(pks, tonics)]
    thetas = [0.15, 0.30, 0.50, 0.70, 0.90] if mode == 'clean' else [0.5]
    hills = [1.0, 2.0, 4.0, 8.0, 16.0] if mode == 'clean' else [4.0]
    return [('M1', pk, t, {'theta': th, 'hill': h})
            for pk, t, th, h in itertools.product(pks, tonics, thetas, hills)]


if __name__ == '__main__':
    models = sys.argv[1].split(',') if len(sys.argv) > 1 else ['M0']
    nproc = int(sys.argv[2]) if len(sys.argv) > 2 else os.cpu_count()
    mode = sys.argv[3] if len(sys.argv) > 3 else 'base'
    jobs = [j for m in models for j in grid(m, mode)]
    print(f'{len(jobs)} cells on {nproc} procs', flush=True)
    with Pool(nproc) as p:
        for i, path in enumerate(p.imap_unordered(cell, jobs), 1):
            r = json.load(open(path))
            m = r['measured']
            print(f"{i}/{len(jobs)} {r['model']} pk={r['p_kill']:.1e} t={r['exhaust_tonic']:.1e} | "
                  f"d7={m['d7_cont']:5.1f} d14={m['d14_cont']:5.1f} d28={m['d28_cont']:5.1f} "
                  f"|| held-out d14tfi={m['d14_tfi']:5.1f} ({'PASS' if r['held_out_pass'] else 'fail'})"
                  f" | fit_loss={r['fit_loss']:8.1f}", flush=True)
