"""Calibrate the THRESHOLD exhaustion mechanism against the same four Philipp targets.

Deliberately a separate file from `calibrate_exhaustion.py` so the linear-mechanism scan and this
one are directly comparable and neither can quietly acquire the other's assumptions. Same assay,
same targets, same loss. The only thing that differs is the functional form under test.
"""
import json, os, itertools, sys
from multiprocessing import Pool
import numpy as np

import calibrate_exhaustion as C
from exhaustion import ThresholdLymphoid

OUT = 'results/calib/threshold'
os.makedirs(OUT, exist_ok=True)


class ThresholdAssay(C.Assay):
    """C.Assay, but the model under test is ThresholdLymphoid."""

    def __init__(self, L=60, seed=0, target_occ=0.40, et_ratio=0.25,
                 theta=0.35, hill=6.0, **kw):
        self.m = ThresholdLymphoid(L=L, seed=seed, t_influx=0.0, t_div=0.0,
                                   p_div=0.0, p_death=0.0,
                                   theta=theta, hill=hill, **kw)
        self.L = L
        self.target_occ = target_occ
        n_target = int(L * L * target_occ)
        self.m.seed_dispersed(n_target, occupancy=target_occ)
        self.m.seed_tcells(int(n_target * et_ratio))
        self.n_target = n_target


def evaluate(tonic, per_kill, theta, hill, recover_tau=10080.0, seed=0, dt=5.0,
             L=60, p_kill=0.0042):
    kw = dict(dt=dt, p_kill=p_kill, exhaust_tonic=tonic, exhaust_per_kill=per_kill,
              recover_tau=recover_tau, theta=theta, hill=hill)
    out = {}
    a = ThresholdAssay(L=L, seed=seed, **kw)
    on = lambda t: 1.0
    a.run_days(7, on);  out['d7_cont'] = a.lysis_probe()
    a.run_days(7, on);  out['d14_cont'] = a.lysis_probe()
    a.run_days(14, on); out['d28_cont'] = a.lysis_probe()
    b = ThresholdAssay(L=L, seed=seed, **kw)
    b.run_days(7, on)
    b.run_days(7, lambda t: 0.0)
    out['d14_tfi'] = b.lysis_probe()
    return out


def cell(job):
    tonic, per_kill, theta, hill = job
    path = f'{OUT}/t{tonic:.2e}_k{per_kill:.4f}_th{theta:.2f}_h{hill:.0f}.json'
    if os.path.exists(path):
        return path
    per_seed = [evaluate(tonic, per_kill, theta, hill, seed=s) for s in (0, 1)]
    meas = {k: float(np.mean([m[k] for m in per_seed])) for k in C.PHILIPP}
    row = dict(exhaust_tonic=tonic, exhaust_per_kill=per_kill, theta=theta, hill=hill,
               measured=meas, target=C.PHILIPP, loss=C.loss(meas),
               err={k: round(meas[k] - C.PHILIPP[k], 1) for k in C.PHILIPP})
    with open(path, 'w') as f:
        json.dump(row, f, indent=1)
    return path


if __name__ == '__main__':
    grid = list(itertools.product(
        [5e-5, 1e-4, 2e-4],       # exhaust_tonic
        [0.0],                    # kill-driven accrual held at 0 for this first pass
        [0.3, 0.5, 0.7],          # theta
        [4.0, 8.0]))              # hill
    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else os.cpu_count()
    print(f'{len(grid)} cells on {nproc} procs', flush=True)
    with Pool(nproc) as p:
        for i, path in enumerate(p.imap_unordered(cell, grid), 1):
            r = json.load(open(path))
            m = r['measured']
            print(f"{i}/{len(grid)} tonic={r['exhaust_tonic']:.0e} theta={r['theta']:.2f} "
                  f"hill={r['hill']:.0f} | d7={m['d7_cont']:5.1f} d14={m['d14_cont']:5.1f} "
                  f"d28={m['d28_cont']:5.1f} d14tfi={m['d14_tfi']:5.1f} | loss={r['loss']:7.1f}",
                  flush=True)
