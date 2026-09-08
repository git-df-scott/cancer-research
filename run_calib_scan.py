"""Parallel driver for the exhaustion calibration scan (Phase 1 gate).

Each grid cell is independent, so cells are farmed out one per process and written individually.
A stall therefore costs one cell, not the grid -- the same pattern rerun_missing.py established
after the L1 out-of-memory incident.
"""
import json, os, itertools, sys
from multiprocessing import Pool
import numpy as np
import calibrate_exhaustion as C

OUT = 'results/calib/cells'
os.makedirs(OUT, exist_ok=True)

TONICS = [2.48e-5, 5e-5, 1e-4, 2e-4]
PER_KILLS = [0.0, 0.002, 0.01, 0.05]
TAUS = [3360.0, 10080.0]
SEEDS = (0, 1)


def cell(job):
    tonic, pk, tau = job
    path = f'{OUT}/t{tonic:.3e}_k{pk:.4f}_tau{tau:.0f}.json'
    if os.path.exists(path):
        return path
    per_seed = [C.evaluate(tonic, pk, tau, seed=s) for s in SEEDS]
    meas = {k: float(np.mean([m[k] for m in per_seed])) for k in C.PHILIPP}
    row = dict(exhaust_tonic=tonic, exhaust_per_kill=pk, recover_tau=tau,
               measured=meas, target=C.PHILIPP, loss=C.loss(meas),
               err={k: round(meas[k] - C.PHILIPP[k], 1) for k in C.PHILIPP},
               seeds=list(SEEDS))
    with open(path, 'w') as f:
        json.dump(row, f, indent=1)
    return path


if __name__ == '__main__':
    jobs = list(itertools.product(TONICS, PER_KILLS, TAUS))
    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else os.cpu_count()
    print(f'{len(jobs)} cells on {nproc} procs', flush=True)
    with Pool(nproc) as p:
        for i, path in enumerate(p.imap_unordered(cell, jobs), 1):
            r = json.load(open(path))
            m = r['measured']
            print(f"{i}/{len(jobs)} tonic={r['exhaust_tonic']:.2e} k={r['exhaust_per_kill']:.4f} "
                  f"tau={r['recover_tau']:.0f} | d7={m['d7_cont']:5.1f} d14={m['d14_cont']:5.1f} "
                  f"d28={m['d28_cont']:5.1f} d14tfi={m['d14_tfi']:5.1f} | loss={r['loss']:8.1f}",
                  flush=True)
