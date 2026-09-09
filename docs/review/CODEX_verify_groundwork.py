"""Small review diagnostics. Requires numpy/scipy and --checkout at commit 5d439e8.

Runs in a temporary directory because reviewed modules create results folders on import.
Prints JSON; no full calibration or tumour experiment is launched.
"""
import argparse
import contextlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--checkout', type=Path, required=True)
args = parser.parse_args()
checkout = args.checkout.resolve()
sha = subprocess.check_output(['git', '-C', str(checkout), 'rev-parse', 'HEAD'], text=True).strip()
assert sha == '5d439e871693c3e5c3ab9f2c32ad30af5443ad9b', sha
sys.path.insert(0, str(checkout))

with tempfile.TemporaryDirectory() as scratch:
    os.chdir(scratch)
    import pk
    import schedules as S
    import exp_tarlatamab as T
    import calibrate_exhaustion as C

    result = {}
    with contextlib.redirect_stdout(io.StringIO()):
        result['pk_selftest'] = bool(pk._selftest())
    result['exposure'] = {n: S.exposure_fraction(f) for n, f in [
        ('reference_cont', S.continuous()), ('L1_cont', S.legacy_l1(0)),
        ('L1_tfi14', S.legacy_l1(14))]}

    a = C.Assay(L=20, seed=7, dt=5)
    a.m.E[a.m.T] = .6
    before = dict(t=a.m.t, nT=int(a.m.T.sum()), meanE=float(a.m.E[a.m.T].mean()),
                  history=len(a.m.history))
    old_T = a.m.T.copy()
    value = a.lysis_probe()
    after = dict(t=a.m.t, nT=int(a.m.T.sum()), meanE=float(a.m.E[a.m.T].mean()),
                 history=len(a.m.history), E_on_nonT=int(np.count_nonzero(a.m.E[~a.m.T])),
                 T_sites_changed=int(np.count_nonzero(old_T != a.m.T)))
    result['probe_mutation'] = dict(seed=7, L=20, before=before, after=after, probe_value=value)

    def auc(m, end=84):
        A = (m.alpha - m.k21) / (m.alpha - m.beta)
        B = 1 - A
        return sum(mg / m.V1 * (
            A * (-np.expm1(-m.alpha * (end-d))) / m.alpha
            + B * (-np.expm1(-m.beta * (end-d))) / m.beta)
            for d, mg in m.doses if d < end)

    result['comparators'] = {}
    for kind in ['q2w', 'infusion', 'weekly_half', 'q4w_double']:
        m = T.regimen(kind, 11.2)
        result['comparators'][kind] = dict(total_mg=m.total_mg(),
            mg_before84=sum(mg for d, mg in m.doses if d < 84), auc_0_84=auc(m))

    rows = []
    for th in [5.8, 11.2]:
        for frac in [.01, .1, .5, .9, .99]:
            m = pk.TwoCompartmentPK(th, v_ratio=pk.feasible_v_ratio_max(th)*frac)
            m.add_regimen(until_day=84)
            A = (m.alpha-m.k21)/(m.alpha-m.beta)
            B = 1-A
            # Closed-form sums over infinitely many 10-mg doses separated by 14 days.
            fast = np.exp(-m.alpha*14)
            slow = np.exp(-m.beta*14)
            peak = 10/m.V1*(A/(1-fast) + B/(1-slow))
            trough = 10/m.V1*(A*fast/(1-fast) + B*slow/(1-slow))
            ec = .105  # 1 nM under the reviewed code's 105-kDa convention
            drop = 100*(1-(trough/(trough+ec))/(peak/(peak+ec)))
            boundary = .25*trough*peak/(.75*peak-trough)/.105 if .75*peak > trough else None
            rows.append(dict(half_life=th, feasible_range_fraction=frac, v_ratio=m.v_ratio,
                day56_70_drop=T.occupancy_stats(m, ec)['drop_pct'],
                steady_state_drop_1nM=drop, EC50_threshold_nM=boundary))
    result['pk_ratio_sensitivity'] = rows
    print(json.dumps(result, indent=2))
