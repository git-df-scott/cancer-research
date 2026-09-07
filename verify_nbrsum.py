"""Equivalence check for the fast nbr_sum against the original scipy implementation.

The Phase E positive control was launched before this optimisation and therefore ran on the
scipy version; every later phase runs on the fast version. That is only acceptable if the two
are equivalent on the code paths that carry results, so this checks both the legacy single-state
path and the new two-state path, on both architectures, and is kept in the repository as evidence.
"""
import json, sys
import numpy as np
import lymphoid
from lymphoid import Lymphoid
import exp_schedule as S


def snapshot(use_reference, arch, seed, twostate):
    fast = lymphoid.nbr_sum
    if use_reference:
        lymphoid.nbr_sum = lymphoid._nbr_sum_reference
    try:
        kw = dict(L=S.L, seed=seed, dt=S.DT, p_kill=S.PARAMS['p_kill'],
                  p_div=S.PARAMS['p_div'], t_div=0.0, t_influx=2.5e-5)
        if twostate:
            kw.update(exhaust_model='twostate', k_exh=1.104e-4, frac_durable=0.93,
                      recover_tau_r=10780.0, c50_exh=0.0)
        m = Lymphoid(**kw)
        if arch == 'follicle':
            m.seed_follicles([(S.L // 2, S.L // 2)], S.R_BIG)
        else:
            m.seed_dispersed(S.N_B, occupancy=0.95)
        m.seed_tcells(S.N_T)
        for i in range(288 * 5):                      # 5 simulated days
            m.step(1.0 if (i // 288) % 2 == 0 else 0.0)
        nT = int(m.T.sum())
        return dict(nB=m.nB, nT=nT, kills=int(m.kills),
                    E=round(float(m.E[m.T].mean()) if nT else 0.0, 12),
                    Ed=round(float(m.Ed[m.T].mean()) if nT else 0.0, 12))
    finally:
        lymphoid.nbr_sum = fast


if __name__ == '__main__':
    allok = True
    for twostate in (False, True):
        for arch in ('follicle', 'dispersed'):
            a = snapshot(True, arch, 3, twostate)
            b = snapshot(False, arch, 3, twostate)
            ok = a == b
            allok &= ok
            print(f'  {"twostate" if twostate else "legacy  "}  {arch:10s}  '
                  f'{"IDENTICAL" if ok else "DIFFERS"}')
            if not ok:
                print(f'      scipy {a}\n      fast  {b}')
    print('\nALL IDENTICAL' if allok else '\nMISMATCH - do not use the fast path')
    sys.exit(0 if allok else 1)
