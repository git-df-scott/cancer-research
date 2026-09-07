"""Locate the parameter regime where continuous engager therapy CONTROLS BUT DOES NOT CLEAR
a follicular tumour, which is the clinically relevant regime (CR rates 40-60%, not 100%).
T cells are a finite supply arriving from the vasculature (border influx), they do not
proliferate freely, and they exhaust with dwell time in contact."""
import sys, json, itertools
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, '.')
from lymphoid import Lymphoid, nbr_sum

DAYS = 21


def run(args):
    pk, inf, pdiv, nT0 = args
    T = Lymphoid(L=130, seed=0, p_kill=pk, t_influx=inf, p_div=pdiv, t_div=0.0)
    T.seed_follicles([(65, 65)], 42)
    n0 = T.nB
    T.seed_tcells(nT0)
    traj = []
    for d in range(DAYS):
        for _ in range(1440):
            T.step(1.0)
        nT = int(T.T.sum())
        traj.append((d + 1, T.nB, nT, round(T.history[-1]['meanE'], 3)))
        if T.nB == 0:
            break
    eng = ((nbr_sum(T.B) > 0) & T.T).sum() / max(int(T.T.sum()), 1)
    return dict(p_kill=pk, influx=inf, p_div=pdiv, nT0=nT0, n0=n0, traj=traj,
                cleared=T.nB == 0, engFrac=float(eng))


if __name__ == '__main__':
    jobs = [(pk, inf, pdiv, nT0)
            for pk in (0.0042, 0.0012)
            for inf in (0.0004, 0.0001)
            for pdiv in (1/2880, 1/1440)
            for nT0 in (200,)]
    with Pool(6) as p:
        res = p.map(run, jobs)
    print(f"{'p_kill':>8}{'influx':>8}{'p_div':>9}{'nT0':>5} | {'nB d7':>7}{'nB d14':>8}{'nB d21':>8}{'nT d21':>8}{'E d21':>7}{'eng':>6}  {'cleared':>7}")
    for r in res:
        tr = {x[0]: x for x in r['traj']}
        g = lambda d, i: (tr[d][i] if d in tr else (0 if i == 1 else -1))
        print(f"{r['p_kill']:>8.4f}{r['influx']:>8.4f}{r['p_div']:>9.5f}{r['nT0']:>5} | "
              f"{g(7,1):>7}{g(14,1):>8}{g(21,1):>8}{g(21,2):>8}{g(21,3):>7}{r['engFrac']:>6.2f}  {str(r['cleared']):>7}   n0={r['n0']}")
    json.dump(res, open('results/scan.json', 'w'))
