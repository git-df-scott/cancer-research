"""
Stage 2 of Phase 2 calibration: verify that the parameters fitted in calib_philipp.py reproduce
the measured curve when run through the ACTUAL lymphoid.py code path, not just the scalar ODE.

The ODE fit assumes a T cell is in antigen contact 100% of the time, which is true of Philipp's
assay (E:T 1:4 with targets replenished) but is a property of the assay, not of the model. This
file rebuilds the assay on the lattice - irradiated targets that neither divide nor die,
replenished on schedule, a closed T-cell pool with no influx, division or death - and checks that
the engaged fraction really is ~1 and that the model then lands on the measured points.

If the lattice replica disagrees with the ODE fit, the ODE fit is not a calibration of the model.
"""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lymphoid import Lymphoid, nbr_sum
from calib_philipp import DATA, TFI_OFF, drug_on, predict_one

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'results')
L = 48
N_B = 1200            # ~52% occupancy: every T cell has a target neighbour with prob >0.99
N_T = N_B // 4        # E:T = 1:4, as in the paper
DT = 5.0
REPLENISH_EVERY = 720.0   # minutes. Paper replenishes targets on day 3 and re-cultures every 7 d;
                          # replenishing more often only makes contact more continuous, which is
                          # the condition being replicated.


def run_arm(arm, seed, params, days=28.0):
    m = Lymphoid(L=L, seed=seed, dt=DT,
                 p_div=0.0, p_death=0.0,        # irradiated OCI-Ly1: no division, no turnover
                 t_influx=0.0, t_div=0.0, t_death=0.0,   # closed culture
                 p_kill=0.0012, exhaust_model='twostate', **params)
    m.seed_dispersed(N_B, occupancy=0.95)
    m.seed_tcells(N_T)
    nsteps = int(days * 1440 / DT)
    out, engfrac = {}, []
    for i in range(nsteps):
        day = i * DT / 1440.0
        m.step(1.0 if drug_on(arm, day) else 0.0)
        engfrac.append(m.history[-1]['engFrac'])
        if (i + 1) % int(REPLENISH_EVERY / DT) == 0:
            _replenish(m)
        d = (i + 1) * DT / 1440.0
        if abs(d - round(d)) < 1e-9 and int(round(d)) in (7, 14, 28):
            out[int(round(d))] = dict(
                f=float(1.0 - m.E[m.T].mean()),
                Er=float(m.Er[m.T].mean()), Ed=float(m.Ed[m.T].mean()),
                C_days=float(m.C[m.T].mean() / 1440.0))
    return out, float(np.mean(engfrac))


def _replenish(m):
    """Restore the target pool to its starting size, as the protocol does with fresh cells."""
    need = N_B - m.nB
    if need <= 0:
        return
    free = np.flatnonzero((~m.B & ~m.T).ravel())
    if not len(free):
        return
    pick = m.rng.choice(free, min(need, len(free)), replace=False)
    flat = m.B.ravel().copy(); flat[pick] = True
    m.B[:] = flat.reshape(m.L, m.L)


if __name__ == '__main__':
    cal = json.load(open(f'{OUT}/calib_philipp.json'))
    reps = cal['representatives']
    allres = {}
    for label, ch in reps.items():
        a = ch['a']
        params = dict(k_exh=ch['k_per_min'], frac_durable=ch['rho'],
                      recover_tau_r=ch['tau_r_min'], c50_exh=ch['c50_min'], n_exh=4.0)
        print(f'=== [{label}]  k={ch["k_per_day"]:.4f}/d  rho={ch["rho"]:.2f}  '
              f'tau_r={ch["tau_r_days"]:.2f}d  C50={ch["c50_days"]:.0f}d  a={a:.2f} ===')
        seeds = [0, 1, 2]
        res = {}
        for arm in ('cont', 'tfi'):
            per = [run_arm(arm, s_, params) for s_ in seeds]
            res[f'{arm}_engfrac'] = float(np.mean([p[1] for p in per]))
            for d in (7, 14, 28):
                if d not in per[0][0]:
                    continue
                res[f'{arm}_{d}'] = dict(
                    f=float(np.mean([p[0][d]['f'] for p in per])),
                    Er=float(np.mean([p[0][d]['Er'] for p in per])),
                    Ed=float(np.mean([p[0][d]['Ed'] for p in per])),
                    C_days=float(np.mean([p[0][d]['C_days'] for p in per])))
        print(f'    mean engaged fraction over the run:  cont {res["cont_engfrac"]:.3f}   '
              f'tfi {res["tfi_engfrac"]:.3f}   (assay condition requires ~1)')
        print('    readout      arm       day  measured   ODE fit   ABM')
        sse_abm = 0.0
        for kind, arm, d, y in DATA:
            if kind == 'lysis':
                f_abm = res[f'{arm}_{d}']['f']
                v = 100.0 * (1 - np.exp(-a * f_abm))
                tag = 'lysis %'
            else:
                v = 100.0 * res[f'cont_{d}']['f'] / max(res[f'tfi_{d}']['f'], 1e-12)
                tag = 'GzB ratio'
            ode = predict_one(ch)[0][(kind, arm, d)]
            sse_abm += (v - y) ** 2
            print(f'    {tag:10s} {arm:9s} d{d:<3d} {y:7.1f}   {ode:7.1f}   {v:7.1f}')
        print(f'    ABM SSE = {sse_abm:.1f}   (ODE fit SSE = {ch["sse"]:.1f})')
        print('    per-cell state at day 28, continuous arm:  '
              f'Er={res["cont_28"]["Er"]:.3f}  Ed={res["cont_28"]["Ed"]:.3f}  '
              f'contact={res["cont_28"]["C_days"]:.1f} d\n')
        allres[label] = dict(res=res, sse_abm=float(sse_abm), sse_ode=ch['sse'])
    json.dump(allres, open(f'{OUT}/calib_abm_check.json', 'w'), indent=1)
    print(f'written {OUT}/calib_abm_check.json')
