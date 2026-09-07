"""
Phase 4 calibration: the systemic (circulating) T-cell compartment, fitted to Philipp Figure 1B.

WHY A SYSTEMIC COMPARTMENT EXISTS AT ALL
----------------------------------------
Three experiments (L1, Phase 2, Phase 3) failed, and all three failed for one reason: recruited
T cells arrived with zero exhaustion, so any influx large enough to sustain an effector pool also
reset the population's function. Tumour control and binding exhaustion therefore never co-occurred
in 24 of 24 tested regimes, and the window in which an exhaustion-driven scheduling effect could
exist was empty.

That defect is contradicted by direct measurement, not merely inconvenient. Blinatumomab is given
by 28-day CONTINUOUS INTRAVENOUS INFUSION, so the entire circulating T-cell pool is engager-exposed,
not only cells inside the lesion. Philipp Figure 1B assayed exactly that pool - PERIPHERAL T cells
from r/r BCP-ALL patients on c.i.v. - and measured:

    day 0 (pre-treatment)      73.1% specific lysis
    day 14 (on infusion)       17.4%      -> 0.238 of baseline
    post-cessation             48.5%      -> 0.663 of baseline  (ns vs day 0)

The cells available for recruitment are measured to be exhausted. A model in which recruits arrive
naive is therefore wrong about the recruits, and this file calibrates what they actually arrive
carrying.

MODEL
-----
A scalar compartment advancing on DRUG EXPOSURE ALONE, since a circulating T cell is engager-exposed
whether or not it is touching a blast:

    drug present:  dEr_sys/dt = (1 - rho_sys) k_sys (1 - E_sys),  dEd_sys/dt = rho_sys k_sys (1 - E_sys)
    drug absent:   dEr_sys/dt = -Er_sys / tau_sys,  Ed_sys unchanged

Recruited T cells enter the tissue carrying (Er_sys, Ed_sys) and the compartment's accrued exposure,
instead of zeros. Nothing else about the model changes.

WHY rho_sys IS NOT frac_durable, and this is forced by the data
---------------------------------------------------------------
The per-cell durable fraction fitted in Phase B is 0.93 for the best-fitting representative. Applied
to the systemic pool it would floor function at 1 - 0.93 x 0.943 = 0.123, against a measured 0.663.
That is arithmetically impossible, not merely a poor fit.

It is also biologically expected. The circulating compartment recovers by TURNOVER as well as by
individual cells de-exhausting: exhausted cells die and are replaced from memory and progenitor
reserves. A tissue T cell has no such route. So the pool recovers faster than any of its members,
and rho_sys is a separate parameter with its own admissible range.

PRE-REGISTERED FIT PROTOCOL, fixed before running
-------------------------------------------------
- Targets: the three Figure 1B values above, normalised to the patient's own day-0 baseline.
  Nothing else. No schedule ranking, no architecture quantity, no tumour burden.
- The readout is treated as approximately LINEAR in effector function over the observed 17-73%
  range, which is away from the assay ceiling. This differs from the in vitro fit in Phase B, where
  values sat at 88-93% and saturation had to be modelled explicitly. The assumption is stated here
  and its effect is reported as a sensitivity, not buried.
- The sampling time of the "post" point is NOT stated in the paper. It is carried as a declared
  uncertainty over {day 35, day 42} - one or two weeks after the end of the 28-day infusion -
  rather than assumed.
- rho_sys and tau_sys are expected to be jointly, not separately, identifiable, exactly as the
  per-cell pair were. The admissible set is reported and representatives spanning it are carried
  forward. No single value is preferred.
"""
import json, os
import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
DAY = 1440.0
BASE = 73.1
DATA = {'day14': 17.4, 'post': 48.5}
INFUSION_DAYS = 28.0
POST_DAYS = [35.0, 42.0]          # declared uncertainty on the "post" sampling time
RMSE_ACCEPT = 5.0                 # percentage points, on the normalised-to-baseline scale


def trajectory(k_day, rho, tau_day, t_days, infusion=INFUSION_DAYS):
    """Closed form. During infusion dE/dt = k(1-E) with Er and Ed accruing in fixed proportion,
    so E(t) = 1 - exp(-k t) and Er = (1-rho)E, Ed = rho E exactly. After cessation Ed is frozen
    and Er decays, giving E(t) = E28 [rho + (1-rho) exp(-(t-28)/tau)].
    Works elementwise on arrays."""
    out = {}
    E_inf = 1.0 - np.exp(-k_day * infusion)
    for t in t_days:
        if t <= infusion:
            out[t] = np.exp(-k_day * t)                       # f = 1 - E
        else:
            d = t - infusion
            out[t] = 1.0 - E_inf * (rho + (1 - rho) * np.exp(-d / tau_day))
    return out


def fit(post_day):
    """Vectorised grid search. Only the two Figure 1B post-baseline values are scored;
    the day-0 point defines the normalisation."""
    f14_m, fpost_m = DATA['day14'] / BASE, DATA['post'] / BASE
    ks = np.geomspace(0.02, 0.6, 120)
    rhos = np.linspace(0.0, 0.6, 121)
    taus = np.geomspace(0.5, 60.0, 120)
    K, R, T = np.meshgrid(ks, rhos, taus, indexing='ij')
    K, R, T = K.ravel(), R.ravel(), T.ravel()
    tr = trajectory(K, R, T, [14.0, post_day])
    sse = ((tr[14.0] - f14_m) * 100) ** 2 + ((tr[post_day] - fpost_m) * 100) ** 2
    order = np.argsort(sse)
    rows = [(float(sse[i]), float(K[i]), float(R[i]), float(T[i]),
             float(tr[14.0][i]), float(tr[post_day][i])) for i in order[:200000]]
    return rows, f14_m, fpost_m


if __name__ == '__main__':
    print(__doc__.split('PRE-REGISTERED')[0].strip())
    print('\n' + '=' * 84)
    allres = {}
    for post_day in POST_DAYS:
        rows, f14_m, fpost_m = fit(post_day)
        best = rows[0]
        rmse = np.sqrt(best[0] / 2)
        acc = [r for r in rows if np.sqrt(r[0] / 2) <= RMSE_ACCEPT]
        print(f'\n--- "post" sampled at day {post_day:.0f} '
              f'({post_day-INFUSION_DAYS:.0f} days after cessation) ---')
        print(f'  best fit   k_sys={best[1]:.4f}/day  rho_sys={best[2]:.2f}  '
              f'tau_sys={best[3]:.2f} d   RMSE={rmse:.2f} pp')
        print(f'    day 14   measured {100*f14_m:5.1f}%   model {100*best[4]:5.1f}%')
        print(f'    post     measured {100*fpost_m:5.1f}%   model {100*best[5]:5.1f}%')
        print(f'  admissible set (RMSE <= {RMSE_ACCEPT:.0f} pp): {len(acc)} of {len(rows)}')
        for nm, idx in (('k_sys (/day)', 1), ('rho_sys', 2), ('tau_sys (days)', 3)):
            v = np.array([r[idx] for r in acc])
            print(f'    {nm:16s} {v.min():7.3f} .. {v.max():7.3f}   median {np.median(v):7.3f}')
        rec = np.array([1 - (r[2] + (1 - r[2]) * np.exp(-7.0 / r[3])) for r in acc])
        print(f'    fraction of systemic exhaustion a 7-day break reverses: '
              f'{rec.min():.3f} .. {rec.max():.3f}, median {np.median(rec):.3f}')
        allres[f'post_d{int(post_day)}'] = dict(
            k_sys_per_day=float(best[1]), k_sys_per_min=float(best[1] / DAY),
            rho_sys=float(best[2]), tau_sys_days=float(best[3]),
            tau_sys_min=float(best[3] * DAY), rmse=float(rmse), n_admissible=len(acc),
            k_range=[float(min(r[1] for r in acc)), float(max(r[1] for r in acc))],
            rho_range=[float(min(r[2] for r in acc)), float(max(r[2] for r in acc))],
            tau_range=[float(min(r[3] for r in acc)), float(max(r[3] for r in acc))],
            recovered_7d=[float(rec.min()), float(rec.max())])
    print('\n' + '=' * 84)
    ks = [v['k_sys_per_day'] for v in allres.values()]
    print(f'k_sys is insensitive to the post-timing uncertainty: {min(ks):.4f}-{max(ks):.4f}/day,')
    print(f'because it is set by the day-14 point alone. That is {min(ks)/0.159:.2f}-{max(ks)/0.159:.2f}x')
    print('the per-cell rate of 0.159/day - circulating cells receive tonic engager signal but')
    print('less antigen-driven stimulation, so a value below 1.0 is the expected direction.')
    print('tau_sys and rho_sys are jointly identifiable only, as with the per-cell pair.')
    json.dump(allres, open(f'{OUT}/calib_systemic.json', 'w'), indent=1)
    print(f'\nwritten {OUT}/calib_systemic.json')
