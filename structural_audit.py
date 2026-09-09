"""
Phase 5: attack M1's limiting behaviour. A structural counterexample kills the model class
regardless of how well it fits the calibration.

A model can reproduce three numbers and still be nonsense outside them. These tests ask whether
M1 behaves sensibly in regimes the calibration never probes -- which is exactly where R and T
would operate.
"""
import numpy as np
import philipp_assay as PA
from exhaustion import ThresholdLymphoid
from lymphoid import Lymphoid

# The current M1 point. Ranges are being re-selected in parallel; these limits should hold for
# any admissible parameterisation, so they are not contingent on that outcome.
M1 = dict(p_kill=1e-3, exhaust_tonic=5e-5, exhaust_per_kill=0.0,
          recover_tau=10080.0, theta=0.5, hill=4.0)

FAILS = []


def check(name, cond, detail=''):
    tag = 'PASS' if cond else 'FAIL'
    if not cond:
        FAILS.append(name)
    print(f'  [{tag}] {name}' + (f'  {detail}' if detail else ''))


def culture(days, drug, antigen=True, dt=5.0, L=40, **over):
    kw = dict(M1); kw.update(over); kw['dt'] = dt
    c = PA.ChronicCulture(L=L, seed=0, model_cls=ThresholdLymphoid, **kw)
    if not antigen:
        c.m.B[:] = False
        c.n_target = 0
    for _ in range(int(days * 1440 / dt)):
        if antigen:
            c._replenish_targets()
        c.m.step(drug)
    E = c.m.E[c.m.T]
    # Return NaN rather than 0.0 on extinction. Returning 0.0 let the "very long stimulation stays
    # bounded" test PASS for the wrong reason: E(56d) read 0.000 not because exhaustion is bounded
    # but because every T cell had died. A test that passes on an empty population is worse than
    # no test.
    return float(E.mean()) if len(E) else float('nan')


print('Phase 5 structural audit of M1\n')

# 1-4. Exposure limits. Exhaustion must be monotone in occupancy and vanish at zero.
e0 = culture(7, 0.0)
eps = culture(7, 1e-6)
half = culture(7, 0.5)
sat = culture(7, 1.0)
check('zero drug -> no drug-attributable exhaustion', e0 < 1e-6, f'E={e0:.2e}')
check('infinitesimal occupancy -> negligible exhaustion', eps < 1e-3, f'E(1e-6)={eps:.2e}')
check('exhaustion monotone in occupancy', e0 <= eps <= half <= sat,
      f'{e0:.3f} <= {eps:.3f} <= {half:.3f} <= {sat:.3f}')
check('half occupancy is intermediate, not saturating', 0.2 * sat < half < 0.9 * sat,
      f'half={half:.3f} sat={sat:.3f}')

# 5. No antigen. Engager present but nothing to engage: exhaustion must not accrue.
noag = culture(7, 1.0, antigen=False)
check('drug without antigen -> no exhaustion', noag < 1e-6, f'E={noag:.2e}')

# 6. Withdrawal with antigen still present must permit recovery.
kw = dict(M1); kw['dt'] = 5.0
c = PA.ChronicCulture(L=40, seed=0, model_cls=ThresholdLymphoid, **kw)
for _ in range(int(7 * 1440 / 5)):
    c._replenish_targets(); c.m.step(1.0)
e_on = float(c.m.E[c.m.T].mean())
for _ in range(int(7 * 1440 / 5)):
    c._replenish_targets(); c.m.step(0.0)
e_off = float(c.m.E[c.m.T].mean())
check('withdrawal with antigen present permits recovery', e_off < e_on,
      f'{e_on:.3f} -> {e_off:.3f}')
check('sustained engagement does NOT permit recovery', culture(14, 1.0) > culture(7, 1.0),
      f'd7={culture(7,1.0):.3f} d14={culture(14,1.0):.3f}')

# 7-8. Long limits. Bounded above, and recovery approaches zero.
long_stim = culture(56, 1.0)
check('very long stimulation: population survives to be measured', not np.isnan(long_stim),
      f'E(56d)={long_stim:.3f} (NaN means the T-cell population went extinct)')
check('very long stimulation stays bounded in [0,1]',
      (not np.isnan(long_stim)) and 0.0 <= long_stim <= 1.0, f'E(56d)={long_stim:.3f}')
c2 = PA.ChronicCulture(L=40, seed=0, model_cls=ThresholdLymphoid, **kw)
for _ in range(int(14 * 1440 / 5)):
    c2._replenish_targets(); c2.m.step(1.0)
for _ in range(int(56 * 1440 / 5)):
    c2._replenish_targets(); c2.m.step(0.0)
e_long_rec = float(c2.m.E[c2.m.T].mean())
check('very long recovery: population survives to be measured', not np.isnan(e_long_rec),
      f'E={e_long_rec:.4f}')
check('very long recovery drives exhaustion toward zero',
      (not np.isnan(e_long_rec)) and e_long_rec < 0.02, f'E after 56d rest={e_long_rec:.4f}')

# 12. STRUCTURAL COUNTEREXAMPLE CHECK. Philipp's chronic culture sustains and proliferates T cells;
#     this model's only loses them. The mismatch is large and, critically, ARM-DEPENDENT in the
#     real experiment, which is what makes it a counterexample rather than a nuisance.
c3 = PA.ChronicCulture(L=60, seed=0, model_cls=ThresholdLymphoid, **dict(M1, dt=5.0))
n_start = int(c3.m.T.sum())
for _ in range(4):
    c3.run_cycle(True)
frac = int(c3.m.T.sum()) / n_start
check('chronic culture retains a usable T-cell population to day 28 (>50%)', frac > 0.5,
      f'{100*frac:.0f}% of starting cells remain; Philipp report proliferation over this window '
      f'(CD2+ fold change d14: 1.1 continuous vs 4.1 TFI)')

# 9. Timestep dependence. dt is a numerical choice, not biology.
e_dt5 = culture(7, 1.0, dt=5.0)
e_dt1 = culture(7, 1.0, dt=1.0)
rel = abs(e_dt5 - e_dt1) / max(e_dt1, 1e-9)
check('exhaustion not strongly dt-dependent (<15%)', rel < 0.15,
      f'dt=5 {e_dt5:.3f} vs dt=1 {e_dt1:.3f}, rel {100*rel:.1f}%')

# 10. Observable bounds.
class _F:
    def __init__(s, E): s._E = E; s.model_kw = dict(M1, dt=5.0); s.model_cls = ThresholdLymphoid
    def harvest(s): return s._E
ys = [PA.readout_specific_lysis(_F(np.full(200, e)), assay_seed=7, L=40)
      for e in (0.0, 0.5, 1.0)]
check('lysis within physical bounds [0,100]', all(0.0 <= y <= 100.0 for y in ys),
      f'{[round(y,1) for y in ys]}')
check('lysis monotone decreasing in exhaustion', ys[0] >= ys[1] >= ys[2],
      f'{[round(y,1) for y in ys]}')

# 11. Pathological parameter sensitivity: a 1% change must not swing the observable wildly.
base = culture(7, 1.0)
pert = culture(7, 1.0, exhaust_tonic=M1['exhaust_tonic'] * 1.01)
swing = abs(pert - base) / max(base, 1e-9)
check('no pathological sensitivity to 1% parameter change (<5%)', swing < 0.05,
      f'{base:.4f} -> {pert:.4f}, {100*swing:.2f}%')

print(f"\nStructural audit: {'PASS' if not FAILS else 'FAIL -> ' + ', '.join(FAILS)}")
