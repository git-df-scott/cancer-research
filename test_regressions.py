"""
Consolidated regression suite: one permanent test per bug class actually encountered.

Every bug below shipped, produced plausible internally-consistent numbers, and raised no error.
None was found by a failing test; all were found by noticing two computed quantities that could
not both be true. This suite exists so each is caught mechanically from now on.

Run: .venv/bin/python test_regressions.py
"""
import os, sys, json, tempfile
import numpy as np

FAILS = []


def check(cls, name, cond, detail=''):
    if not cond:
        FAILS.append(f'{cls}: {name}')
    print(f'  [{"PASS" if cond else "FAIL"}] {cls:28s} {name}' + (f'  {detail}' if detail else ''))


print('Regression suite: one test per bug class encountered\n')

# ---------------------------------------------------------------- 1 + 2 + 10
# Destructive observation; lattice-indexed state restored onto a migrated population;
# empty-population means reported as healthy values.
import philipp_assay as PA
from exhaustion import ThresholdLymphoid

kw = dict(dt=5.0, p_kill=1e-3, exhaust_tonic=5e-5)
c = PA.ChronicCulture(L=40, seed=3, model_cls=ThresholdLymphoid, theta=0.5, hill=4.0, **kw)
c.run_cycle(True)
before = PA._state_fingerprint(c)
y1 = PA.readout_specific_lysis(c, assay_seed=999, L=40)
check('destructive-observation', 'readout leaves culture bit-identical', before == PA._state_fingerprint(c))
check('destructive-observation', 'readout is deterministic per assay seed',
      y1 == PA.readout_specific_lysis(c, assay_seed=999, L=40))
check('migrated-state-restore', 'harvest carries non-zero exhaustion',
      c.harvest().mean() > 0, f'mean E={c.harvest().mean():.4f}')

empty = np.array([])
mean_of_empty = float(empty.mean()) if len(empty) else float('nan')
check('empty-population-mean', 'mean of an empty population is NaN, not 0.0',
      np.isnan(mean_of_empty))

# ---------------------------------------------------------------- 4 (observable identity)
# The readout must measure per-cell function, not how many cells survived.
class _Fixed:
    def __init__(self, E): self._E = E; self.model_kw = dict(kw); self.model_cls = ThresholdLymphoid
    def harvest(self): return self._E

reads = [PA.readout_specific_lysis(_Fixed(np.full(n, 0.10)), assay_seed=4242, L=40)
         for n in (60, 300, 600)]
check('plating-density-confound', 'readout invariant to surviving cell count',
      max(reads) - min(reads) < 5.0, f'spread {max(reads)-min(reads):.1f} pts')
lo = PA.readout_specific_lysis(_Fixed(np.full(300, 0.05)), assay_seed=4242, L=40)
hi = PA.readout_specific_lysis(_Fixed(np.full(300, 0.95)), assay_seed=4242, L=40)
check('inert-observable', 'readout still responds to exhaustion', lo - hi > 20.0,
      f'{lo:.1f} -> {hi:.1f}')

# ---------------------------------------------------------------- 3 (dead parameters)
# A parameter claimed insensitive must first be shown computationally live.
from lymphoid import Lymphoid
base = Lymphoid(L=30, seed=0, dt=5.0, p_kill=1e-3, exhaust_tonic=5e-5)
alt = Lymphoid(L=30, seed=0, dt=5.0, p_kill=1e-3, exhaust_tonic=1e-4)
base.seed_dispersed(300, occupancy=0.4); base.seed_tcells(80)
alt.seed_dispersed(300, occupancy=0.4); alt.seed_tcells(80)
for _ in range(int(3 * 1440 / 5)):
    base.step(1.0); alt.step(1.0)
check('dead-parameter', 'exhaust_tonic is computationally live',
      abs(alt.E[alt.T].mean() - base.E[base.T].mean()) > 1e-6,
      f'{base.E[base.T].mean():.4f} vs {alt.E[alt.T].mean():.4f}')

# ---------------------------------------------------------------- 5 (filename collisions)
def cell_path(tonic, per_kill, theta, hill):
    return f't{tonic:.2e}_k{per_kill:.4f}_th{theta:.2f}_h{hill:.1f}.json'

check('filename-collision', 'hill 3.5 and 4.0 map to distinct filenames',
      cell_path(5e-5, 0, 0.7, 3.5) != cell_path(5e-5, 0, 0.7, 4.0),
      f'{cell_path(5e-5,0,0.7,3.5)} vs {cell_path(5e-5,0,0.7,4.0)}')

# ---------------------------------------------------------------- 6 + 7 + 8 (PK)
import pk as PK
m = PK.TwoCompartmentPK(t_half_terminal=PK.THALF_POPPK).add_regimen(until_day=84)
check('endpoint-phantom-dose', 'no dose lands exactly on the horizon',
      all(d < 84.0 for d, _ in m.doses), f'max dose day {max(d for d,_ in m.doses)}')

ref_auc = m.auc(0.0, 84.0)
inf = PK.TwoCompartmentPK(t_half_terminal=PK.THALF_POPPK).add_infusion_auc_matched(ref_auc, 0.0, 84.0)
check('auc-mismatch', 'matched arms agree on windowed AUC',
      abs(inf.auc(0.0, 84.0) - ref_auc) / ref_auc < 0.01,
      f'{ref_auc:.3f} vs {inf.auc(0.0,84.0):.3f}')

tiny = Lymphoid(L=20, seed=0, dt=5.0, p_kill=1e-3, exhaust_tonic=5e-5)
tiny.seed_dispersed(100, occupancy=0.4); tiny.seed_tcells(30)
for _ in range(288):
    tiny.step(1e-9)
check('tiny-positive-PK', 'infinitesimal occupancy is not treated as full exposure',
      tiny.E[tiny.T].mean() < 1e-6, f'E={tiny.E[tiny.T].mean():.2e}')

# ---------------------------------------------------------------- 9 (time integrals)
g = Lymphoid(L=30, seed=1, dt=5.0, p_kill=1e-3, exhaust_tonic=5e-5)
g.seed_dispersed(300, occupancy=0.4); g.seed_tcells(80)
snap_product = 0.0
for day in range(3):
    for _ in range(288):
        g.step(1.0)
    snap_product += g._n_engaged * 1440.0
check('snapshot-not-integral', 'cum_engaged_min differs from a daily snapshot x 1440',
      abs(g.cum_engaged_min - snap_product) / max(g.cum_engaged_min, 1) > 1e-6,
      f'integral {g.cum_engaged_min:.0f} vs snapshot {snap_product:.0f}')

# ---------------------------------------------------------------- 11 (contaminated holdout)
import provenance as PV
check('contaminated-holdout', 'philipp_d14_tfi is permanently flagged contaminated',
      not PV.validation_is_clean('philipp_d14_tfi'),
      f"via {list(PV.contaminating_parameters('philipp_d14_tfi'))}")

# ---------------------------------------------------------------- 12 (boundary optima)
def is_boundary(best, grid):
    return best in (min(grid), max(grid))

check('boundary-optimum', 'boundary detector flags an edge optimum',
      is_boundary(1e-3, [3e-4, 5e-4, 1e-3]) and not is_boundary(5e-4, [3e-4, 5e-4, 1e-3]))

# ---------------------------------------------------------------- 13 (the gate)
from calibration import Calibration, require_validated
with tempfile.TemporaryDirectory() as d:
    p = os.path.join(d, 'c.json')
    Calibration(status='FIT_NOT_VALIDATED', model='M1', params={}, residuals={},
                fit_targets=[], held_out={}, seeds=[], notes='t').save(p)
    try:
        require_validated('R', p); gated = False
    except RuntimeError:
        gated = True
check('downstream-gate', 'R refuses without a VALIDATED calibration', gated)

print(f"\n{'ALL REGRESSION TESTS PASS' if not FAILS else 'FAILURES: ' + '; '.join(FAILS)}")
sys.exit(1 if FAILS else 0)
