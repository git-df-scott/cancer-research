"""
Exhaustion calibration harness: replicate Philipp's ASSAY, not the tumour.

WHY THE EXISTING CALIBRATION IS NOT WRONG, AND THE MODEL IS STILL UNDER-EXHAUSTED
--------------------------------------------------------------------------------
`lymphoid.py` sets exhaust_tonic = 2.48e-5 per minute of contact, documented as "28 d contact
-> E=1". That arithmetic is right, and it matches the external target:

    Philipp et al., Blood 2022 (PMID 35878001), continuous CD19xCD3 exposure:
        specific lysis 88.4% (day 7) -> 8.6% (day 28)
        functional retention d28/d7 = 0.097, i.e. required E(d28) ~ 0.903
    Model under 100% contact: E(d7) = 0.250, E(d28) = 1.000.

So the per-minute rate is defensible. The model nevertheless reached only 12% mean exhaustion in
L1, and a smoke run of the faithful reference configuration (1:4 E:T, 50% occupancy, no influx)
reached 4.9% -- because the model never sustains contact:

  - In L1, engaged fraction was 0.098 (follicle) to 0.463 (dispersed), so T cells accrued
    exhaustion for a fraction of the elapsed time, and influx kept adding E=0 arrivals.
  - In the faithful configuration there is no influx, but the tumour is ERADICATED before day 16
    (n0 9600 -> nB16 0). After clearance there is no contact, so exhaustion stops accruing and
    every schedule endpoint sits on the floor. The experiment cannot discriminate.

Philipp's assay does not have that problem, because it is a serial re-challenge: effectors are
kept in the presence of a maintained excess of targets. Contact is continuous BY CONSTRUCTION.

That is the mapping error. The rate was calibrated against an assay whose conditions the model
does not reproduce. This harness reproduces the assay conditions instead, so the exhaustion
parameters can be fitted to the measured curve in the setting the curve was measured in, and then
transplanted unchanged into the tumour setting.

THE SECOND STRUCTURAL DIFFERENCE
--------------------------------
`lymphoid.py` sets exhaust_per_kill = 0.0, commented "exhaustion is dwell-time driven, not
kill-count driven". The reference model does the opposite: T cells "accumulate PD-1 expression
during tumor killing under TCE exposure at higher rates than baseline" -- accrual is tied to
killing events.

This matters for the time constant, not just the magnitude. Dwell-driven exhaustion is linear in
elapsed contact time and takes weeks to bind. Kill-driven exhaustion ramps fastest exactly when
targets are abundant, which is the early phase where this model instead eradicates the tumour.
That is a candidate explanation for why the reference reaches >80% exhaustion by day 16 and this
model does not, and it is testable here rather than assumable.

Both parameters are therefore fitted jointly against the same external curve.

WHAT IS FITTED, AND TO WHAT
---------------------------
Free:   exhaust_tonic (per contact-minute), exhaust_per_kill (per kill event)
Fixed:  everything else, including p_kill

Targets, all external and all from Philipp et al. (PMID 35878001):
    T1  specific lysis at day  7, continuous exposure          88.4%
    T2  specific lysis at day 28, continuous exposure           8.6%
    T3  specific lysis at day 14, continuous exposure          34.9%
    T4  specific lysis at day 14, with a treatment-free interval 93.4%

T3 and T4 together pin the RECOVERY behaviour (recover_tau), which T1/T2 alone cannot. A fit that
matches T1/T2 but not T3/T4 has calibrated decay without calibrating reinvigoration, and would
silently break every treatment-free-interval comparison downstream. All four are reported.

HARD RULE
---------
Calibration targets the measured lysis curve, which is fixed and external. It NEVER targets a
schedule ranking. Any fit is accepted or rejected on T1-T4 alone. If a fit that matches T1-T4
still fails to reproduce the reference's schedule ranking, that is a result to report, not a
reason to refit.
"""
import json, os, itertools
import numpy as np

from lymphoid import Lymphoid

OUT = 'results/calib'
os.makedirs(OUT, exist_ok=True)

MIN_PER_DAY = 1440

# Philipp et al., Blood 2022, PMID 35878001. Specific lysis, percent.
PHILIPP = dict(d7_cont=88.4, d14_cont=34.9, d28_cont=8.6, d14_tfi=93.4)


class Assay:
    """In-vitro replica: fixed effectors, targets replenished to hold density constant.

    This is the serial-re-challenge condition Philipp's numbers were measured under. Targets are
    topped back up every step, so contact never lapses because the effectors won the fight. That
    is the whole point: it isolates exhaustion from tumour-clearance dynamics.
    """

    def __init__(self, L=60, seed=0, target_occ=0.40, et_ratio=0.25, **kw):
        self.m = Lymphoid(L=L, seed=seed, t_influx=0.0, t_div=0.0,
                          p_div=0.0, p_death=0.0, **kw)
        self.L = L
        self.target_occ = target_occ
        n_target = int(L * L * target_occ)
        self.m.seed_dispersed(n_target, occupancy=target_occ)
        self.m.seed_tcells(int(n_target * et_ratio))
        self.n_target = n_target

    def _replenish(self):
        """Hold target number at its set point by refilling empty sites at random."""
        deficit = self.n_target - int(self.m.B.sum())
        if deficit <= 0:
            return
        free = np.flatnonzero(~self.m.B & ~self.m.T)
        if len(free) == 0:
            return
        pick = self.m.rng.choice(free, min(deficit, len(free)), replace=False)
        flat = self.m.B.reshape(-1)
        flat[pick] = True

    def lysis_probe(self, probe_min=240):
        """Specific lysis over a short probe window, as a percentage of the naive-effector rate.

        Philipp measures lysis in a fresh short co-culture against a reference. Here the probe is
        the kill count over `probe_min` minutes at the current exhaustion state, normalised by the
        same probe run at E=0 with an identical configuration.
        """
        saved_E = self.m.E.copy()
        saved_kills = self.m.kills
        # exhausted-state rate
        k0 = self.m.kills
        for _ in range(int(probe_min / self.m.dt)):
            self._replenish()
            self.m.step(1.0)
        rate_now = self.m.kills - k0
        # naive reference rate, same cells, E zeroed
        self.m.E = np.zeros_like(self.m.E)
        k1 = self.m.kills
        for _ in range(int(probe_min / self.m.dt)):
            self._replenish()
            self.m.step(1.0)
        rate_naive = self.m.kills - k1
        # restore
        self.m.E = saved_E
        self.m.kills = saved_kills
        return 100.0 * rate_now / rate_naive if rate_naive else 0.0

    def run_days(self, days, drug_fn):
        steps = int(days * MIN_PER_DAY / self.m.dt)
        for i in range(steps):
            self._replenish()
            self.m.step(drug_fn(i * self.m.dt))


def evaluate(exhaust_tonic, exhaust_per_kill, recover_tau=10080.0,
             seed=0, dt=5.0, L=60, p_kill=0.0042):
    """Run the four Philipp conditions and return measured specific lysis for each."""
    kw = dict(dt=dt, p_kill=p_kill, exhaust_tonic=exhaust_tonic,
              exhaust_per_kill=exhaust_per_kill, recover_tau=recover_tau)
    out = {}

    # Continuous exposure, probed at days 7, 14, 28
    a = Assay(L=L, seed=seed, **kw)
    on = lambda t: 1.0
    a.run_days(7, on);  out['d7_cont'] = a.lysis_probe()
    a.run_days(7, on);  out['d14_cont'] = a.lysis_probe()
    a.run_days(14, on); out['d28_cont'] = a.lysis_probe()

    # Treatment-free interval arm, probed at day 14.
    # Philipp's comparison is a 7-day exposure followed by a 7-day drug-free window.
    b = Assay(L=L, seed=seed, **kw)
    b.run_days(7, on)
    b.run_days(7, lambda t: 0.0)
    out['d14_tfi'] = b.lysis_probe()
    return out


def loss(meas):
    """Sum of squared error against Philipp, in percentage points. Reported per-target too."""
    return sum((meas[k] - PHILIPP[k]) ** 2 for k in PHILIPP)


def scan(tonics, per_kills, recover_taus=(10080.0,), seeds=(0,), **kw):
    """Coarse grid over the two free exhaustion parameters. Writes every cell, keeps every seed."""
    rows = []
    for tonic, pk, tau in itertools.product(tonics, per_kills, recover_taus):
        per_seed = [evaluate(tonic, pk, tau, seed=s, **kw) for s in seeds]
        meas = {k: float(np.mean([m[k] for m in per_seed])) for k in PHILIPP}
        row = dict(exhaust_tonic=tonic, exhaust_per_kill=pk, recover_tau=tau,
                   measured=meas, target=PHILIPP, loss=loss(meas),
                   err={k: round(meas[k] - PHILIPP[k], 1) for k in PHILIPP})
        rows.append(row)
        print(f"tonic={tonic:.2e} per_kill={pk:.4f} tau={tau:.0f} "
              f"d7={meas['d7_cont']:5.1f} d14={meas['d14_cont']:5.1f} "
              f"d28={meas['d28_cont']:5.1f} d14tfi={meas['d14_tfi']:5.1f} "
              f"loss={row['loss']:8.1f}", flush=True)
    with open(f'{OUT}/scan.json', 'w') as f:
        json.dump(rows, f, indent=1)
    return rows


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == 'baseline':
        # What the shipped parameters actually produce in the assay condition.
        m = evaluate(2.48e-5, 0.0)
        print('shipped parameters in Philipp assay conditions:')
        for k in PHILIPP:
            print(f'  {k:9s} model {m[k]:6.1f}   Philipp {PHILIPP[k]:6.1f}   '
                  f'err {m[k]-PHILIPP[k]:+6.1f}')
        print(f'  loss = {loss(m):.1f}')
    else:
        scan(tonics=[1e-5, 2.48e-5, 5e-5, 1e-4],
             per_kills=[0.0, 0.002, 0.01, 0.05],
             recover_taus=[3360.0, 10080.0])
