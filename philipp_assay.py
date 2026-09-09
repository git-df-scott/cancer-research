"""
Philipp assay replica: chronic culture and readout are SEPARATE experiments.

See PHILIPP_ASSAY.md for the protocol this implements and its sources. The short version is that
the previous rig in `calibrate_exhaustion.py` conflated two distinct experiments and measured the
wrong quantity:

  chronic stimulation   E:T 1:4, irradiated targets, 28 d in 4 weekly cycles, drug per arm,
                        medium/targets/drug replenished on day 3, T cells re-isolated and
                        recultured with fresh targets and drug on day 7
  cytotoxicity readout  E:T 1:1, FRESH targets, 72 h, AMG 562 at 5 ng/mL, normalised against a
                        control construct

Exhaustion accrues in the first condition and is read out in the second. The old probe ran a 4-hour
kill-rate measurement inside the chronic culture itself and normalised against the same cells with
exhaustion zeroed. Every part of that is a different quantity from Equation 1 of the paper:

    % specific lysis = (1 - target count WITH engager / target count WITH control construct) x 100

The denominator is a NO-ENGAGER control run with the same (exhausted) T cells, not a naive-effector
reference. So the observable is target depletion attributable to the drug, over 72 hours, at 1:1.

DESIGN RULE, NON-NEGOTIABLE
---------------------------
The readout must not mutate the culture it measures. That is enforced structurally here rather than
by save/restore: the readout never touches the culture object at all. It reads out a vector of
per-cell exhaustion values and instantiates entirely fresh simulations from them. The previous bug
-- restoring a lattice-indexed exhaustion field onto a population that had migrated, silently
rejuvenating every cell -- is not merely fixed, it is unrepresentable in this design.

WHAT IS ASSUMED AND MARKED AS SUCH
----------------------------------
- The paper does not state whether the day-3 replenishment tops up the existing well or replaces
  its contents. Implemented as a top-up to the set point. Recorded in PHILIPP_ASSAY.md.
- Absolute seeding densities are not given. Lattice occupancy is chosen so the readout is not
  crowding-limited, and `target_occ` is exposed so the choice can be swept.
- Irradiated targets do not divide: p_div = 0. This one IS stated.
"""
import copy
import numpy as np

from lymphoid import Lymphoid

MIN_PER_DAY = 1440

# Philipp et al., Blood 2022 (PMID 35878001). All four are 72 h specific lysis, n=6, mean +/- SEM.
PHILIPP = dict(d7_cont=88.4, d14_cont=34.9, d28_cont=8.6, d14_tfi=93.4)
PHILIPP_FIG = dict(d7_cont='Fig 2E', d14_cont='Fig 3E', d28_cont='Fig 2E', d14_tfi='Fig 3E')


class ChronicCulture:
    """28-day chronic stimulation at E:T 1:4 with irradiated (non-dividing) targets."""

    def __init__(self, L=60, seed=0, target_occ=0.40, et_ratio=0.25, model_cls=Lymphoid, **kw):
        kw.setdefault('p_div', 0.0)      # targets are irradiated - stated in the methods
        kw.setdefault('p_death', 0.0)
        self.m = model_cls(L=L, seed=seed, t_influx=0.0, t_div=0.0, **kw)
        self.L = L
        self.n_target = int(L * L * target_occ)
        self.m.seed_dispersed(self.n_target, occupancy=target_occ)
        self.m.seed_tcells(int(self.n_target * et_ratio))
        self.model_kw = dict(kw)
        self.model_cls = model_cls

    def _replenish_targets(self):
        """Top up targets to the set point. Whether the paper tops up or replaces is UNAVAILABLE."""
        deficit = self.n_target - int(self.m.B.sum())
        if deficit <= 0:
            return
        free = np.flatnonzero(~self.m.B & ~self.m.T)
        if not len(free):
            return
        pick = self.m.rng.choice(free, min(deficit, len(free)), replace=False)
        self.m.B.reshape(-1)[pick] = True

    def _run(self, days, drug):
        for _ in range(int(days * MIN_PER_DAY / self.m.dt)):
            self.m.step(drug)

    def run_cycle(self, drug_on):
        """One 7-day stimulation cycle: replenish at day 3, reculture with fresh targets at day 7."""
        d = 1.0 if drug_on else 0.0
        self._run(3, d)
        self._replenish_targets()          # day 3: medium, targets and drug replenished
        self._run(4, d)
        self._replenish_targets()          # day 7: recultured with fresh targets
        return self

    def harvest(self):
        """Per-cell exhaustion of the living T cells. A vector of cell states, not a lattice field.

        This is what leaves the culture. Nothing about the culture is modified, and the readout
        cannot write back, which is why the old restore-a-moved-field bug cannot recur.
        """
        return self.m.E[self.m.T].copy()


def readout_specific_lysis(chronic, assay_seed, hours=72, L=60, et_ratio=1.0,
                           target_occ=0.20, model_cls=None, n_plate=250, **kw):
    """Equation 1: 72 h, E:T 1:1, fresh targets, normalised to a no-engager control construct.

    Both arms use the SAME harvested T cells and the same assay RNG seed, so the comparison is
    paired and the ratio is not contaminated by between-arm sampling noise.
    """
    harvested = chronic.harvest()
    if len(harvested) == 0:
        return 0.0

    # STANDARDISE THE PLATING DENSITY. Philipp plate a defined number of effectors at E:T 1:1;
    # they do not plate however many happened to survive the chronic culture. Reading out at the
    # surviving count confounds per-cell function with population attrition, and the confound is
    # large: at a fixed exhaustion of 0.10 this readout gives 100% lysis with 600 cells plated and
    # 51.7% with 60. The chronic culture falls from 229 to 137 to 51 T cells across the three
    # timepoints, so each readout was being taken at a different density.
    #
    # That artefact masked the recovery mechanism entirely. The day-14 TFI arm reaches mean
    # exhaustion 0.048 against the continuous arm's 0.318 -- reinvigoration works in the model --
    # yet read out at 65% because only 134 cells were plated. It also made recover_tau look
    # computationally dead: a 14-fold change moved the reading by 1.5 points.
    #
    # Sampling with replacement when the culture has fewer than n_plate survivors preserves the
    # exhaustion DISTRIBUTION, which is the quantity the assay is meant to interrogate, while
    # holding density fixed. Recorded as an assumption: the real experiment has a finite cell
    # yield, and a culture that cannot supply n_plate cells is being represented optimistically.
    rng = np.random.default_rng(assay_seed)
    n_T = int(n_plate)
    E_vals = rng.choice(harvested, size=n_T, replace=(len(harvested) < n_T))
    model_cls = model_cls or chronic.model_cls
    params = dict(chronic.model_kw)
    params.update(kw)
    params.setdefault('p_div', 0.0)
    params.setdefault('p_death', 0.0)
    n_target = int(n_T / et_ratio)

    def arm(drug):
        m = model_cls(L=L, seed=assay_seed, t_influx=0.0, t_div=0.0, **params)
        m.seed_dispersed(n_target, occupancy=target_occ)
        idx = np.flatnonzero(~m.B)
        pick = m.rng.choice(idx, min(n_T, len(idx)), replace=False)
        flat_T = m.T.reshape(-1); flat_E = m.E.reshape(-1)
        flat_T[pick] = True
        flat_E[pick] = E_vals[:len(pick)]        # harvested cells carry their exhaustion
        n0 = int(m.B.sum())
        for _ in range(int(hours * 60 / m.dt)):
            m.step(drug)
        return int(m.B.sum()), n0

    n_with, _ = arm(1.0)
    n_ctrl, _ = arm(0.0)                          # control construct: no engager
    return 100.0 * (1.0 - n_with / n_ctrl) if n_ctrl else 0.0


def evaluate(model_cls=Lymphoid, assay_seed=12345, seed=0, dt=5.0, L=60, **params):
    """The four Philipp conditions. Continuous arm and TFI arm, probed at the published timepoints.

    TFI arm: drug absent during cycles 2 (days 7-14) and 4 (days 21-28), targets still present.
    """
    kw = dict(dt=dt, **params)
    out = {}

    cont = ChronicCulture(L=L, seed=seed, model_cls=model_cls, **kw)
    cont.run_cycle(True)                                    # days 0-7
    out['d7_cont'] = readout_specific_lysis(cont, assay_seed, L=L, model_cls=model_cls)
    cont.run_cycle(True)                                    # days 7-14
    out['d14_cont'] = readout_specific_lysis(cont, assay_seed, L=L, model_cls=model_cls)
    cont.run_cycle(True); cont.run_cycle(True)              # days 14-28
    out['d28_cont'] = readout_specific_lysis(cont, assay_seed, L=L, model_cls=model_cls)

    tfi = ChronicCulture(L=L, seed=seed, model_cls=model_cls, **kw)
    tfi.run_cycle(True)                                     # cycle 1: drug on
    tfi.run_cycle(False)                                    # cycle 2: drug off, targets present
    out['d14_tfi'] = readout_specific_lysis(tfi, assay_seed, L=L, model_cls=model_cls)
    return out


def loss(meas, targets=PHILIPP):
    return sum((meas[k] - targets[k]) ** 2 for k in targets)


# --------------------------------------------------------------------- Phase C invariants
def _state_fingerprint(c):
    """Everything persistent. If a readout changes any of this, the readout is destructive."""
    m = c.m
    return (m.t, int(m.kills), m.B.tobytes(), m.T.tobytes(), m.E.tobytes(),
            m.hits.tobytes(), m.drug, float(m.cum_engaged_min),
            repr(m.rng.bit_generator.state))


def _invariants():
    ok = True

    def check(name, cond):
        nonlocal ok
        ok &= bool(cond)
        print(f'  [{"PASS" if cond else "FAIL"}] {name}')

    c = ChronicCulture(L=40, seed=3, dt=5.0, p_kill=0.0042, exhaust_tonic=5e-5)
    c.run_cycle(True)

    # 1. The readout must leave the culture bit-identical in every persistent variable,
    #    including the RNG state. This is the regression test for the destructive-probe bug.
    before = _state_fingerprint(c)
    y1 = readout_specific_lysis(c, assay_seed=999, L=40)
    after = _state_fingerprint(c)
    check('readout leaves culture bit-identical (incl. RNG state)', before == after)

    # 2. Same state and same assay seed must give the same reading.
    y2 = readout_specific_lysis(c, assay_seed=999, L=40)
    check(f'readout is deterministic given assay seed ({y1:.4f} == {y2:.4f})', y1 == y2)

    # 3. Exhaustion must actually leave the culture. If harvest returned zeros the readout would
    #    silently measure naive cells -- which is exactly what the old bug did.
    Ev = c.harvest()
    check(f'harvest carries non-zero exhaustion (mean {Ev.mean():.4f})', Ev.mean() > 0)

    # 4. Specific lysis is bounded and, for a lightly exhausted population, substantial.
    check(f'0 <= lysis <= 100 (got {y1:.1f})', 0.0 <= y1 <= 100.0)

    # 5. A fully exhausted population must read out near zero lysis. Guards against the readout
    #    being insensitive to the very quantity it is supposed to measure.
    c_dead = ChronicCulture(L=40, seed=3, dt=5.0, p_kill=0.0042, exhaust_tonic=5e-5)
    c_dead.m.E[c_dead.m.T] = 1.0
    y_dead = readout_specific_lysis(c_dead, assay_seed=999, L=40)
    check(f'fully exhausted population reads ~0 lysis (got {y_dead:.2f})', y_dead < 5.0)

    # 6. DENSITY INVARIANCE. The readout must depend on per-cell function, not on how many cells
    #    the chronic culture happened to have left. This is the regression test for the
    #    plating-density confound (see superseded/calib2_density_bug/README.md), which passed
    #    invariants 1-5 while measuring the wrong quantity entirely: at fixed exhaustion 0.10,
    #    plating 600 cells gave 100% lysis and plating 60 gave 51.7%.
    #
    #    Checks 1-5 verify that measurement does not MUTATE state. This one verifies that it
    #    measures the intended QUANTITY. They are different classes of check and the first does
    #    not imply the second.
    class _Fixed:
        def __init__(self, E, kw, cls): self._E = E; self.model_kw = kw; self.model_cls = cls
        def harvest(self): return self._E

    kwd = dict(dt=5.0, p_kill=1e-3, exhaust_tonic=5e-5)
    readings = [readout_specific_lysis(_Fixed(np.full(n, 0.10), kwd, Lymphoid),
                                       assay_seed=4242, L=40)
                for n in (60, 120, 300, 600)]
    spread = max(readings) - min(readings)
    check(f'readout invariant to surviving cell count at fixed E '
          f'(spread {spread:.1f} pts over 60-600 cells: {[round(r,1) for r in readings]})',
          spread < 5.0)

    # 7. The readout must still RESPOND to exhaustion. A readout made density-invariant by being
    #    insensitive to everything would pass check 6 and be useless.
    lo = readout_specific_lysis(_Fixed(np.full(300, 0.05), kwd, Lymphoid), assay_seed=4242, L=40)
    hi = readout_specific_lysis(_Fixed(np.full(300, 0.95), kwd, Lymphoid), assay_seed=4242, L=40)
    check(f'readout still responds to exhaustion (E=0.05 -> {lo:.1f}, E=0.95 -> {hi:.1f})',
          lo - hi > 20.0)

    print('\nPhase C invariants:', 'PASS' if ok else 'FAIL')
    return ok


if __name__ == '__main__':
    import sys
    if not _invariants():
        sys.exit(1)
