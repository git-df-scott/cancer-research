"""
Threshold exhaustion: a mechanism, not a rate.

WHY A NEW MECHANISM IS NEEDED, AND WHY NO RATE WILL DO
------------------------------------------------------
`lymphoid.py` makes cytotoxic function decline linearly with an exhaustion scalar:
kill_efficiency = 1 - E, with E accrued linearly in contact time. In the replenished assay that
reproduces Philipp's conditions, contact and kill rate are both approximately constant, so E is
linear in time and so is lysis. A straight line can pass through any TWO of three points. It
cannot pass through all three.

Philipp et al. (Blood 2022, PMID 35878001), continuous exposure:

    day  7   specific lysis 88.4%
    day 14   specific lysis 34.9%
    day 28   specific lysis  8.6%

Inverting lysis = naive x (1 - E) gives the exhaustion trajectory the data demand:

    naive lysis assumed    E(7)     E(14)    E(28)    E(14)/E(7)
         88.4%             0.000    0.605    0.903      infinite
         95.0%             0.069    0.633    0.909          9.11
        100.0%             0.116    0.651    0.914          5.61

**Linear accrual forces E(14)/E(7) = 2.00 exactly.** The measured ratio is 5.6 to 9.1, or
unbounded. The curve is flat, then collapses, then saturates -- a sigmoid. This is a structural
mismatch, and the calibration scan confirms it empirically rather than by argument: the cell
tonic=5e-5, per_kill=0 hits day 28 exactly on target (8.6) while giving 55.7 at day 7 against a
target of 88.4. Nailing one end throws the other.

That is why the calibration scan cannot succeed by searching harder over rates, and why the
reference model does not have this problem: it uses a hard PD-1 threshold. Cells accumulate
signal, cross a threshold, and only then lose function and begin to die. A population of cells
with dispersed accrual histories crossing a threshold at different times produces exactly the
observed sigmoid.

`README.md` and `HANDOFF.md` frame the blocking task as recalibrating exhaustion. On this evidence
that framing is incomplete in the same way `FINDINGS.md`'s was: the *rate* is not the only thing
wrong, the *functional form* is. This module supplies the alternative form so the two can be
compared on the same external curve.

WHAT IS CHANGED, AND WHAT IS NOT
--------------------------------
`Lymphoid.E` continues to be the accumulating latent signal, accrued, recovered, moved with
migrating cells, reset on influx and cleared on death by the parent class exactly as before. Only
the MAPPING from that signal to cytotoxic function changes, via the `kill_efficiency` hook. So
nothing about the parent's bookkeeping is duplicated or second-guessed.

    Lymphoid          kill_efficiency = 1 - E                    (linear)
    ThresholdLymphoid kill_efficiency = 1 / (1 + (E/theta)^n)    (Hill, sigmoid)

A Hill function is used rather than a hard step for two reasons. It is a strict generalisation --
large n approaches a hard threshold, n = 1 with theta = 1 is close to the linear default -- so the
comparison is a nested one rather than two unrelated models. And a hard step makes the fit
discontinuous in theta, which turns calibration into a search over a staircase.

FREE PARAMETERS, AND THE HONEST COST OF ADDING THEM
---------------------------------------------------
This adds two free parameters (theta, n) to the two that already exist (exhaust_tonic,
exhaust_per_kill). Four parameters against four external targets is not a comfortable margin, and
a fit that succeeds must be reported with that stated plainly. It is not evidence that the
threshold mechanism is correct; it is evidence that it is *capable* of the observed shape where
the linear form provably is not.

The discriminating test is therefore NOT the calibration fit. It is whether parameters fitted to
Philipp's assay, then transplanted unchanged into the tumour setting, reproduce the reference's
schedule ranking -- a quantity nothing here was fitted to. That is experiment R.

EXHAUSTED-CELL DEATH
--------------------
The reference also has exhausted cells die, calibrated so ~40% die within 10 days. That is a
separate mechanism from loss of function and it changes the standing population, so it is
implemented here as an option (`exhausted_death_per_min`) and reported separately. It is disabled
by default so that the functional-form change can be assessed on its own first.
"""
import numpy as np

from lymphoid import Lymphoid

# Philipp et al., Blood 2022 (PMID 35878001). Retained here so the arithmetic above is checkable.
PHILIPP_LYSIS = {7: 88.4, 14: 34.9, 28: 8.6}


def required_E(naive_lysis=95.0, lysis=PHILIPP_LYSIS):
    """Exhaustion trajectory implied by the measured lysis curve, given an assumed naive lysis."""
    return {d: 1.0 - v / naive_lysis for d, v in lysis.items()}


def linear_is_impossible(naive_lysis=95.0):
    """Ratio E(14)/E(7) demanded by the data versus the 2.00 that linear accrual forces."""
    E = required_E(naive_lysis)
    return dict(required_ratio=(np.inf if E[7] == 0 else E[14] / E[7]),
                linear_ratio=2.0, naive_lysis=naive_lysis, E=E)


class ThresholdLymphoid(Lymphoid):
    """Lymphoid with sigmoid rather than linear loss of function.

    theta : latent signal at which function is halved
    hill  : steepness. hill -> large approaches a hard threshold; hill = 1 is soft.
    exhausted_death_per_min : optional death hazard for cells past threshold (default off)
    """

    def __init__(self, *a, theta=0.35, hill=6.0, exhausted_death_per_min=0.0, **kw):
        super().__init__(*a, **kw)
        self.theta = float(theta)
        self.hill = float(hill)
        self.exhausted_death_per_min = float(exhausted_death_per_min)

    def kill_efficiency(self):
        # Clipped so the power is well behaved at E = 0 and cannot overflow at large hill.
        x = np.clip(self.E / self.theta, 0.0, 50.0)
        return 1.0 / (1.0 + x ** self.hill)

    def exhausted_mask(self):
        """T cells past the half-function point. Reported, and used for death if enabled."""
        return self.T & (self.E >= self.theta)

    def step(self, drug=None):
        super().step(drug)
        if self.exhausted_death_per_min > 0:
            m = self.exhausted_mask()
            if m.any():
                die = m & (self.rng.random(self.T.shape)
                           < self.exhausted_death_per_min * self.dt)
                self.T[die] = False
                self.E[die] = 0.0


if __name__ == '__main__':
    for naive in (88.4, 95.0, 100.0):
        r = linear_is_impossible(naive)
        print(f"naive={naive:5.1f}%  E(7)={r['E'][7]:.3f} E(14)={r['E'][14]:.3f} "
              f"E(28)={r['E'][28]:.3f}  required ratio={r['required_ratio']:.2f} "
              f"vs linear 2.00")
