"""
Dosing-schedule generators.

WHY THIS FILE EXISTS
--------------------
`exp_schedule.py` (experiment L1) encodes a treatment-free interval as ONE interruption at the
end of a 28-day cycle:

    on for (28 - tfi) days, then off for tfi days

The reference this project compares itself against does something different. Obertopp, Froid,
Pilon-Thomas & Basanta (bioRxiv 2025.11.17.688873, published PMC12667981) state in their methods:

    "Intermittent schedules TFI_2 through TFI_7 incorporate recurring treatment-free intervals
     (2-7 days, respectively), alternating with periods of TCE dosing"

and separately that the 42-day cycle is 28 days of dosing followed by a uniform 14-day rest
applied to EVERY arm, continuous included.

Those are not the same regimen, and the difference is not cosmetic. Under the L1 encoding TFI_k
is exactly "continuous minus the last k days of drug", so it is dominated by continuous by
construction unless exhaustion is severe enough to repay the lost exposure. That single fact
predicts the monotone ordering L1 actually reported (TFI0 > TFI2 > TFI4 > TFI7 > TFI14) in all
three architectures, with no appeal to exhaustion calibration at all. The reference's TFI_k is a
DUTY CYCLE: many short recovery windows spread through treatment, not one long one at the end.

L1 additionally dosed its continuous arm for all 42 days, where the reference doses for 28 and
rests for 14. So the L1 continuous arm received 1.5x the drug-time the reference gave it.

This module provides both encodings so the difference can be measured rather than argued about.

WHAT IS ASSUMED, AND WHY IT IS SWEPT
------------------------------------
The reference does not tabulate how many days are ON between its OFF-intervals. That number is
required to build the duty cycle and it is not recoverable from the text. It is therefore a free
assumption, and `exp_replicate.py` sweeps it rather than picking one. If the ranking is stable
across the sweep the assumption does not matter; if it flips, the sweep locates where.
"""

MIN_PER_DAY = 1440


def continuous(dose_days=28, total_days=42):
    """Reference CONT arm: drug on for dose_days, then a rest phase to total_days.

    Note this is NOT what L1's tfi=0 arm did. L1 dosed continuously for the entire 42 days.
    Use `legacy_l1(0)` for that.
    """
    on_min = dose_days * MIN_PER_DAY
    return lambda s: 1.0 if s.t < on_min else 0.0


def duty_cycle(on_days, off_days, dose_days=28, total_days=42):
    """Reference TFI arm: repeating on_days ON / off_days OFF, through the dosing phase only.

    After `dose_days` every arm rests, so the arms differ only in what happened during dosing.
    """
    on_min = on_days * MIN_PER_DAY
    period_min = (on_days + off_days) * MIN_PER_DAY
    dose_min = dose_days * MIN_PER_DAY

    def sched(s):
        if s.t >= dose_min:
            return 0.0
        return 1.0 if (s.t % period_min) < on_min else 0.0
    return sched


def legacy_l1(tfi, cycle_days=28):
    """The encoding experiment L1 actually used. Preserved verbatim so the ablation can toggle it.

    on for (cycle_days - tfi) days, then off for tfi days, repeating with no rest phase.
    """
    if tfi == 0:
        return lambda s: 1.0
    on_min = (cycle_days - tfi) * MIN_PER_DAY
    cyc_min = cycle_days * MIN_PER_DAY
    return lambda s: 1.0 if (s.t % cyc_min) < on_min else 0.0


def exposure_fraction(sched, total_days=42, dt=5.0):
    """Fraction of the horizon with drug present. Dose-matching controls need this.

    Uses a stub with only the attribute the schedules read, so it never touches a real model.
    """
    class _Stub:
        t = 0
    stub = _Stub()
    steps = int(total_days * MIN_PER_DAY / dt)
    on = 0
    for i in range(steps):
        stub.t = i * dt
        if sched(stub) > 0:
            on += 1
    return on / steps
