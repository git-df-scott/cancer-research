"""
Schedule definitions shared by the positive control and the architecture comparison.
Pre-registered in PHASE2_PREREG.md section 4, before any Phase 2 result was generated.

Family A reproduces experiment L1's design: a single 28-day cycle with the last k days off.
Family B is the repeating family in which the reference effect actually lives. Obertopp/Basanta
(bioRxiv 2025, PPR1121269) report that a 7-day interval loses its advantage by day 42 but that
"shorter TFIs consistently outperformed both 7-day and continuous schedules", and that a
"Monday-through-Friday" regimen achieved comparable benefits. MO_FR is 5 days on, 2 days off,
every week. L1 never tested any repeating schedule, so it never sampled that family at all.
"""
CYCLE = 28

FAMILY_A = [('A_cont', 28, 0), ('A_tfi2', 26, 2), ('A_tfi4', 24, 4),
            ('A_tfi7', 21, 7), ('A_tfi14', 14, 14)]
FAMILY_B = [('B_6on1off', 6, 1), ('B_MO_FR', 5, 2), ('B_4on3off', 4, 3),
            ('B_12on2off', 12, 2), ('B_7on7off', 7, 7)]
ALL = FAMILY_A + FAMILY_B
CONTINUOUS = 'A_cont'


def by_name(name):
    for n, on, off in ALL:
        if n == name:
            return on, off
    raise KeyError(name)


def schedule_fn(on_days, off_days):
    """Repeating cycle: engager present for on_days, absent for off_days."""
    if off_days == 0:
        return lambda s: 1.0
    on_min = on_days * 1440.0
    cyc_min = (on_days + off_days) * 1440.0
    return lambda s: 1.0 if (s.t % cyc_min) < on_min else 0.0


def duty_cycle(on_days, off_days):
    """Fraction of time the engager is present. A schedule with duty cycle d is dose-matched by a
    continuous arm at level d, because the kill hazard is linear in drug."""
    return on_days / float(on_days + off_days)
