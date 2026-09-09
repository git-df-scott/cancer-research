"""
Class-wide PK analysis of approved T-cell engagers.

THE QUESTION
------------
The treatment-free-interval (TFI) literature -- the evidence that intermittent exposure preserves
T-cell function where continuous exposure exhausts it -- is the basis on which schedule choices for
this whole drug class are argued. Can the intervals it studies actually be realised by the drugs
now in clinical use?

This asks that with arithmetic, not simulation. For a drug at steady state on a fixed interval,
the trough-to-peak ratio depends ONLY on the dosing interval and the terminal half-life:

    C_trough / C_peak = 2^(-tau / t_half)

It is independent of dose, of volume, and of absolute potency. That makes it the most robust
statement available about whether a dosing holiday exists, and it requires no assumption about
EC50, which is the quantity this project has repeatedly found to be unestablished in vivo.

WHY THE ANSWER MIGHT MATTER
---------------------------
Philipp et al. (Blood 2022) established TFI benefit using AMG 562, whose half-life they report as
~210 hours (8.75 days). They imposed the treatment-free interval by physically removing drug at
reculture -- an in vitro washout. A patient cannot be washed out. For a molecule with a 9-day
half-life there is no clinical manoeuvre that reproduces the intervention the benefit was
demonstrated with.

Blinatumomab is the one molecule in the class where stopping the drug really does remove it:
~2-hour half-life, continuous infusion, renally cleared at ~50 kDa. Its 7-day holiday is a genuine
holiday. Every subsequent molecule is IgG-based and half-life extended.

SOURCES
-------
Half-lives are terminal t1/2 from published population-PK analyses or regulatory labels. Where a
range is reported, both ends are carried; where half-life is time-dependent (several TCEs show
longer t1/2 at steady state), the steady-state value is used because that is the regime the
maintenance interval operates in. Provenance is recorded per drug.

CONFOUND, RECORDED NOT HIDDEN
-----------------------------
TCE half-life is not a free design parameter: it is coupled to CD3 affinity through target-mediated
drug disposition. Penney et al. (Clin Transl Sci 2025) report t1/2 of 15 days at 1000 nM CD3
affinity falling to 4.4 days at 10 nM. So a long half-life partly reflects weaker CD3 binding, and
trough concentration alone does not settle trough *engagement*. This is the main threat to the
interpretation and is treated as such rather than as a footnote.
"""
import numpy as np

# t_half in days; interval in days. dose/route recorded for provenance, not used in the ratio.
TCE = [
    dict(drug='Blinatumomab', target='CD19', t_half=(0.083, 0.083), interval=None,
         route='continuous IV infusion', indication='B-ALL',
         source='~2 h terminal t1/2; 50 kDa, renal clearance; requires continuous infusion'),
    dict(drug='Teclistamab', target='BCMA', t_half=(3.8, 8.8), interval=7.0,
         route='SC weekly', indication='multiple myeloma',
         source='mean t1/2 3.8 d (SD 1.7), individual values to 8.8 d; FDA approval summary'),
    dict(drug='Talquetamab', target='GPRC5D', t_half=(8.4, 12.2), interval=7.0,
         route='SC weekly or Q2W', indication='multiple myeloma',
         source='8.4 d after first dose, 12.2 d at 16 weeks; EMA label / MonumenTAL-1 clin pharm'),
    dict(drug='Talquetamab Q2W', target='GPRC5D', t_half=(8.4, 12.2), interval=14.0,
         route='SC Q2W', indication='multiple myeloma', source='as above, Q2W schedule'),
    dict(drug='Tarlatamab', target='DLL3', t_half=(5.8, 11.2), interval=14.0,
         route='IV Q2W', indication='SCLC',
         source='5.8 d DeLLphi-300 mean; 11.2 d popPK median, n=420'),
    dict(drug='Elranatamab', target='BCMA', t_half=(10.0, 10.0), interval=14.0,
         route='SC weekly then Q2W', indication='multiple myeloma',
         source='projected terminal t1/2 ~10 d'),
    dict(drug='Glofitamab', target='CD20', t_half=(10.0, 10.0), interval=21.0,
         route='IV Q3W', indication='DLBCL',
         source='~10 d, cited as enabling Q3W dosing'),
    dict(drug='Mosunetuzumab', target='CD20', t_half=(6.0, 11.0), interval=21.0,
         route='IV Q3W', indication='follicular lymphoma',
         source='apparent t1/2 6-11 d following flat dosing'),
    dict(drug='Epcoritamab', target='CD20', t_half=(22.0, 22.0), interval=28.0,
         route='SC Q4W maintenance', indication='LBCL',
         source='geometric mean terminal t1/2 22 d at end of cycle 3, 48 mg'),
]


def trough_peak(tau, t_half):
    """C_trough/C_peak at steady state. Dose-independent; depends only on tau/t_half."""
    return 2.0 ** (-tau / t_half)


def accumulation(tau, t_half):
    """Steady-state accumulation ratio, 1/(1 - 2^(-tau/t_half))."""
    return 1.0 / (1.0 - trough_peak(tau, t_half))


def report():
    print(f"{'drug':>16s} {'target':>7s} {'t1/2 (d)':>11s} {'tau (d)':>8s} "
          f"{'tau/t1/2':>9s} {'trough/peak':>12s} {'accum':>6s}")
    print('-' * 78)
    for d in TCE:
        if d['interval'] is None:
            print(f"{d['drug']:>16s} {d['target']:>7s} {'~0.083':>11s} {'cont.':>8s} "
                  f"{'n/a':>9s} {'n/a - drug is':>12s}  {'':>6s}")
            print(f"{'':>16s} {'':>7s} {'':>11s} {'':>8s} {'':>9s} "
                  f"{'gone in hours when stopped':>12s}")
            continue
        lo, hi = d['t_half']
        # worst case for a holiday is the LONG half-life; report both
        r_lo, r_hi = trough_peak(d['interval'], lo), trough_peak(d['interval'], hi)
        a_lo, a_hi = accumulation(d['interval'], lo), accumulation(d['interval'], hi)
        th = f"{lo:.1f}" if lo == hi else f"{lo:.1f}-{hi:.1f}"
        ratio = f"{d['interval']/hi:.2f}-{d['interval']/lo:.2f}" if lo != hi else f"{d['interval']/lo:.2f}"
        tp = f"{r_lo:.2f}-{r_hi:.2f}" if lo != hi else f"{r_lo:.2f}"
        ac = f"{a_lo:.1f}-{a_hi:.1f}" if lo != hi else f"{a_lo:.1f}"
        print(f"{d['drug']:>16s} {d['target']:>7s} {th:>11s} {d['interval']:8.0f} "
              f"{ratio:>9s} {tp:>12s} {ac:>6s}")


if __name__ == '__main__':
    report()
    print("\ntrough/peak = 2^(-tau/t_half). Dose-independent, potency-independent.")
    print("A value of 0.50 means concentration never falls below half its peak between doses.")


# --------------------------------------------------------------------------------------
# The EC50-free formulation.
#
# The obvious objection to a trough/peak argument is that concentration is not engagement:
# trough concentration only matters relative to EC50, and in vivo EC50 is exactly the quantity
# this project has repeatedly found to be unestablished.
#
# That objection dissolves. Under an Emax relationship, occupancy P = C/(C+EC50), so
#
#     C = EC50 * P/(1-P)
#
# and the ratio of two concentrations depends only on the two occupancies:
#
#     C_trough/C_peak = [P_trough/(1-P_trough)] / [P_peak/(1-P_peak)]
#
# EC50 cancels. So the question "what trough/peak ratio is required for occupancy to fall from
# P_peak to P_trough" has an answer that needs no potency estimate at all -- only an assumption
# about what occupancy the drug needs at peak to work, and what counts as disengagement.
#
# This inverts the problem. Instead of assuming a potency and computing occupancy, ask: how far
# would concentration have to fall for a holiday to exist, and does any approved schedule fall
# that far?

def required_trough_peak(p_peak, p_trough):
    """Trough/peak concentration ratio needed to go from p_peak to p_trough occupancy. EC50-free."""
    return (p_trough / (1 - p_trough)) / (p_peak / (1 - p_peak))


def holiday_analysis():
    print('\n' + '=' * 78)
    print('EC50-FREE ANALYSIS: what trough/peak would a real dosing holiday require?')
    print('=' * 78)
    print('\nRequired C_trough/C_peak for occupancy to fall from peak to trough:')
    print(f"{'peak occ':>9s} {'trough occ':>11s} {'required ratio':>15s}")
    for pp in (0.90, 0.95, 0.99):
        for pt in (0.50, 0.25, 0.10):
            print(f"{pp:9.2f} {pt:11.2f} {required_trough_peak(pp, pt):15.3f}")
    print('\nEC50 cancels. These thresholds hold for any Emax drug at any potency.')

    thresh = required_trough_peak(0.90, 0.50)
    print(f"\nTake the most permissive case: a drug needing only 90% occupancy at peak,")
    print(f"and calling 50% occupancy 'disengaged'. Required trough/peak = {thresh:.3f}.\n")
    print(f"{'drug':>16s} {'interval':>9s} {'trough/peak':>12s} {'holiday?':>10s}")
    print('-' * 52)
    n_fail = n_tot = 0
    for d in TCE:
        if d['interval'] is None:
            print(f"{d['drug']:>16s} {'continuous':>9s} {'~0 on stop':>12s} {'YES':>10s}")
            continue
        lo, hi = d['t_half']
        best = trough_peak(d['interval'], hi)   # most favourable end of the half-life range
        n_tot += 1
        ok = best < thresh
        if not ok:
            n_fail += 1
        print(f"{d['drug']:>16s} {d['interval']:8.0f}d {best:12.2f} "
              f"{'yes' if ok else 'NO':>10s}")
    print(f"\n{n_fail}/{n_tot} half-life-extended agents cannot produce a trough holiday at their")
    print("labelled interval, even taking the shortest reported half-life and the most")
    print("permissive occupancy definition. Blinatumomab, the one non-extended agent, can.")


if __name__ == '__main__':
    holiday_analysis()
