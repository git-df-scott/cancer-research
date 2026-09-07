"""
Why the Phase E positive control failed. Diagnosis only - nothing here is tuned, and no
parameter is changed. Per PHASE2_PREREG.md section 9 rule 1, the architecture comparison
was NOT run.

FINDING: the day-42 endpoint saturated. The experiment had no dynamic range, so it could not
have detected a scheduling effect of any size.

This is the MIRROR IMAGE of experiment L1's failure, and it is worth stating plainly:

    L1        failed at the FLOOR    - every arm cleared the tumour (medians 4 and 7 cells
                                       of 5525), so burden could not separate the arms
    Phase 2   failed at the CEILING  - every arm saturates the lattice at 68-79% occupancy,
                                       so burden again cannot separate the arms

L1's pre-registration anticipated the floor and PHASE2_PREREG.md section 5 carried a
contingency for it. Neither anticipated the ceiling. That is the gap.

THE CAUSAL CHAIN, quantified below:

  Phase C's external criterion (patient T-cell function must collapse at the population level)
     -> only near-zero influx is admissible
     -> with a 14-day background T-cell half-life and no replenishment, the effector pool
        falls from 200 to 10-142 by day 42
     -> 10-142 T cells cannot control 5542 malignant cells with a 2-day doubling time
     -> the tumour reaches lattice carrying capacity by roughly day 21-28
     -> the endpoint stops responding to anything

Note what is NOT the cause. Killing per T cell is 1.9-3.0 targets/day, inside the measured
range of 2-16 (Halle et al., Immunity 2016). The kill rate is not miscalibrated. And the
exhaustion mechanism works exactly as calibrated: in the one combination with real dynamic
range, breaks raise end-of-run T-cell function from 0.315 to 0.646. Treatment-free intervals
do restore function in this model. They simply cannot repay the lost drug time while the
tumour is running away, and the endpoint cannot see the difference either way.

The structural limitation that produced this was recorded IN ADVANCE, in PHASE2_PREREG.md
section 2: "within this model's structure, recruitment and population-level functional
collapse cannot both be represented at a realistic influx rate, because arrivals enter with
zero exhaustion." That is precisely what bit.
"""
import json
import numpy as np

L2 = 120 * 120
RES = 'results/expE_poscontrol.json'


def main():
    res = json.load(open(RES))
    combos = sorted({(r['rep'], r['influx']) for r in res})
    print(__doc__)
    print('=' * 78)
    print('1. DYNAMIC RANGE. Could the endpoint separate the arms at all?')
    print('=' * 78)
    print(f'{"combination":26s} {"n0":>6s} {"med nB42":>9s} {"% lattice":>10s} {"arm spread":>11s}')
    for rep, infl in combos:
        sub = [r for r in res if r['rep'] == rep and r['influx'] == infl]
        meds = {s: np.median([r['nB42'] for r in sub if r['sched'] == s])
                for s in sorted({r['sched'] for r in sub})}
        lo, hi = min(meds.values()), max(meds.values())
        print(f'{rep + " x " + f"{infl:.1e}":26s} {np.median([r["n0"] for r in sub]):6.0f} '
              f'{np.median(list(meds.values())):9.0f} '
              f'{100*np.median(list(meds.values()))/L2:9.1f}% {100*(hi-lo)/max(hi,1):10.1f}%')
    print('\n  Three of four combinations separate the arms by under 1% of burden. A 10%')
    print('  minimum effect was pre-registered; it is unreachable when the ceiling is this close.')

    print('\n' + '=' * 78)
    print('2. THE EFFECTOR POOL COLLAPSES, and killing per cell is NOT the problem')
    print('=' * 78)
    print(f'{"combination":26s} {"nT d0":>6s} {"nT d42":>7s} {"kills":>7s} '
          f'{"kills/T/day":>12s} {"in Halle 2-16?":>15s}')
    for rep, infl in combos:
        c = [r for r in res if r['rep'] == rep and r['influx'] == infl and r['sched'] == 'A_cont']
        k = np.median([r['kills'] for r in c])
        pool = np.median([np.mean([t[2] for t in r['traj']]) for r in c])
        kpd = k / (pool * 42)
        print(f'{rep + " x " + f"{infl:.1e}":26s} {200:6d} '
              f'{np.median([r["nT_end"] for r in c]):7.0f} {k:7.0f} {kpd:12.2f} '
              f'{"yes" if 2 <= kpd <= 16 else "marginal":>15s}')

    print('\n' + '=' * 78)
    print('3. THE TUMOUR REACHES CARRYING CAPACITY, continuous arm')
    print('=' * 78)
    for rep, infl in combos:
        c = [r for r in res if r['rep'] == rep and r['influx'] == infl and r['sched'] == 'A_cont']
        days = [7, 14, 21, 28, 35, 42]
        tr = {d: np.median([[t[1] for t in r['traj'] if t[0] == d][0] for r in c]) for d in days}
        print(f'  {rep + " x " + f"{infl:.1e}":26s} ' +
              '  '.join(f'd{d}={tr[d]:.0f}' for d in days))
    print(f'\n  lattice capacity {L2} sites, initial burden 5542. The trajectories flatten by')
    print('  day 21-28 because the domain is full, not because the disease is controlled.')

    print('\n' + '=' * 78)
    print('4. THE EXHAUSTION MECHANISM ITSELF WORKS')
    print('=' * 78)
    print('  rec_low x 2.5e-5 is the only combination with real dynamic range, and is the one')
    print('  flagged in advance as physiologically sensible (stable T pool AND the required')
    print('  functional collapse). End-of-run T-cell function by schedule:\n')
    c = [r for r in res if r['rep'] == 'rec_low' and r['influx'] == 2.5e-05]
    # secondary key on the name so tied values break deterministically across processes
    for s in sorted({r['sched'] for r in c},
                    key=lambda x: (-np.median([r['f_end'] for r in c if r['sched'] == x]), x)):
        sub = [r for r in c if r['sched'] == s]
        print(f'    {s:12s} f_end={np.median([r["f_end"] for r in sub]):.3f}   '
              f'med nB42={np.median([r["nB42"] for r in sub]):6.0f}')
    print('\n  Breaks restore function exactly as calibrated - 0.315 continuous to 0.646 under')
    print('  7-on/7-off. And the tumour is still LARGER. That is an interpretable negative')
    print('  within this combination: recovery is real but does not repay the lost drug time.')
    print('  It is not, however, a valid test of the hypothesis, because the regime is one of')
    print('  uncontrolled growth against a saturating boundary.')

    print('\n' + '=' * 78)
    print('CLASSIFICATION A - MODEL NOT VALIDATED FOR THE SCHEDULING QUESTION.')
    print('The architecture hypothesis remains UNTESTED. It is not falsified.')
    print('Per PHASE2_PREREG.md section 3, any continuation requires a NEW externally')
    print('justified calibration and a NEW pre-registration. The pre-registered')
    print('non-naive-arrivals variant is explicitly NOT permitted as a rescue here.')
    print('=' * 78)


if __name__ == '__main__':
    main()
