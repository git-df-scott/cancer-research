"""
Why the Phase 3 positive control failed, and what that failure actually establishes.

Diagnosis only. Nothing is tuned and no parameter is changed. Per PHASE3_PREREG.md section 4,
two failed positive controls under two pre-registrations ends this line of attack; the model is
not modified a third time.

THE HEADLINE
------------
Unlike Phase 2, this failure is informative, because this time the endpoint had full dynamic
range: the continuous arm left 22-1296 malignant cells and the worst schedule left 7546, a
spread of two orders of magnitude, with no arm anywhere near the lattice ceiling. The experiment
could have detected a scheduling benefit of almost any size. There was none to detect. Every
non-continuous schedule was worse than continuous, monotonically in duty cycle, in all four
combinations.

THE REASON, and it is structural rather than parametric
-------------------------------------------------------
For a treatment-free interval to help, exhaustion has to be the binding constraint. In this model
it never is, in any regime where the drug also works:

  - At a low effector:target ratio (Phase 2, 1:28 to 1:7) exhaustion binds hard - terminal T-cell
    function 0.001 to 0.36 - but the tumour is never controlled. It grows past its starting burden
    and saturates the domain. The endpoint cannot measure anything.

  - At a high effector:target ratio (Phase 3, 1:4 and above) the tumour is controlled, but it is
    killed so fast that the antigen disappears before exhaustion accrues. Terminal function is
    0.53 to 1.00. There is nothing for a rest period to restore.

Those two regimes do not overlap. The sweep covered a 14-fold range of effector:target ratio, and
in 24 of 24 (N_T x parameter combination) cells, tumour control and binding exhaustion never
co-occur. The window in which an exhaustion-driven scheduling effect could exist is EMPTY.

That is why both positive controls failed, at opposite ends of the same axis, and it is one
finding rather than two accidents.

WHAT THE MODEL IS MISSING, stated precisely
-------------------------------------------
Real patients occupy the regime this model cannot reach. Philipp's Figure 1B measured exactly it:
patients on continuous blinatumomab with persistent disease AND peripheral T-cell function down
to 0.238 of baseline. Partial control with persistent antigen and a progressively exhausting
effector pool is the clinically relevant state, and it is precisely the state that is unreachable
here.

It is unreachable for one identifiable reason. Sustaining an engaged effector pool against
persistent antigen requires ongoing recruitment; but in this model recruited T cells arrive with
zero exhaustion, so any influx rate high enough to sustain the pool also resets the population
mean and prevents functional collapse. That is what Phase C measured and excluded, and it is what
PHASE2_PREREG.md section 2 recorded as a structural limitation BEFORE any of these runs:

    "within this model's structure, recruitment and population-level functional collapse cannot
     both be represented at a realistic influx rate, because arrivals enter with zero exhaustion"

Every thread in this project converges on that single sentence.

WHAT WOULD BE REQUIRED, and it is not attempted here
----------------------------------------------------
A model able to test the architecture hypothesis needs recruited T cells that are not naive - an
exhaustion state that is partly systemic rather than purely per-cell, so that a sustained pool can
still lose function. That is a different model, requiring its own external calibration and its own
pre-registration. PHASE3_PREREG.md section 4 forbids reaching for it now, and it is not run.
"""
import json
import numpy as np

CAP = 120 * 120


def main():
    print(__doc__)
    sw = json.load(open('results/expNT_sweep.json'))
    p3 = json.load(open('results/expE3_poscontrol.json'))

    print('=' * 84)
    print('1. THE ENDPOINT HAD DYNAMIC RANGE THIS TIME. The failure is real, not an artefact.')
    print('=' * 84)
    print(f'  {"combination":26s} {"cont nB42":>10s} {"worst arm":>10s} {"ratio":>7s} '
          f'{"peak %cap":>10s}  ceiling?')
    for rep, infl in sorted({(r['rep'], r['influx']) for r in p3}):
        sub = [r for r in p3 if r['rep'] == rep and r['influx'] == infl]
        cont = np.median([r['nB42'] for r in sub if r['sched'] == 'A_cont'])
        worst = max(np.median([r['nB42'] for r in sub if r['sched'] == s])
                    for s in {r['sched'] for r in sub})
        peak = max(max(t[1] for t in r['traj']) for r in sub) / CAP
        print(f'  {rep + " x " + f"{infl:.0e}":26s} {cont:10.0f} {worst:10.0f} '
              f'{worst/max(cont,1):6.0f}x {100*peak:9.1f}%  '
              f'{"HIT" if peak >= 0.70 else "no"}')
    print('\n  No combination approaches the 70% ceiling. The pre-registered ceiling contingency')
    print('  does not trigger; all four combinations count, and all four failed.')

    print('\n' + '=' * 84)
    print('2. CONTROL AND EXHAUSTION NEVER CO-OCCUR. 24 of 24 cells, 14-fold range of E:T.')
    print('=' * 84)
    print(f'  {"N_T":>5s} {"E:T":>6s} {"nB42/n0":>8s} {"f_end":>7s}   controlled   exhausted   both')
    both_any = False
    for n_t in sorted({r['n_t'] for r in sw}):
        rows = []
        for rep, infl in sorted({(r['rep'], r['influx']) for r in sw}):
            s = [r for r in sw if r['n_t'] == n_t and r['rep'] == rep and r['influx'] == infl]
            if not s:
                continue
            n0 = np.median([r['n0'] for r in s])
            rows.append((np.median([r['nB42'] for r in s]) / n0,
                         np.median([r['f_end'] for r in s]), n0))
        rat = np.median([r[0] for r in rows]); f = np.median([r[1] for r in rows])
        ctl, exh = rat < 1.0, f < 0.5
        both_any |= (ctl and exh)
        print(f'  {n_t:5d} {"1:"+f"{rows[0][2]/n_t:.0f}":>6s} {rat:8.2f} {f:7.3f}   '
              f'{"yes" if ctl else "no ":>10s}  {"yes" if exh else "no ":>10s}   '
              f'{"** YES **" if (ctl and exh) else "no"}')
    print(f'\n  Any regime with both tumour control and binding exhaustion: {both_any}')
    print('  The window in which an exhaustion-driven scheduling effect could exist is EMPTY.')

    print('\n' + '=' * 84)
    print('3. EVERY SCHEDULE IS WORSE THAN CONTINUOUS, MONOTONICALLY IN DUTY CYCLE')
    print('=' * 84)
    sub = [r for r in p3 if r['rep'] == 'rec_low' and r['influx'] == 2.5e-05]
    print('  rec_low x 2.5e-5, the combination flagged in advance as physiologically sensible:\n')
    print(f'  {"schedule":12s} {"duty":>5s} {"med nB42":>9s} {"f_end":>7s}')
    # secondary key on the name: three schedules share duty 0.86, and set iteration order
    # varies between processes under hash randomisation, so ties must break deterministically
    for s in sorted({r['sched'] for r in sub},
                    key=lambda x: (-np.median([r['duty'] for r in sub if r['sched'] == x]), x)):
        ss = [r for r in sub if r['sched'] == s]
        print(f'  {s:12s} {ss[0]["duty"]:5.2f} {np.median([r["nB42"] for r in ss]):9.0f} '
              f'{np.median([r["f_end"] for r in ss]):7.3f}')
    print('\n  Terminal T-cell function is 0.90 under continuous dosing. Exhaustion is simply not')
    print('  the binding constraint at this effector:target ratio, so a rest period has nothing')
    print('  to restore and only costs drug time.')

    print('\n' + '=' * 84)
    print('CLASSIFICATION A - MODEL NOT VALIDATED FOR THE SCHEDULING QUESTION.')
    print('The architecture hypothesis remains UNTESTED. It is not falsified.')
    print('Two positive controls have now failed under two pre-registrations. Per')
    print('PHASE3_PREREG.md section 4 this line of attack ends here and is written up as a')
    print('negative. The model is NOT modified a third time.')
    print('=' * 84)


if __name__ == '__main__':
    main()
