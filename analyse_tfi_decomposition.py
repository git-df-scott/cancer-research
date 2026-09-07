"""
How much of the measured treatment-free-interval benefit is REINVIGORATION, and how much is
simply DELIVERING LESS DRUG?

This is Phase H's dose-matching logic (PHASE2_PREREG.md 7.3) applied to Philipp's own experiment
rather than to the lattice model, using the parameters fitted in calib_philipp.py. It is a
diagnostic on the calibration, not a test of the architecture hypothesis, and it uses no new data.

THE DECOMPOSITION
-----------------
Philipp's TFI arm receives the engager on days 0-7 and 14-21 only: 14 days of cumulative exposure
against the continuous arm's 28. So its day-28 advantage has two possible sources, and they are
separable in the model because recovery can be switched off without touching anything else:

    CONT        28 days of exposure, recovery irrelevant             -> f_cont
    TFI         14 days of exposure, reversible component recovers   -> f_tfi
    TFI_norec   14 days of exposure, recovery DISABLED               -> f_norec

    exposure component  = f_norec - f_cont      (would happen with no reinvigoration at all)
    recovery component  = f_tfi   - f_norec     (what the rest periods actually buy back)

Because accrual acts on remaining function, exhaustion after a given total exposure is independent
of how that exposure is arranged when recovery is off. TFI_norec is therefore also exactly the
dose-matched continuous control, which is what makes the split clean.

WHY THIS MATTERS BEYOND THIS PROJECT
------------------------------------
The clinical rationale for treatment-free intervals in engager therapy is reinvigoration: rest the
T cells and they come back. If most of the measured benefit is instead explained by reduced
cumulative exposure, that rationale is weaker than it looks, and the same conclusion would apply to
any engager schedule, not only to lymphoma geometry.

This is reported whatever it shows. It was not run to support a conclusion.
"""
import json, os
import numpy as np
from calib_philipp import simulate, DATA, OUT, DAY

BIG = 1e9      # tau_r large enough that the reversible component never decays


def decompose(ch, day=28):
    k, rho, tau, c50 = ch['k_per_min'], ch['rho'], ch['tau_r_min'], ch['c50_min']
    f_cont = simulate('cont', k, rho, tau, c50)[day]
    f_tfi = simulate('tfi', k, rho, tau, c50)[day]
    f_norec = simulate('tfi', k, rho, BIG, c50)[day]
    total = f_tfi - f_cont
    expo = f_norec - f_cont
    reco = f_tfi - f_norec
    return dict(f_cont=f_cont, f_tfi=f_tfi, f_norec=f_norec,
                total=total, exposure=expo, recovery=reco,
                pct_exposure=100.0 * expo / total if total > 0 else float('nan'),
                pct_recovery=100.0 * reco / total if total > 0 else float('nan'),
                lysis_cont=100 * (1 - np.exp(-ch['a'] * f_cont)),
                lysis_tfi=100 * (1 - np.exp(-ch['a'] * f_tfi)),
                lysis_norec=100 * (1 - np.exp(-ch['a'] * f_norec)),
                matched_ratio=f_tfi / f_norec if f_norec > 0 else float('nan'))


if __name__ == '__main__':
    cal = json.load(open(f'{OUT}/calib_philipp.json'))
    print(__doc__.split('THE DECOMPOSITION')[0].strip())
    print('\nMeasured day 28: CONT 8.6% specific lysis, TFI 58.7%.')
    print('Measured matched-exposure comparison: TFI at day 21 and CONT at day 14 have both had')
    print('14 days of engager. Granzyme B ratio at those points gives a timing effect of ~1.4x.\n')
    out = {}
    for label, ch in cal['representatives'].items():
        d = decompose(ch)
        out[label] = d
        print(f'[{label}]  7-day break reverses {ch.get("recovered_7d", float("nan"))*100:5.1f}% '
              f'of accrued exhaustion   (rho={ch["rho"]:.2f}, tau_r={ch["tau_r_days"]:.2f} d)')
        print(f'    remaining function at day 28   CONT {d["f_cont"]:.4f}   '
              f'TFI {d["f_tfi"]:.4f}   TFI with recovery off {d["f_norec"]:.4f}')
        print(f'    day-28 TFI advantage           {d["total"]:.4f}   '
              f'= exposure {d["exposure"]:.4f} + recovery {d["recovery"]:.4f}')
        print(f'    share of the benefit           {d["pct_exposure"]:5.1f}% from delivering less '
              f'drug, {d["pct_recovery"]:5.1f}% from reinvigoration')
        print(f'    dose-matched timing effect     {d["matched_ratio"]:.2f}x  '
              f'(model, vs ~1.4x from measured granzyme B)')
        print(f'    as specific lysis              CONT {d["lysis_cont"]:.1f}%  '
              f'TFI {d["lysis_tfi"]:.1f}%  TFI-no-recovery {d["lysis_norec"]:.1f}%\n')
    json.dump(out, open(f'{OUT}/tfi_decomposition.json', 'w'), indent=1)
    lo = min(v['pct_recovery'] for v in out.values())
    hi = max(v['pct_recovery'] for v in out.values())
    print(f'Across the admissible parameter class, reinvigoration accounts for '
          f'{lo:.0f}%-{hi:.0f}% of the day-28 TFI benefit;')
    print(f'the remainder is reduced cumulative exposure. Both ends of that range are consistent')
    print(f'with the five measured points to within the fit quality reported in calib_philipp.py.')
    print(f'\nCAVEAT: the durable/reversible split is not identifiable (see PHASE2_PREREG.md 1.2),')
    print(f'so this range is wide by construction. It is a statement about what the data can and')
    print(f'cannot distinguish, not a precise estimate.')
    print(f'written {OUT}/tfi_decomposition.json')
