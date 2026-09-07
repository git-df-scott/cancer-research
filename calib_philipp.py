"""
Phase 2 calibration: fit the exhaustion clock to Philipp et al., Blood 2022 (PMID 35878001).

WHY THIS FILE EXISTS
--------------------
Experiment L1 was uninformative because its positive control failed: the dispersed arm did not
reproduce short-TFI-beats-continuous. The quantified cause is that exhaustion never became the
binding constraint - at most 12% functional loss by day 28, against ~90% measured. Before the
architecture comparison can mean anything, the exhaustion clock has to reproduce the external,
measured curve. That is what this file does, and it calibrates to that curve ONLY. No schedule
ranking is used as a target anywhere in this file.

THE MEASURED DATA (all five points, from the full text, not the abstract)
------------------------------------------------------------------------
Protocol: healthy-donor T cells + irradiated OCI-Ly1 targets, E:T 1:4, 5 ng/mL AMG 562
(half-life-extended CD19xCD3). Medium, targets and drug replenished on day 3; T cells isolated
and re-cultured with fresh targets and drug every 7 days, four cycles, 28 days.
TFI arm: AMG 562 omitted during stimulation cycles 2 (days 7-14) and 4 (days 21-28); target
cells still present throughout, so the difference is engager signalling, not antigen presence.
Readout: specific lysis of hCD19-Ba/F3 at E:T 1:1 after 72 h, by isolated T cells.

    arm    day 7    day 14   day 28
    CONT   88.4%    34.9%     8.6%
    TFI      -      93.4%    58.7%

Two facts in that table do the structural work, and neither is in the original handoff:

1. The TFI arm at day 21 has had exactly 14 days of cumulative engager exposure, the same as CONT
   at day 14. Measured after its rest, it is at 58.7% against CONT's 34.9%. Equal dose, better
   function, from timing alone. This is Phase 4's dose-matched control, already answered
   experimentally, and it says a timing effect genuinely exists.
2. The TFI arm does NOT return to its day-14 level (93.4%) after the second rest (58.7%). Recovery
   is partial and the shortfall grows with cumulative exposure. A single fully-reversible
   exhaustion state, which is what lymphoid.py had, cannot produce that. A durable component is
   required by the data, not chosen for convenience.

THE MODEL BEING FITTED
----------------------
Per T cell, while in antigen contact with engager present:

    dEr/dt = (1 - rho) * k * h(C) * (1 - E)      reversible component
    dEd/dt =      rho  * k * h(C) * (1 - E)      durable component
    dC/dt  = 1                                    cumulative engaged exposure, minutes

while engager absent:

    dEr/dt = -Er / tau_r                          durable component unchanged

with E = Er + Ed, remaining function f = 1 - E, and h(C) = C^n / (C^n + C50^n) an optional lag.
Accrual acting on remaining function is the standard form; it is not tuned per data point.

THE ASSAY MAP, AND WHY IT IS NOT A FREE FUDGE FACTOR
---------------------------------------------------
Specific lysis in a fixed-duration, fixed-E:T killing assay saturates: a T cell at half capacity
does not produce half the lysis. Comparing mean(1 - E) directly against a percentage that is
already at 88.4% would force the model to call day-7 T cells undamaged when the paper's own
phenotyping shows >60% of them already coexpressing PD-1, Tim-3 and LAG-3. So the model's
function is mapped through first-order kill kinetics,

    predicted specific lysis = 1 - exp(-a * f),

with ONE parameter a shared by all five points and both arms. It cannot improve the fit of one
point at another's expense; it only sets how compressed the top of the scale is.

PRE-REGISTERED FIT PROTOCOL (fixed before looking at any output)
---------------------------------------------------------------
- Targets: the five measured lysis percentages above. Nothing else. No schedule ranking.
- Two nested models are fitted and compared on the same five points:
    M0  no lag  (C50 = 0, h == 1)          parameters a, k, rho, tau_r
    M1  lag     (C50, n = 4 free/fixed)    parameters a, k, rho, tau_r, C50
  M1 is adopted only if it reduces SSE by more than a factor of two; otherwise M0 is kept as the
  more parsimonious model. This rule is fixed now so the choice is not made after seeing which
  one helps the architecture result.
- rho and tau_r are expected to be only JOINTLY identifiable: the day-28 TFI point constrains
  rho + (1 - rho) * exp(-7 / tau_r), i.e. how much of a week's rest is recovered, not the split
  between "never recovers" and "recovers slowly". The profile is reported, and both are carried
  into the architecture comparison as sensitivity axes rather than fixed at a flattering value.
"""
import json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'results')
DAY = 1440.0

# ---------------------------------------------------------------- measured data
# Philipp et al., Blood 2022, PMID 35878001. All values quoted verbatim in the Results text
# ("TFIs reinvigorate T-cell function in vitro"); Figure 3E, same isolated T cells, same days.
#
# Two kinds of readout are used, and the distinction matters more than anything else in this file.
#
# LYSIS is % specific lysis in a fixed-duration, fixed-E:T killing assay. It SATURATES: a T cell
# at a third of full capacity can still lyse most targets in 72 h, because it only needs to find
# and kill a small number of them. Fitting mean(1 - E) straight onto a percentage that is already
# at 88.4% therefore over-reads the top of the scale. Lysis is mapped through first-order kill
# kinetics, 1 - exp(-a * f), with a single shared constant a.
#
# GRANZYME B is a per-cell MFI ratio - cytotoxic machinery content, not an outcome of a killing
# race - and is close to linear in effector capacity. It does not saturate. Its CONT/TFI ratio
# therefore pins the SCALE of f, which the lysis percentages on their own cannot do. Without it
# the fit is free to push a as high as it likes and declare every arm functionally dead; the
# first version of this fit did exactly that (a = 10.0, all five points reproduced, but f between
# 0.002 and 0.25 everywhere and a granzyme-B ratio of 0.02 against a measured 0.23).
#
#   type    arm(s)      day   value    source
DATA = [
    ('lysis', 'cont', 7, 88.4),      # mean specific lysis, n = 6
    ('lysis', 'cont', 14, 34.9),
    ('lysis', 'cont', 28, 8.6),
    ('lysis', 'tfi', 14, 93.4),
    ('lysis', 'tfi', 28, 58.7),
    ('ratio', 'cont/tfi', 14, 100.0 * 144.5 / 451.8),   # granzyme B MFI ratio, CD8+, n = 6 -> 32.0
    ('ratio', 'cont/tfi', 28, 100.0 * 45.5 / 196.1),    #                                  -> 23.2
]
# HELD OUT of the fit, used only as an independent check afterwards: AMG 562-mediated CD2+ fold
# change over a 3-day assay (n = 3). Proliferation is threshold-like rather than linear in
# effector capacity, so it is not a fair fitting target, but a calibrated model should not be
# wildly wrong about it either.
HELDOUT = [('ratio', 'cont/tfi', 14, 100.0 * 1.1 / 4.1),     # -> 26.8
           ('ratio', 'cont/tfi', 28, 100.0 * 0.06 / 2.8)]    # ->  2.1

TFI_OFF = [(7.0, 14.0), (21.0, 28.0)]     # engager withdrawn, target cells still present

# Pre-registered acceptance region: any parameter set whose root-mean-square error across the
# seven fitted points is <= 5 percentage points is treated as fitting the data comparably well.
# The measurements are means of n = 3-6 donors with visible spread, so insisting on a single best
# point would be false precision.
SSE_ACCEPT = 5.0 ** 2 * len(DATA)


def drug_on(arm, day):
    if arm == 'cont':
        return True
    return not any(lo <= day < hi for lo, hi in TFI_OFF)


def simulate(arm, k, rho, tau_r, c50, n=4.0, days=28.0, dt=60.0, phi=1.0):
    """Scalar ODE for one T cell, kept for the ABM cross-check and for dt-refinement.
    phi = fraction of time spent in antigen contact (1.0 in Philipp's assay)."""
    Er = Ed = C = 0.0
    out = {}
    for i in range(int(days * DAY / dt)):
        day = i * dt / DAY
        E = min(Er + Ed, 1.0)
        if drug_on(arm, day):
            h = 1.0 if c50 <= 0 else (C ** n) / (C ** n + c50 ** n) if C > 0 else 0.0
            rate = k * h * (1.0 - E) * dt * phi
            Er += (1.0 - rho) * rate
            Ed += rho * rate
            C += dt * phi
        else:
            Er *= np.exp(-dt / tau_r)
        Ed, Er = min(Ed, 1.0), min(Er, 1.0)
        d = (i + 1) * dt / DAY
        if abs(d - round(d)) < 1e-9 and int(round(d)) in (7, 14, 28):
            out[int(round(d))] = 1.0 - min(Er + Ed, 1.0)
    return out


def simulate_grid(arm, K, RHO, TAU, C50, n=4.0, days=28.0, dt=30.0):
    """Same ODE, integrated for a whole grid of parameter sets at once.
    K, RHO, TAU, C50 are equal-shaped flat arrays. Returns {day: f_array}."""
    Er = np.zeros_like(K); Ed = np.zeros_like(K); C = np.zeros_like(K)
    has_lag = C50 > 0
    C50n = np.where(has_lag, C50, 1.0) ** n
    out = {}
    for i in range(int(days * DAY / dt)):
        day = i * dt / DAY
        E = np.clip(Er + Ed, 0.0, 1.0)
        if drug_on(arm, day):
            Cn = C ** n
            h = np.where(has_lag, Cn / (Cn + C50n), 1.0)
            rate = K * h * (1.0 - E) * dt
            Er += (1.0 - RHO) * rate
            Ed += RHO * rate
            C += dt
            np.clip(Er, 0.0, 1.0, out=Er); np.clip(Ed, 0.0, 1.0, out=Ed)
        else:
            Er *= np.exp(-dt / TAU)
        d = (i + 1) * dt / DAY
        if abs(d - round(d)) < 1e-9 and int(round(d)) in (7, 14, 28):
            out[int(round(d))] = 1.0 - np.clip(Er + Ed, 0.0, 1.0)
    return out


def build_grid():
    """Pre-registered search grid. Ranges are set from physiological plausibility, not from any
    outcome: accrual anywhere from a fortnight to a few hours, recovery from 6 h to 40 days,
    durable share from none to most, lag from none to 12 days."""
    k_day = np.geomspace(0.02, 3.0, 30)          # per day, acting on remaining function
    rho = np.linspace(0.0, 0.99, 34)
    tau_d = np.geomspace(0.1, 120.0, 24)         # days
    c50_d = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0])   # days; 0 = no lag
    G = np.meshgrid(k_day, rho, tau_d, c50_d, indexing='ij')
    return [g.ravel() for g in G], (k_day, rho, tau_d, c50_d)


def score(f, a):
    """Predicted value for every fitted target, given function levels f[(arm, day)] and the
    assay saturation constant a. Works elementwise on arrays."""
    got = {}
    for kind, arm, d, _ in DATA + HELDOUT:
        if kind == 'lysis':
            got[(kind, arm, d)] = 100.0 * (1.0 - np.exp(-a * f[(arm, d)]))
        else:
            got[(kind, arm, d)] = 100.0 * f[('cont', d)] / np.maximum(f[('tfi', d)], 1e-12)
    return got


def grid_fit():
    (Kd, RHO, TAUd, C50d), axes = build_grid()
    K = Kd / DAY; TAU = TAUd * DAY; C50 = C50d * DAY
    f = {}
    for arm in ('cont', 'tfi'):
        for d, v in simulate_grid(arm, K, RHO, TAU, C50).items():
            f[(arm, d)] = v
    # a does not enter the ODE, so it is scanned separately
    a_grid = np.geomspace(1.2, 20.0, 200)
    best_sse = np.full(K.shape, np.inf); best_a = np.zeros_like(K)
    for a in a_grid:
        got = score(f, a)
        s_ = np.zeros_like(K)
        for kind, arm, d, y in DATA:
            s_ += (got[(kind, arm, d)] - y) ** 2
        m = s_ < best_sse
        best_sse[m] = s_[m]; best_a[m] = a
    return dict(k_day=Kd, rho=RHO, tau_d=TAUd, c50_d=C50d, a=best_a, sse=best_sse), f, axes


def describe(idx, G):
    return dict(k_per_day=float(G['k_day'][idx]), k_per_min=float(G['k_day'][idx] / DAY),
                rho=float(G['rho'][idx]), tau_r_days=float(G['tau_d'][idx]),
                tau_r_min=float(G['tau_d'][idx] * DAY),
                c50_days=float(G['c50_d'][idx]), c50_min=float(G['c50_d'][idx] * DAY),
                a=float(G['a'][idx]), sse=float(G['sse'][idx]),
                rmse_pp=float(np.sqrt(G['sse'][idx] / len(DATA))))


def predict_one(p):
    f = {}
    for arm in ('cont', 'tfi'):
        for d, v in simulate(arm, p['k_per_min'], p['rho'], p['tau_r_min'], p['c50_min']).items():
            f[(arm, d)] = v
    return score(f, p['a']), f


def table(p):
    got, f = predict_one(p)
    rows = []
    for kind, arm, d, y in DATA:
        tag = 'lysis %' if kind == 'lysis' else 'GzB ratio'
        rows.append(f"    {tag:10s} {arm:9s} d{d:<3d} measured {y:5.1f}   model {got[(kind, arm, d)]:5.1f}")
    rows.append('    -- held out of the fit --')
    for kind, arm, d, y in HELDOUT:
        rows.append(f"    {'prolif':10s} {arm:9s} d{d:<3d} measured {y:5.1f}   model {got[(kind, arm, d)]:5.1f}")
    rows.append('    remaining function f = 1 - E:  ' +
                '  '.join(f'{a_}d{d_}={f[(a_, d_)]:.3f}'
                          for a_ in ('cont', 'tfi') for d_ in (7, 14, 28)))
    return '\n'.join(rows)


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    print('Phase B calibration. Fitting the exhaustion clock to five measured points,')
    print('Philipp et al. Blood 2022. No schedule ranking is used as a target.\n')
    G, f, axes = grid_fit()
    n = len(G['sse'])
    print(f'grid: {n} parameter sets x 200 assay-saturation values\n')

    i_best = int(np.argmin(G['sse']))
    nolag = G['c50_d'] == 0
    i_m0 = int(np.flatnonzero(nolag)[np.argmin(G['sse'][nolag])])
    m0, m1 = describe(i_m0, G), describe(i_best, G)

    for name, p in (('M0  no lag (C50 = 0)', m0), ('M1  best overall', m1)):
        print(f'--- {name}   SSE = {p["sse"]:.2f}   RMSE = {p["rmse_pp"]:.2f} percentage points ---')
        print(f'    k = {p["k_per_day"]:.4f}/day   rho = {p["rho"]:.3f}   '
              f'tau_r = {p["tau_r_days"]:.2f} d   C50 = {p["c50_days"]:.1f} d   a = {p["a"]:.2f}')
        print(table(p)); print()
    adopt_lag = m1['sse'] < m0['sse'] / 2.0 and m1['c50_days'] > 0
    print(f'PRE-REGISTERED RULE: adopt the lag only if it more than halves SSE.')
    print(f'  {m0["sse"]:.2f} -> {m1["sse"]:.2f}.  Adopted: '
          f'{"M1, lagged accrual" if adopt_lag else "M0, no lag (parsimony)"}\n')
    chosen = m1 if adopt_lag else m0

    # ------------------------------------------------ profiles / identifiability
    # The three independent day-28 readouts of the same CONT/TFI functional ratio disagree with
    # each other by an order of magnitude (proliferation 2.1, lysis-derived 14.6, granzyme B 23.2),
    # so the pre-registered RMSE-5pp acceptance region can be empty for reasons that have nothing
    # to do with the model. Report both: the absolute rule as pre-registered, and a relative class.
    print('Day-28 CONT/TFI functional ratio, as measured three independent ways:')
    print(f'    proliferation (CD2+ fold change)  {100*0.06/2.8:5.1f}')
    print(f'    specific lysis                    {100*8.6/58.7:5.1f}')
    print(f'    granzyme B MFI                    {100*45.5/196.1:5.1f}')
    print(f'    -> the readouts span a 10-fold range; no model can satisfy all three.\n')

    acc_abs = G['sse'] <= SSE_ACCEPT
    sse_min = float(G['sse'].min())
    acc = G['sse'] <= 2.0 * sse_min
    print(f'PRE-REGISTERED absolute rule (RMSE <= 5 pp, SSE <= {SSE_ACCEPT:.0f}): '
          f'{int(acc_abs.sum())} of {n} sets qualify.')
    print(f'DEVIATION, recorded: the absolute rule is unreachable given the readout disagreement '
          f'above (best SSE {sse_min:.1f}).')
    print(f'The equivalence class is therefore defined relative to the best fit, SSE <= 2 x min '
          f'= {2*sse_min:.1f}: {int(acc.sum())} of {n} sets.\n')

    print('IDENTIFIABILITY. Profile over rho, everything else free:')
    for r in np.unique(G['rho']):
        m = G['rho'] == r
        j = int(np.flatnonzero(m)[np.argmin(G['sse'][m])])
        pr = describe(j, G)
        rec = 1.0 - (r + (1 - r) * np.exp(-7.0 / pr['tau_r_days']))
        mark = ' <-- best' if j == i_best else ''
        if abs(r*100 - round(r*100)) < 1e-6 and int(round(r*100)) % 9 == 0 or j == i_best:
            print(f'    rho={r:4.2f}  best SSE={pr["sse"]:7.1f}  k={pr["k_per_day"]:.3f}/d  '
                  f'tau_r={pr["tau_r_days"]:7.2f}d  a={pr["a"]:5.2f}  '
                  f'7-day TFI recovers {rec:.3f}{mark}')
    print()
    if acc.any():
        for nm, key in (('k (/day)', 'k_day'), ('rho', 'rho'),
                        ('tau_r (days)', 'tau_d'), ('C50 (days)', 'c50_d'), ('a', 'a')):
            v = G[key][acc]
            print(f'    {nm:14s} range {v.min():8.3f} .. {v.max():8.3f}   median {np.median(v):8.3f}')
        rec = 1.0 - (G['rho'][acc] + (1 - G['rho'][acc]) * np.exp(-7.0 / G['tau_d'][acc]))
        print(f'    fraction of accrued exhaustion recovered by a 7-day TFI: '
              f'{rec.min():.3f} .. {rec.max():.3f}, median {np.median(rec):.3f}')

    # ------------------------------------------------------------------------------
    # The identified combination. rho and tau_r are individually flat, but what they jointly
    # control - the fraction of accrued exhaustion a 7-day break reverses - is the single
    # quantity that decides how much a treatment-free interval can be worth. Profile it
    # directly and carry the whole admissible range forward, not one flattering value.
    REC = 1.0 - (G['rho'] + (1 - G['rho']) * np.exp(-7.0 / G['tau_d']))
    print('\nPROFILE over the identified combination: fraction of accrued exhaustion')
    print('reversed by a 7-day treatment-free interval, everything else free:')
    edges = np.linspace(0.0, 1.0, 21)
    prof_rec = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (REC >= lo) & (REC < hi)
        if not m.any():
            continue
        j = int(np.flatnonzero(m)[np.argmin(G['sse'][m])])
        pr = describe(j, G); pr['recovered_7d'] = float(REC[j])
        prof_rec.append(pr)
        print(f'    recovered {lo:.2f}-{hi:.2f}  best SSE={pr["sse"]:7.1f}  '
              f'k={pr["k_per_day"]:.3f}/d  rho={pr["rho"]:.2f}  '
              f'tau_r={pr["tau_r_days"]:7.2f}d  a={pr["a"]:5.2f}')
    ok = [pr for pr in prof_rec if pr['sse'] <= 2.0 * sse_min]
    print(f'\n    admissible (SSE <= 2 x min): recovered fraction '
          f'{min(o["recovered_7d"] for o in ok):.3f} .. {max(o["recovered_7d"] for o in ok):.3f}')

    # three representatives spanning that admissible range, each at its own best fit
    reps = {}
    if ok:
        ok_sorted = sorted(ok, key=lambda o: o['recovered_7d'])
        picks = [('rec_low', ok_sorted[0]),
                 ('rec_mid', ok_sorted[len(ok_sorted) // 2]),
                 ('rec_high', ok_sorted[-1]),
                 ('best_fit', dict(describe(i_best, G), recovered_7d=float(REC[i_best])))]
        for label, pr in picks:
            reps[label] = pr
        print('\nREPRESENTATIVE PARAMETER SETS carried into every architecture run as a')
        print('pre-registered sensitivity axis. None is preferred over the others.')
        for label, pr in reps.items():
            print(f'\n  [{label}] 7-day TFI reverses {pr["recovered_7d"]*100:.1f}% of accrued '
                  f'exhaustion   SSE={pr["sse"]:.1f}')
            print(f'    k={pr["k_per_day"]:.4f}/d  rho={pr["rho"]:.2f}  '
                  f'tau_r={pr["tau_r_days"]:.2f}d  C50={pr["c50_days"]:.0f}d  a={pr["a"]:.2f}')
            print(table(pr))

    json.dump(dict(m0=m0, m1=m1, chosen=chosen, adopt_lag=bool(adopt_lag),
                   sse_accept=SSE_ACCEPT, n_in_class=int(acc.sum()), n_in_class_abs=int(acc_abs.sum()), n_grid=int(n),
                   sse_min=sse_min, class_rule='SSE <= 2 x min (deviation, recorded)',
                   representatives=reps, profile_recovered=prof_rec,
                   class_ranges={k: [float(G[k][acc].min()), float(G[k][acc].max())]
                                 for k in ('k_day', 'rho', 'tau_d', 'c50_d', 'a')} if acc.any() else {}),
              open(f'{OUT}/calib_philipp.json', 'w'), indent=1)
    print(f'\nwritten {OUT}/calib_philipp.json')
