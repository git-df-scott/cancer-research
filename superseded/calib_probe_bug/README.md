# Superseded: calibration scans invalidated by a probe bug

**Every result in this directory is wrong. Nothing here should be cited.** It is kept because the
project's practice is to preserve failed attempts rather than delete them, and because the failure
mode is instructive.

## What was wrong

`calibrate_exhaustion.Assay.lysis_probe` measured the exhausted-state kill rate, then zeroed `E` to
measure a naive reference rate, then restored the saved array:

```python
saved_E = self.m.E.copy()
...run 240 min...            # exhausted rate
self.m.E = np.zeros_like(self.m.E)
...run 240 min...            # naive rate
self.m.E = saved_E           # <- the bug
```

`E` is indexed by lattice position, and T cells migrate at roughly one site per minute. Over a
480-minute probe their positions decorrelate completely. Restoring the saved array therefore gave
each surviving T cell whatever exhaustion value happened to sit at its **new** coordinates, which
was almost always zero.

So every probe silently reset the population's exhaustion. Measured directly:

```
before d7 probe:  E.mean = 0.473
after  d7 probe:  E.mean = 0.024
```

## What that did to the results

The day-14 reading never measured 14 days of exhaustion. It measured 7 days, twice. The day-28
reading measured 14 days. Every simulated decay curve was flattened, and the flattening grew with
the number of probes already taken.

The same parameters, before and after the fix:

| | d7 | d14 | d28 |
|---|---|---|---|
| buggy probe | 55.7 | 50.9 | 8.6 |
| correct probe | 61.4 | 4.6 | 0.0 |
| Philipp target | 88.4 | 34.9 | 8.6 |

## What it invalidated

- The 32-cell linear scan (`linear_cells/`).
- The 69-cell threshold scan (`threshold_cells/`), across all three grids.
- The empirical claim that "neither mechanism fits", and the specific residual pattern that
  claim rested on — a curve "too flat between days 7 and 14" — which was the bug's signature, not
  the model's behaviour.
- The follow-on hypothesis that lattice contact heterogeneity was flattening the transition. That
  was already disproved by direct measurement before the bug was found: exhaustion across the
  T-cell population is highly uniform, CV ≈ 0.011.

## What survived

The arithmetic argument that a linear exhaustion mechanism cannot reproduce Philipp's curve does
**not** depend on any simulation. Inverting `lysis = naive x (1 - E)` requires E(14)/E(7) between
5.6 and unbounded; linear accrual forces exactly 2.0. That result stands, and the corrected probe
confirms the model does behave as the arithmetic says: at `exhaust_tonic = 5e-5`, measured
E(7) = 0.473 and E(14) = 0.948, a ratio of 2.00.

## How it was found

Not by a test. By instrumenting the assay to measure the exhaustion distribution while chasing a
different question, and noticing that a reported day-14 lysis of 50.9% was inconsistent with a
measured mean exhaustion of 0.948, which implies about 5%. The two numbers could not both be true.

The lesson worth keeping: the probe mutated the state it was measuring, and produced plausible
numbers while doing so. It never raised an error. The corrected version runs on deep copies, so it
cannot perturb what it measures, and a regression check asserts `E` is identical before and after.
