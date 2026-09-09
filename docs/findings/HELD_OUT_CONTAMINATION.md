# Phase F integrity failure: the held-out observation was already inside the model

Recorded before the grid extension finishes, so it cannot be presented as a reaction to the
outcome.

## What I designated

`fit_calibration.py` declares `d14_tfi` (Philipp Fig 3E, 93.4% specific lysis at day 14 after a
7-day drug-free interval) as **held out**. The fit uses only the three continuous-arm points. The
stated justification was that "nothing in the fit informs `recover_tau`, so it is a genuine
out-of-sample prediction."

## Why that justification is false

`lymphoid.py:91`:

```python
recover_tau=10080.0,     # 7 d; TFI reinvigoration (Philipp 2022, Weber Science 2021)
```

The recovery time constant was chosen with explicit reference to Philipp 2022 — the same paper,
and specifically its treatment-free-interval result, that `d14_tfi` is drawn from. The observation
had already influenced model development before I designated it as held out.

So the fit does not inform `recover_tau`, but the *default value* of `recover_tau` was informed by
the held-out observation. That is the leak. A parameter set by looking at a number is not
independent of that number merely because a later fitting procedure left it alone.

This is the failure mode the Phase F brief names directly: an outcome that has already influenced
model development cannot serve as a validation target.

## It fails anyway, which makes this worse rather than better

Across all 16 base-grid cells, **0/16 pass** the predeclared ±15 percentage-point tolerance. The
best-fitting cell predicts 67.5 against 93.4, an error of −25.9. Every cell under-predicts recovery,
and the sign is systematic.

So the model fails a test that was tilted in its favour. Had it passed, the pass would have been
uninterpretable; failing, the failure is informative — it says the recovery behaviour is wrong in a
direction that a favourably-set parameter could not rescue.

## Consequences, stated plainly

1. **No claim of held-out validation can be made from `d14_tfi`.** Not now, and not if a later grid
   happens to pass it.
2. **Classification A is unreachable with the current target set.** There is no uncontaminated
   held-out observation available among the four Philipp values, because all four come from the
   same paper that informed the model's recovery parameter.
3. A genuine held-out test requires an observation from a source that demonstrably did not inform
   any model default. Candidates worth pursuing, none yet obtained:
   - a different chronic-stimulation dataset, ideally a different engager or target,
   - a dose or schedule arm from Philipp not used in any parameter choice, if one exists in the
     supplement,
   - a prospective prediction registered before retrieving the corresponding data.

## What should have caught this

Nothing in the pipeline checked parameter provenance against the held-out target. The check is
cheap and is now owed: for every parameter default, record the source that set it; for every
held-out observation, assert no parameter default cites that source.

That test does not exist yet and is listed as required work rather than claimed as done.

---

## Addendum: the cells that "pass" do so by saturation

The extension grid brackets the optimum (the best-fit cell is now interior in both dimensions), and
it produces six cells that pass the held-out tolerance. None of them is evidence of anything.

| p_kill | tonic | d7 | d14 | d28 | fit SSE | d14_tfi |
|---|---|---|---|---|---|---|
| 1.0e-03 | 2.5e-05 | 79.0 | 45.6 | 9.7 | **204.2** | 67.5 FAIL |
| 2.0e-03 | 2.5e-05 | 96.8 | 87.9 | 28.5 | 3271.6 | 90.3 PASS |
| 4.0e-03 | 1.2e-05 | 99.8 | 98.1 | 66.0 | 7421.0 | 99.0 PASS |
| target | | 88.4 | 34.9 | 8.6 | | 93.4 |

Every passing cell sits at p_kill 2e-3 or 4e-3, where the readout is saturated: day 7, day 14 and
the TFI arm all read 87–99% because the effectors clear the targets regardless of exhaustion. They
pass the recovery test by being unable to represent loss of function at all. Their fit SSE is 16–36×
worse than the best-fitting cell.

**No cell both fits and passes.** The two criteria select disjoint regions of parameter space, and
they do so for a comprehensible reason: fitting the decay requires enough exhaustion sensitivity to
drop day 28 to 8.6, while passing the recovery test requires enough kill capacity to reach 93.4
after a 7-day break, and this observation model cannot do both at one kill rate.

That is a model-adequacy signature, not a search problem. The optimum is bracketed; there is no
unexplored region left to appeal to within this grid.
