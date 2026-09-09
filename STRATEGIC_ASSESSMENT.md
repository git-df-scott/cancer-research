# Strategic assessment: what four instrument bugs in a row actually mean

Written while the M0/M1 refit runs, before its results are known, so it cannot be a rationalisation
of them.

## 1. The bug pattern is structural, not bad luck

Four instrument bugs so far in this calibration:

| # | Bug | How it was found |
|---|---|---|
| 1 | Destructive probe restoring lattice-indexed state onto migrated cells | day-14 lysis of 50.9% was inconsistent with measured mean exhaustion 0.948 |
| 2 | Filename rounding collision silently skipping a cell | cell count was 68 where 69 were expected |
| 3 | Wrong observable entirely (4 h kill ratio vs 72 h control-normalised lysis) | reading the source methods rather than our summary of them |
| 4 | Plating-density confound masking a working recovery mechanism | a mandated parameter-liveness test showed `recover_tau` inert |

None raised an error. All four produced plausible, internally consistent numbers. Every one was
caught by **an inconsistency between two independently computed quantities** — never by something
failing.

That is a property of the enterprise, not a run of bad luck. A simulator of this kind has no
independent check on its outputs: it is flexible enough that almost any number it emits looks
biologically plausible. Plausibility carries no information about correctness here.

**Consequence for how effort should be spent.** Additional grid cells add essentially zero
information about correctness. Cross-checks between two independently computed quantities add a
great deal. The rate-limiting resource in this project is not compute, it is *independent
constraints*. Six of the seven invariants in `philipp_assay.py` now exist because of a bug that had
already shipped; that ratio should be inverted.

## 2. The inferential chain is longer than the evidence supporting it

What the calibration is actually being asked to support:

```
AMG 562, CD19, in vitro, B-lymphoma line, 2D lattice, 4 numbers from 1 paper
  -> exhaustion accrual and recovery parameters
    -> tarlatamab, DLL3, in vivo, SCLC nests, PK-driven graded occupancy
      -> a schedule ranking
```

Every arrow is a substantial extrapolation, and the first one is being fitted from three numbers
after the fourth was found contaminated. Even a perfect fit would leave the downstream claim
resting on cross-molecule, cross-target, cross-compartment and in-vitro-to-in-vivo transfers, none
of which the calibration constrains at all.

This is worth stating plainly because the calibration passing would not make the lung prediction
sound. It would remove one of several reasons it is currently unsound.

## 3. The observation model is compressive, which bounds what any fit can achieve

Measured on the corrected readout: at the fitted kill rate the observable is near-saturated across
a wide band of exhaustion, and collapses only near the top of the range. Where `dY/dE ~ 0`, the
assay contains essentially no information about latent `E`.

So the four Philipp numbers cannot identify per-cell exhaustion across its range. They constrain it
only in the narrow band where the readout is responsive. Any parameter reported from this fit is
identified in that band and extrapolated outside it — and R and T would operate substantially
outside it, because a spatially structured tumour spends most of its time at contact fractions and
occupancies the assay never probes.

That is an argument for classification C independent of how the refit lands.

## 4. The most defensible result so far required no simulation

The PK occupancy analysis is the strongest thing this project has produced:

- terminal half-life 5.8–11.2 d against a 14-day interval,
- steady-state trough occupancy drop 17–22% at EC50 ~1 nM, 37–45% at 3 nM,
- so the hypothesis holds or dies on a quantity not established in vivo.

It depends on published population-PK parameters and arithmetic, is checkable by anyone, and has
survived every correction round untouched. The ABM has not yet added anything to it.

**What the ABM would add, if calibrated,** is whether T cells are in *contact* while drug is
present — occupancy is about drug, contact is what exhaustion integrates. But that requires the
SCLC nest geometry, which is unobtainable so far, and the exhaustion calibration, which is fragile.

The honest strategic read: the ABM may not be the right instrument for the lung question. The
PK-level statement plus a well-posed identifiability result may be the real deliverable, and that
would be a legitimate outcome rather than a failure.

## 5. What I expect, and what would change my mind

Expectation before seeing the refit: **C — underidentified**. Both M0 and M1 will likely reproduce
three continuous-arm points, because the observation model is compressive enough to be forgiving,
while implying different behaviour where the assay is uninformative. No provenance-clean validation
exists — every Philipp value is contaminated through `recover_tau` and `exhaust_tonic` — so A is
unreachable regardless of fit quality.

What would move me to **D**: a systematic residual that survives across both mechanisms and the
whole bracketed parameter region, in a regime where the readout is demonstrably responsive. The
previous move toward D was withdrawn precisely because its residual came from a compressed,
artefact-driven regime.

What would move me to **B** rather than C: M0 and M1 agreeing closely on diagnostic quantities
across the admissible region, which would mean the underdetermination is harmless.

## 6. The discriminating experiment, if C

M0 and M1 differ in the *shape* of the mapping from latent state to function — linear versus
threshold. They agree wherever the readout saturates and disagree where it is responsive, so the
discriminating experiment must sit in the responsive band and vary the quantity the two shapes
treat differently.

The natural axis is **drug-free interval duration**, titrated rather than binary. Philipp used one
7-day rest. A linear recovery mechanism predicts function returning roughly in proportion to rest
duration; a threshold mechanism predicts little return until the latent variable crosses back below
threshold, then a sharp jump.

Sketch, to be specified properly if C is the classification: fixed 14-day chronic stimulation to
reach an intermediate exhaustion state, then rest arms of 0, 1, 2, 3, 5, 7, 10 days, each read out
by the standard 72 h assay at fixed plating density. The prediction that separates the models is
the *shape* of recovery versus rest duration, not any single point — which is exactly what one
7-day observation cannot distinguish.
