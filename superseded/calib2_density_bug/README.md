# Superseded: M0 and M1 fits invalidated by a plating-density confound in the readout

**Every fit in this directory is wrong.** Preserved, not deleted.

## The bug

`readout_specific_lysis` plated however many T cells survived the chronic culture. Philipp plate a
**defined number** of effectors at E:T 1:1. Because the chronic culture loses T cells with no
influx or proliferation, each timepoint was read out at a different density:

| timepoint | surviving T cells | mean exhaustion |
|---|---|---|
| d7 | 229 | 0.131 |
| d14 | 137 | 0.318 |
| d28 | 51 | 0.779 |
| d14 TFI | 134 | 0.048 |

Density alone drives the readout, at fixed exhaustion E = 0.10:

| T cells plated | lysis % |
|---|---|
| 600 | 100.0 |
| 300 | 99.0 |
| 120 | 88.3 |
| 60 | 51.7 |

So the observable confounded per-cell function with population attrition.

## What it hid

The day-14 TFI arm reaches mean exhaustion **0.048** against the continuous arm's **0.318**.
Reinvigoration works in the model. It read out at 65% purely because only 134 cells were plated,
and the resulting −25.9 residual was reported as a systematic failure of the recovery mechanism.

It also made `recover_tau` appear computationally **dead**: a 14-fold change moved d14_tfi by 1.5
points. After standardising density the same sweep spans 84.8–89.6, and the target of 93.4 is
inside the predeclared ±15 tolerance.

## Conclusions invalidated

- M0's best fit (SSE 204.2) and M1's (311.3), and the comparison between them.
- "No cell both fits and passes" — the passing/failing structure was density-driven.
- "Both models under-predict recovery by ~25 points, systematically."
- The inference that the recovery residual is a model-adequacy signature.
- Any move toward classification D on that residual.

## How it was found

By the parameter-liveness test: `recover_tau` had been held fixed in every cell, and a rejection
claim cannot rest on a residual produced by a mechanism never shown to be computationally live.
The liveness test showed it was inert, and asking *why* exposed the confound.

The bug produced plausible, internally consistent numbers and passed all five non-destructiveness
invariants — those check that measurement does not mutate state, not that it measures the right
quantity.

## Fix

The readout now samples `n_plate` cells (default 250) from the harvested exhaustion distribution,
with replacement when the culture cannot supply that many, holding density fixed across timepoints
and arms. Sampling preserves the exhaustion distribution, which is what the assay interrogates.

Recorded assumption: a real culture has finite yield, so representing a depleted culture at full
plating density is optimistic. That is a deviation from the protocol and is listed as such.
