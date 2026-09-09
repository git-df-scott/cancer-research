# Phase 5 verdict on the proliferation counterexample

# **STRUCTURAL FAILURE NARROWED**

I attacked my own counterexample and most of it did not survive. What remains is real, but it is a
different and much more specific claim than the one I made.

## What I claimed, and why it was wrong

`STRUCTURAL_COUNTEREXAMPLE.md` asserted:

> M1's chronic culture has no T-cell proliferation. Philipp's does, and it differs about four-fold
> between the two arms. The real culture sustains and expands its T cells; this model's only loses
> them, monotonically, to extinction.

**Retracted.** The 4.1-versus-1.1 figures are not the chronic culture's population trajectory. The
paper defines:

```
Fold change = CD2+ cell count day 3 / CD2+ cell count day 0
```

measured in a **3-day proliferation assay** on cells harvested from the chronic culture. It is a
functional capacity endpoint, exactly like the 72 h lysis assay — not a measurement of how many
cells the chronic culture contains. Correctly read, it says: *T cells harvested at chronic day 14
expand 4.1-fold over 3 days if rested, versus 1.1-fold if continuously stimulated.*

I misread an assay readout as a culture trajectory and built a counterexample on it.

**Philipp report no absolute viable T-cell counts for the chronic culture at any intermediate
timepoint.** So the model's 360 → 64 attrition is **unconstrained by data, not contradicted by it**.
It may be wrong; nothing measured says so. Claiming a contradiction where no measurement exists is
the same error class as the bugs this project has been finding — asserting more than the instrument
supports.

**Population size does not enter the observable anyway.** The methods specify E:T = 1:1 for the
cytotoxicity assay. A ratio cannot be set without counting effectors, so equal-number replating is
entailed by the stated design. Chronic-culture population size therefore cannot propagate into the
lysis measurement. The model's `n_plate = 250` standardisation already matches this, and turns out
to be protocol-correct for a reason I had not identified when I made it.

## What survives, stated precisely

**Philipp's two arms diverge on two functional axes, and M1 has only one latent state.**

At chronic day 14, harvested cells differ in:

| axis | continuous | TFI | ratio |
|---|---|---|---|
| specific lysis (72 h) | 34.9% | 93.4% | 2.7× |
| proliferative capacity (3 d) | 1.1 | 4.1 | 3.7× |

M1 carries a single latent variable `E`, mapped to killing through the Hill function. It has no
representation of proliferative capacity at all. So M1 cannot distinguish:

- a population restored in **both** killing and proliferation (one underlying state recovering), from
- a population restored in killing but **not** proliferation, or vice versa (two separable states).

M1 fits the lysis axis and is silent on the other. Since both axes were measured and both moved, a
fit to one of them does not establish that the single-state description is right. **The
decomposition remains unidentified — not because a channel is missing from the culture, but because
the model compresses two measured functional dimensions into one latent variable.**

This is weaker than "the model omits measured population biology" and stronger than "no problem".

## Answers to the Phase 5 questions

1. **Equal numbers plated?** Strongly implied by the stated E:T 1:1; the counting step is not
   explicitly described. Recorded as implied, not stated.
2. **Does arm-dependent proliferation still matter?** Not through population size, which the
   replating removes. Potentially through composition and division history — but neither is
   measured, so this is a possibility, not a finding.
3. **Is CD2+ fold change proliferation?** It is a net count ratio over 3 days. It cannot separate
   division from survival. Calling it "proliferation" overstates it, and I did.
4. **Independent proliferation markers?** None. No Ki-67, CFSE or CellTrace anywhere in the paper.
5. **Does the model contradict a measurement?** **No.** It omits a process whose magnitude was never
   measured in the chronic culture.
6. **Can M1 still represent the replated population's functional state?** Yes, for the killing axis.
   That is precisely what it fits, and the replating means it does not need the chronic population
   size to be right in order to do so.

## Consequence for the verdict

The counterexample does not kill M1. It also does not clear it. It converts a claimed structural
contradiction into a demonstrated **identifiability limit with a measured basis**: two functional
axes, one latent variable, both axes observed to move.

`STRUCTURAL_COUNTEREXAMPLE.md` is retained unedited as the superseded claim, with a pointer here.
