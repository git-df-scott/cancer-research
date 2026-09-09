# The single highest-value next experiment

**Question:** does one latent functional state drive both engager-induced loss of killing *and*
loss of expansion capacity, or are these two separable states that merely correlate?

This is the identifiability failure that current evidence cannot resolve, quantified below. The
experiment is designed to break it, not to support either answer.

## Why existing data cannot answer it

Philipp measured both axes, at two arms and two timepoints:

| condition | 72 h specific lysis | 3-day CD2+ expansion |
|---|---|---|
| d14 continuous | 34.9% | 1.1× |
| d14 TFI | 93.4% | 4.1× |
| d28 continuous | 8.6% | 0.06× |
| d28 TFI | 58.7% | 2.8× |

Sorted by lysis, expansion increases monotonically (0.06, 1.1, 2.8, 4.1), so **a single monotone
mapping is consistent with every point**. The one-state hypothesis is not falsified.

It is also barely tested:

- A perfectly monotone ordering arises by chance in 1 of 4! = **4.2%** of orderings.
- The four points are not independent: day 28 is below day 14 on *both* axes in *both* arms simply
  because time passes, so a shared time trend manufactures correlation.
- The only time-controlled comparisons are the two within-day ones, both concordant. Chance
  concordance is **0.5² = 0.25**.

**Two concordant comparisons cannot distinguish one latent state from two correlated states.** That
is the identifiability failure, stated with its power.

## The design

**Principle.** One state predicts the two axes recover with the *same* time constant, so their
normalised recovery curves superimpose. Two states permit different time constants — plausibly
fast for cytolytic function, which may need only granzyme reloading, and slow for proliferative
capacity, which is transcriptionally and epigenetically gated. Titrating rest duration separates
them; a single 7-day rest cannot, which is exactly why Philipp's design leaves it open.

| | |
|---|---|
| **System** | Healthy-donor T cells + irradiated OCI-Ly1, as in Philipp, so the chronic arm is directly comparable |
| **Engager** | AMG 562 (or blinatumomab, stated either way) at the established chronic concentration |
| **Chronic phase** | Continuous stimulation, E:T 1:4, to **day 14** — chosen because both axes are mid-range there (34.9%, 1.1×) and so have dynamic range in both directions |
| **Rest arms** | Drug removed for **0, 1, 2, 3, 5, 7, 10, 14 days**, targets still present, matching Philipp's TFI condition |
| **Replicates** | n = 6 donors, matching the power of the calibration data |

**At every rest duration, measure all four:**

1. **72 h specific lysis**, E:T 1:1 vs fresh targets, **equal-number replating**, normalised to a
   no-engager control construct — identical to Philipp's readout so the arms are comparable.
2. **3-day CD2+ expansion**, identical to Philipp's readout.
3. **Absolute viable T-cell count** at harvest. Never reported by Philipp, and its absence is why
   the model's chronic attrition is currently unconstrained rather than tested.
4. **Division history (CellTrace/CFSE) plus a viability stain.** CD2+ fold change is a net count
   ratio and cannot separate division from survival — this resolves that ambiguity directly.

**Controls:** unstimulated T cells carried in parallel; a no-engager chronic arm; a continuously
stimulated arm sampled at matched calendar times, so recovery is measured against elapsed time
rather than against day 14 alone.

## Predictions and falsification

Normalise each axis to its own maximum across rest durations and fit a half-recovery time,
t½(lysis) and t½(expansion).

| | prediction | falsified by |
|---|---|---|
| **H1 — one latent state** | the two normalised curves superimpose; ratio t½(expansion)/t½(lysis) ≈ 1 | ratio outside **[0.5, 2.0]** |
| **H2 — two separable states** | curves separate, with expansion recovering more slowly | ratio inside [0.5, 2.0] and curves superimposing within replicate error |

**The decisive quantity is the ratio of half-recovery times, not either curve alone.** It is
dimensionless, insensitive to the absolute scale of either assay, and unaffected by the plating
normalisation that removes population abundance.

A ratio outside [0.5, 2.0] **kills the single-latent-state description**, and with it M1's claim to
represent transferable exhaustion biology rather than one fitted axis.

## What it costs and why it is worth it

Eight rest arms × four readouts × six donors, in an assay system already established and published.
No new reagents, no animal work, no new molecule. It reuses Philipp's exact chronic culture and
lysis assay, adding rest durations and two measurements that were already being made in adjacent
figures.

It is the smallest experiment that can distinguish the surviving explanations, and its result is
decisive in either direction — which is the property that makes it worth doing rather than
another simulation.
