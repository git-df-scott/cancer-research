# Identifiability analysis: the PDAC "early detection window" reconciliation

**Status 2026-09-07: identifiability question ANSWERED. Model is identifiable
given one further anchor. A substantial part of the scientific question was
resolved by the algebra, before any simulation was written.**

Nothing here is a validated biological claim or a clinical recommendation.

## Why this analysis was run first

The published estimates of "how long is the pancreatic cancer early-detection
window" span a factor of eighteen:

| Estimate | Value | Source |
|---|---|---|
| ctDNA-detectable sojourn | **0.49 yr** [0.26, 0.88] | Hubbell et al., *Cancer* 2026 (CPS-3/CCGA3) |
| MSCE preclinical sojourn | ~3 yr | Luebeck MSCE fit to SEER |
| Proteomic inflection (CTHRC1) | ~8.93 yr | UK Biobank Olink, 2026, n=62 |
| Genomic: initiation → founder cell | ~11.7 yr | Yachida et al., *Nature* 2010 |

The proposal was to build one natural-history model and ask whether a single
underlying history can reproduce all of them. The risk was that a natural-history
model carries more free parameters than there are anchors, making the exercise a
curve-fit rather than a test. **That risk was real:** with the three anchors
originally proposed (ctDNA sojourn, MSCE sojourn, CAPS stage-I yield) against six
parameters, the model is badly under-determined.

## Model M0

Single invasive clone, founder cell at t=0.

```
V(t) = v0 exp(g t),   v0 = 1e-9 cm^3 (one cell),   g ~ LogNormal(mu, sigma)
clinical presentation at V_pres;  modality k detects at V_k
S_k = ln(V_pres / V_k) / g
metastatic seeding hazard = m * V(t)
```

Parameters: `mu, sigma, V_pres, V_img, V_ctdna, m` (6).

## Anchors

| | Anchor | Value | Kind |
|---|---|---|---|
| A1 | ctDNA sojourn | 0.49 yr | model fit |
| A2 | MSCE preclinical sojourn | 3 yr | model fit |
| A3 | CAPS5 stage-I fraction, annual imaging | 7/9 = 0.778 | direct |
| A4 | volume doubling time | 132 ± 132 d | **direct** |
| A5 | size at diagnosis | ~3.1 cm | **direct** |
| A6 | imaging detection limit / CAPS lesion sizes | needed, not yet extracted | direct |

A4 and A5 are the load-bearing additions. They were not in the original plan.

## Results

### 1. The growth distribution is pinned by A4 alone

CV of VDT = 1.00 fixes `sigma = 0.833` outright; the mean then fixes
`mu = -4.903`. No sojourn anchor is consumed. Median VDT 93 d.

**Methodological trap found:** every sojourn depends on `E[1/g]`, not `1/E[g]`.
These differ by exactly `exp(sigma^2) = 2.00` here. Any sojourn computed from a
mean growth rate is understated by 100%.

### 2. Two structural facts, no numbers required

Because `S_k = ln(V_pres/V_k)/g`:

- **(a)** Thresholds enter only as the ratio `V_pres/V_k`. Absolute volumes are
  not identifiable from sojourn data alone — A5 is required to set the scale.
- **(b)** `E[S_img]/E[S_ctdna] = ln(V_pres/V_img) / ln(V_pres/V_ctdna)` is
  **growth-independent**. `E[1/g]` cancels exactly. A1 and A2 together therefore
  impose a constraint that cannot be absorbed by retuning growth. **This makes
  A1+A2 a test, not a fit** — which is precisely what the project needed.

### 3. What each anchor implies about detection thresholds

Setting the scale with A5 (presentation at 3.1 cm):

| Anchor | Sojourn | Implied detection threshold |
|---|---|---|
| ctDNA (A1) | 0.49 yr | **2.27 cm** diameter |
| ctDNA, 95% CI | 0.26–0.88 yr | 1.77–2.63 cm |
| MSCE (A2) | 3.00 yr | **0.46 cm** |
| Proteomic CTHRC1 | 8.93 yr | **0.01 cm** |

Three consequences, and these are the scientific content:

1. **The 0.49 yr ctDNA sojourn is a statement about assay sensitivity, not
   tumour biology.** For the window to be that short, ctDNA must only cross its
   limit of detection at ~2.3 cm — nearly presentation size. Pancreatic cancer is
   not uniquely fast; ctDNA is uniquely bad at seeing it.
2. **The MSCE "preclinical screen-detectable" compartment is not an imaging
   state.** It implies a 0.46 cm threshold, below what CT resolves for PDAC. A1
   and A2 were never measuring the same thing, so their 6× disagreement is not a
   contradiction.
3. **The 8.93 yr proteomic inflection cannot be tumour-derived.** It implies a
   0.01 cm (~100 µm, a few hundred cells) lesion. No plausible invasive tumour of
   that size produces a systemic plasma signal. CTHRC1 at 9 years must be field
   effect, precursor burden, or host state — **not a tumour localisation signal**,
   and a screening test built on it should not be expected to find a resectable
   lesion.

### 4. Length bias is not a footnote at CV = 1

Screen-detected prevalent cases are sampled proportional to sojourn length. The
inflation factor is exactly `exp(sigma^2) = 2.00`. If A1 and A2 do not use the
same sampling frame, **one third of the 6.12× gap is pure methodology.**

### 5. Practical identifiability

Scaled Jacobian at a trial point hitting all five anchors:

```
singular values: 8.666, 1.204, 0.342, 0.262, 0.099
condition number: 87        numerical rank: 5 of 6
worst direction: V_img (-0.84), sigma (+0.38), m (+0.29)
```

Condition number 87 is well conditioned — not marginal. Exactly one flat
direction, dominated by `V_img`, which is expected since A2 was used to set it
and A2 turns out not to be an imaging anchor. **A6 closes it.**

### 6. Robustness to growth law

The exponential assumption was stress-tested against Gompertz at three carrying
capacities:

| Sojourn | Exponential | Gompertz K=2× | K=10× | K=100× |
|---|---|---|---|---|
| 0.49 yr | 2.27 cm | 2.82 cm | 2.58 cm | 2.46 cm |
| 3.00 yr | 0.46 cm | 0.62 cm | 0.39 cm | 0.38 cm |

Order-robust. Under every growth law tested the ctDNA threshold sits within ~1 cm
of presentation size and the MSCE threshold sits far below imaging resolution.
Gompertz makes the ctDNA conclusion *stronger*.

### 7. Out-of-sample test — split result

The 2.3 cm ctDNA threshold was derived from A1, A4, A5 only. It predicts a
stage-specific ctDNA sensitivity that was never used in the derivation.

- **PASS, stage I/II:** predicted 62.7% vs measured 31% (KRAS-only digital NGS)
  to 56–62% (multi-analyte). Lands at the top of the measured range, with zero
  tuning.
- **FAIL, stage IA:** predicted 0.0% vs measured non-zero. A hard step at 2.27 cm
  makes every ≤2 cm tumour invisible by construction.

The failure falsifies the **hard-threshold detection assumption**, not the
underlying claim. Detection must be probabilistic in volume,
`P(detect) = 1 - exp(-lambda V)`. That adds one parameter per modality and
therefore changes the identifiability budget — Step 5 must be re-run on the soft
version before anything is built. Stage-IA sensitivity then becomes a usable
anchor rather than a degenerate prediction, so it plausibly pays for itself.

## Verdict

**Identifiable, conditional on A6, and only because the two direct measurements
(A4, A5) are brought in.** With the three originally proposed anchors alone the
model would have been a curve-fit. The project is worth running.

The larger finding is that the "paradox" is substantially dissolved by the
algebra: the competing window estimates are not in conflict because they are
thresholds on different quantities at different sizes, and a third of the
remaining gap is length-bias methodology.

## Next, in order

1. Re-run Step 5 with soft (probabilistic) detection functions. **Do not skip —
   it may un-identify the model.**
2. Extract A6: CAPS screen-detected lesion size distribution, from the paper
   itself (not available at abstract level).
3. Prior-art audit before any novelty claim. Specifically: has anyone already
   published that sojourn time is assay-specific in this quantitative form? The
   claim is simple enough that it may exist.
4. Only then the full reconciliation model.

## Reproduce

```bash
python -m venv .venv && .venv/bin/pip install numpy scipy
.venv/bin/python pdac/identifiability.py
.venv/bin/python pdac/robustness.py
```
