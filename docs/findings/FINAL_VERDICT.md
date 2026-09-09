# Final verdict: **M1 SURVIVES CALIBRATION — NOT VALIDATED, MECHANISM UNIDENTIFIED**

The clean reselection is complete and M1 was not killed. It is also not validated, and the thing
it would need to be validated *for* — a causal account of treatment-free-interval recovery — is
not identified by any available evidence.

## Phase 1 result: the contamination did not drive the choice

225 cells, grid defined from the reachable domain and structural limits **before** running, with
theta and hill ladders deliberately not centred on the historical values.

| theta | hill | p_kill | tonic | d7 | d14 | d28 | SSE |
|---|---|---|---|---|---|---|---|
| **0.50** | **4.0** | **1e-3** | **5e-5** | **87.0** | **33.4** | **10.8** | **9.0** |
| 0.30 | 4.0 | 2e-3 | 5e-5 | 97.6 | 43.6 | 3.6 | 185.3 |
| 0.15 | 4.0 | 2e-3 | 2.5e-5 | 97.6 | 43.6 | 1.2 | 215.1 |
| *target* | | | | *88.4* | *34.9* | *8.6* | |

- **Exactly one cell of 225** meets both frozen criteria (SSE < 50 **and** every |residual| ≤ 10).
- Residuals −1.4 / −1.5 / +2.2 against seed SDs 0.2 / 3.8 / 1.6 — within the simulation's own noise.
- **Interior on all four axes.** Properly bracketed, unlike every previous optimum in this project.
- Next best is **20× worse**.
- The clean search returns **theta = 0.50, hill = 4.0** — the same point the superseded
  destructive-probe scan chose.

**The selection contamination did not determine M1's shape.** An independent search over a grid
built from first principles rediscovers it. That specific worry is retired.

## Why this is not validation

`provenance.py` asserts it permanently: `recover_tau = 10080` cites *"TFI reinvigoration
(Philipp 2022)"*, the source of every calibration value. The held-out point reads −6.8 against a
±15 tolerance and **that pass is uninterpretable** — a favourable result on a target the model was
partly built from carries no information.

The search for a clean external target returned **none available**. Eight candidates, seven fixed
criteria, assessed before any prediction was computed. The only provenance-clean candidate fails
because M1 has no TCR-signalling term and an abstract carries no replicate structure.

## The failure that survives everything

M1 has **one** latent state. Philipp measured **two** functional axes, and both moved:

| | continuous | TFI | ratio |
|---|---|---|---|
| 72 h specific lysis (d14) | 34.9% | 93.4% | 2.7× |
| 3-day CD2+ expansion (d14) | 1.1× | 4.1× | 3.7× |

A single monotone mapping fits all four points, so the one-state hypothesis is **not falsified**.
It is also **barely tested**: monotone by chance in 1/4! = 4.2%; the points share a time trend; only
two time-controlled comparisons exist, both concordant, chance concordance **0.25**.

Parameters are sharply identified. The **mechanism is not**. Those are different things, and only
the second is what R and T would depend on.

## Honest limits on the Phase 1 result

- **Identified to grid resolution, not finer.** theta was sampled at 0.15/0.30/0.50/0.70/0.90 and
  hill at 1/2/4/8/16. A basin narrower than that spacing would look identical. One cell with a 20×
  gap is consistent with sharp identification *and* with a grid too coarse to resolve the basin.
- **Fitted to three numbers from one paper**, in one molecule, one target line, in vitro.
- **`recover_tau` was never fitted**, only shown live. Its value remains the contaminated default.

## Classification

**B — FIT BUT NOT VALIDATED**, unchanged, now with the selection-contamination objection removed
and the mechanism-identifiability failure quantified.

Not **D**: M1 reproduces the calibration within simulation noise on a bracketed grid.
Not **A**: no provenance-clean external observation exists.
Not parameter-**C**: exactly one parameterisation survives on this grid.
But **mechanism-C**: one latent state versus two is unresolved, with 0.25 chance concordance.

## What this licenses, and what it does not

**Licensed:** a threshold functional mapping reproduces Philipp's chronic-stimulation cytotoxicity
trajectory within simulation noise, at an interior optimum of an independently specified grid,
where a linear mapping cannot.

**Not licensed:** that the fitted state is exhaustion rather than a correlate; that the same state
explains expansion capacity; any chronic population dynamics; transfer to SCLC, to tarlatamab, to
in vivo, or to any schedule comparison. R and T remain blocked and `calibration.json` still reads
`FIT_NOT_VALIDATED`.
