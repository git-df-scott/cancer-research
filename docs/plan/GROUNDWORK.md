# Groundwork: what was built, and what it already established

Status 2026-09-08. Four new modules, all executable and self-tested. **No experiment has been run
to completion.** Everything below comes from validation runs and from the PK layer, which needs no
simulation. Read `docs/plan/PLAN.md` for the phase structure and `docs/recon/LUNG_CANCER_RECON.md` for why SCLC.

## What exists now

| File | What it is | State |
|---|---|---|
| `schedules.py` | Duty-cycle, continuous and legacy-L1 schedule generators, plus exposure accounting | Verified |
| `exp_replicate.py` | Experiment R: faithful reference replication + 4-factor ablation. Pre-registered R1–R4 | Runs; 280 jobs defined |
| `calibrate_exhaustion.py` | Philipp-assay replica (serial re-challenge) and 2-parameter fit harness | Runs; baseline measured |
| `pk.py` | Two-compartment PK for a half-life-extended engager | **10/10 self-tests pass** |
| `exp_tarlatamab.py` | Experiment T: PK-driven occupancy into the ABM. Pre-registered T1–T5 | Runs; T1/T5 already answered |

## Three things the groundwork settled before any experiment

### 1. The L1 schedule encoding is definitively wrong, and the exposure numbers prove it

`schedules.py` computes drug-time over the 42-day horizon:

| Arm | Exposure |
|---|---|
| Reference CONT (28 d dosing, 14 d rest, all arms) | 0.667 |
| L1 `tfi=0` (dosed all 42 days) | **1.000** |
| L1 `tfi=14` | 0.667 |

L1's continuous arm received **1.5× the drug-time the reference gives its continuous arm**, and
L1's *most* interrupted arm received exactly what the reference's *continuous* arm receives. The
reference's methods state TFIs are "recurring treatment-free intervals … alternating with periods
of TCE dosing" — a duty cycle. L1 encoded one interruption at the end of a cycle, i.e.
`TFI_k` = continuous minus k days, which is dominated by continuous by construction.

This is no longer an inference. It is arithmetic, and it is reproducible from `schedules.py`.

### 2. But fixing the schedule and E:T is *not* sufficient — and that partly vindicates `docs/plan/HANDOFF.md`

A validation run of the faithful reference configuration (1:4 E:T, 50% occupancy, no influx,
proliferation on, duty-cycle schedules):

```
ref  ('cont',)       n0=9600  nB16=0  nB28=0  nB42=0  meanE16=0.049
ref  ('duty',5,2)    n0=9600  nB16=1  nB28=0  nB42=0  meanE16=0.028
```

**The tumour is eradicated before day 16 and every endpoint sits on the floor.** With exhaustion
this weak, killing never collapses, so no schedule can be distinguished from any other. The
reference gets dynamic range precisely *because* its T cells exhaust hard enough (>80% by day 16)
for killing to fail and the tumour to regrow.

So the correction to `docs/findings/FINDINGS.md` is now two-sided, and I got part of my earlier read wrong:

- The schedule encoding and E:T ratio **are** misconfigured, decisively, and L1's ranking is
  substantially an artefact of that. `docs/findings/FINDINGS.md`'s framing of L1 as a clean test is wrong.
- **Exhaustion recalibration is genuinely required as well.** `docs/plan/HANDOFF.md` was right that it
  blocks. My earlier claim that it was probably *not* the cause was too strong.

Both are true at once. Phase 1 is not optional.

### 3. The exhaustion rate is not wrong — the *mapping* to the assay is

`lymphoid.py` sets `exhaust_tonic = 2.48e-5`/min, documented as "28 d contact → E=1". That
arithmetic checks out against the external target: Philipp's d28/d7 retention of 0.097 implies
E(d28) ≈ 0.903; the model gives 1.000 under continuous contact.

But run the model **in Philipp's actual assay conditions** — serial re-challenge, targets
replenished so contact never lapses — and the shipped parameters give:

| | model | Philipp | error |
|---|---|---|---|
| d7 continuous | 87.1 | 88.4 | −1.3 |
| d14 continuous | 74.0 | 34.9 | **+39.1** |
| d28 continuous | 62.5 | 8.6 | **+53.9** |
| d14 with TFI | 96.4 | 93.4 | +3.0 |

Decay is far too slow; **recovery is already right**. That is a much more precise diagnosis than
"exhaustion is miscalibrated", and it says the fit should move the decay parameters and leave
`recover_tau` alone.

`calibrate_exhaustion.py` also identifies a candidate mechanism rather than just a rate. The
reference accrues exhaustion on *killing events* ("T cells accumulate PD-1 during tumor killing…
at higher rates than baseline"); `lymphoid.py` sets `exhaust_per_kill = 0.0` and accrues on dwell
time only. Dwell-driven exhaustion is linear in elapsed contact and takes weeks; kill-driven
exhaustion ramps hardest exactly when targets are abundant. Both parameters are fitted jointly
against the same four external targets.

## The lung result that needs no simulation

`pk.py` is validated against published quantities (terminal half-life at both 5.8 d and 11.2 d,
AUC = dose/CL, C(0) = dose/V₁, dose proportionality, superposition, and the feasibility bound on
V₂/V₁). Steady-state occupancy under the clinical 10 mg Q2W regimen:

| EC50 | drop from peak, t½ 5.8 d | t½ 11.2 d | T1 (<25%) |
|---|---|---|---|
| 0.1 nM | 2.9% | 2.1% | holds |
| 0.3 nM | 8.1% | 6.0% | holds |
| **1.0 nM** | **22.3%** | **17.3%** | **holds, narrowly** |
| 3.0 nM | 44.7% | 37.1% | fails |
| 10 nM | 69.1% | 62.1% | fails |
| 100 nM | 87.4% | 83.8% | fails |

**The hypothesis holds at the in-vitro-anchored EC50 (~1 nM) and dies between 1 and 3 nM.**

That anchor is the weakest link and must not be oversold: ~1 nM is an *in vitro co-culture* EC50
(PMC13082099), which is a different measurement in a different system from in vivo potency. Any
eventual claim is at most **classification D, conditional on EC50**, with the conditionality in
the claim itself rather than a footnote.

What the ABM still has to add: occupancy is a statement about *drug*, not about whether T cells
are *in contact* while drug is present — and contact is what exhaustion integrates. In spatially
structured disease a T cell can be drug-saturated and unengaged. T2/T3 are not implied by T1.

## Known deviations, recorded not hidden

- **Grid.** Reference is 150×160 = 24,000 sites; `Lymphoid` is square, so L=155 (24,025) is used.
- **T-cell division.** Reference gates on a 5 h refractory clock per cell; `Lymphoid` uses a
  per-minute hazard gated on drug and target adjacency. An approximation, not an equivalence.
- **Exhaustion mechanism.** Reference uses a hard PD-1 threshold plus death; `Lymphoid` uses a
  continuous scalar scaling kill hazard. R2 is the pre-registered check on whether this matters.
- **V₂ and Q are not published.** `pk.py` solves Q analytically for the target terminal half-life
  given an explicit V₂/V₁ ratio, and the feasible range is computed, not assumed — at t½ 5.8 d the
  ratio must be < 0.579, at 11.2 d < 2.048.
- **Validation runs were deleted** rather than kept, because they were smoke tests at n=1 and
  keeping them would imply a result. All real runs write per-job files and are preserved.

## Handoff: what to run, in what order

1. `.venv/bin/python calibrate_exhaustion.py` — 32-cell grid over `exhaust_tonic` ×
   `exhaust_per_kill` × `recover_tau`. **Accept on the Philipp loss alone; never on a ranking.**
   Report all four targets for the chosen cell, not just the two that fit.
2. `.venv/bin/python exp_replicate.py rep` — R1/R2/R4 under the recalibrated parameters.
   **Gate: if R1 fails here, stop and report.** Everything downstream is meaningless otherwise.
3. `.venv/bin/python exp_replicate.py abl` — the 4-factor attribution table. Run this regardless
   of whether R1 passes; it is the deliverable that L1 never produced.
4. `.venv/bin/python exp_tarlatamab.py run` — T2/T3/T4 in the ABM, at EC50 1 nM and 10 nM
   (either side of the boundary).
5. Only then: correct `docs/findings/FINDINGS.md` with the ablation as evidence, and classify.

Runtime is ~50 s per 42-day run at dt=5. Experiment R is 280 runs; T is 80 at 84 days. Every job
writes its own file and skips if present, so a stall costs one run, not the batch.

## Standing constraints

No parameter fishing. Calibrate to external curves, never to rankings. Report every pre-registered
endpoint regardless of outcome. Preserve seeds, raw outputs and deviations. Final classification is
exactly one of A–E. No cure claim, no clinical recommendation, no proposed dosing change.
