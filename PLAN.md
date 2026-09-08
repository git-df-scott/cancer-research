# Plan: fix the control, then ask the lung question

Written 2026-09-08, after the lung recon (`LUNG_CANCER_RECON.md`) and a close read of the
reference model. Two tracks run in parallel by decision of the project owner.

---

## What planning turned up, and why it changes everything

`FINDINGS.md` and `HANDOFF.md` both attribute the failed positive control to **exhaustion
miscalibration**, and set recalibration as the blocking first task. Reading the reference model's
published methods ([PMC12667981](https://pmc.ncbi.nlm.nih.gov/articles/PMC12667981/)) against
`exp_schedule.py` shows something different and much cheaper to fix.

**The dispersed arm was never a replication of the reference.** It differs in three structural
ways, any one of which could produce the observed ranking.

| | Reference (Obertopp/Basanta) | This project's `exp_schedule.py` | Impact |
|---|---|---|---|
| **Schedule shape** | "recurring treatment-free intervals (2–7 days) **alternating with periods of TCE dosing**" | `schedule_fn`: one 28-day cycle, on for `(28 − tfi)` days, then off for `tfi` days | **Decisive** |
| **E:T ratio** | 1:4 (2,400 T : 9,600 tumour) | 1:27.7 (200 T : 5,542 B) | **Decisive** |
| **T-cell proliferation** | on — active T cells divide after killing, 5 h refractory | `t_div=0.0`, disabled | Large |
| Occupancy | 50% | 95% (dispersed arm) | Moderate |
| T-cell influx | none | `t_influx=1e-4`, 200 → ~860 over a run | Large |

### Why the schedule encoding alone can produce the result

Under `schedule_fn`, `TFI_k` means *continuous dosing minus the last k days*. It is strictly
dominated by continuous by construction unless exhaustion is severe enough that the recovery pays
for the lost drug-time. That guarantees the monotone ordering actually observed —
`TFI0 > TFI2 > TFI4 > TFI7 > TFI14` — in **all three architectures**, which is exactly what L1
reported. The reference's `TFI_k` is a *duty cycle*: repeated short off-periods interleaved
throughout treatment, which lets T cells recover many times, not once at the end.

These are not the same regimen, and the model was never asked the reference's question.

### Why the E:T ratio compounds it

At 1:27.7 the system is target-rich and effector-poor: every T cell has a target in contact
essentially always, killing is effector-limited, and losing 25% of drug-time costs ~25% of
killing directly. Exhaustion cannot become the binding constraint in that regime no matter how it
is calibrated. At the reference's 1:4 it can.

**Consequence for the plan:** exhaustion recalibration may still be needed, but it is no longer
the *first* thing to try, and `FINDINGS.md`'s root-cause attribution is probably wrong. Phase 0
settles it in hours rather than weeks.

**Honest caveat.** The reference does not tabulate the ON-days between its OFF-intervals, so the
duty cycle must be assumed. Phase 0 states the assumption explicitly and sweeps it (see A3).

---

## Track A — Phase 0: is the model actually broken?

**Question.** Does the existing model reproduce short-TFI-beats-continuous when configured as the
reference actually configured it?

**New file:** `exp_replicate.py`. Nothing existing is modified.

### A1. Faithful replication arm

Reference conditions, matched: 150×160 lattice, 50% occupancy, 1:4 E:T (2,400 T : 9,600 tumour),
random placement, **no influx**, **T-cell proliferation on**, 42-day horizon (28 dosing + 14 rest).

### A2. Correct schedule encoding

Replace the single-interruption `schedule_fn` with a duty-cycle generator:

```
tfi_cycle(on_days, off_days) -> drug is 1.0 for on_days, 0.0 for off_days, repeating
```

`CONT` = always on for days 1–28. Every arm rests days 29–42.

### A3. The duty-cycle assumption, swept not guessed

The reference does not state ON-days. Run `on_days ∈ {2, 3, 5, 7}` crossed with
`off_days ∈ {2, 3, 7}` plus `CONT`. Report the whole grid. If the ranking is stable across the
sweep the assumption does not matter; if it flips, say exactly where.

### A4. Endpoints, pre-registered before running

- Primary: median tumour burden at day 28 and day 42, per arm.
- Secondary: mean exhaustion trajectory; peak exhaustion under `CONT`; T-cell count.
- Diagnostic: does `CONT` exceed 80% exhaustion by day 16, as the reference reports?

### A5. Decision gates

| Outcome | Meaning | Next |
|---|---|---|
| Some short TFI beats CONT at day 28 **and** CONT exhaustion >80% by d16 | Control passes. The model was fine; the experiment was misconfigured. | Correct `FINDINGS.md`, re-run L1 properly, proceed |
| Ranking correct but exhaustion far below 80% | Schedule was the bug, calibration still off | Phase 1 recalibration, narrowed |
| CONT still wins under faithful conditions | Genuine model failure | Phase 1 recalibration, full scope |

**Ablation, run regardless:** toggle each of the four differences (schedule, E:T, proliferation,
influx) one at a time from the L1 configuration toward the reference, to attribute the failure
quantitatively rather than asserting it. This is the deliverable even if the control passes.

**Cost:** ~15 configurations × 10 seeds. At the recorded ~54 s per 42-day run at dt=5, roughly
2–3 core-hours; parallelisable. Per-result file writes, per `rerun_missing.py`, so a stall costs
one run.

---

## Track B — Phase 3: the tarlatamab PK question

Runs in parallel. It is a **schedule** question, so it does not depend on the FL *architecture*
result — but it does depend on exhaustion being calibrated, so its final answer is gated on
Track A. Build and validate now; interpret after.

**New files:** `pk.py`, `exp_tarlatamab.py`.

### B1. The PK layer

Two-compartment linear elimination, from the published population PK
([PMID 40261494](https://pubmed.ncbi.nlm.nih.gov/40261494/)): CL 0.649 L/day, Vc 3.44 L for a
73 kg subject, terminal half-life 5.8 d (DeLLphi-300) to 11.2 d (popPK median).

Replace the binary `self.drug` with a concentration-driven occupancy term. `drug` becomes
`C(t) / (C(t) + EC50)`, an Emax form, rather than a switch. Killing hazard and exhaustion accrual
both read that occupancy, as they already read `self.drug`.

**Unit-test the PK independently** against the published half-life and dose-proportionality
before it is wired into the ABM. A PK bug would silently invalidate everything downstream.

### B2. Schedules compared

- Clinical: 1 mg step-up C1D1, 10 mg C1D8, C1D15, then 10 mg Q2W
- Continuous-equivalent infusion at matched AUC
- Fractionated: same total exposure, lower peak, engineered trough (e.g. 5 mg weekly)
- Genuine drug holidays: Q2W with a real off-cycle, at matched and unmatched exposure

### B3. The pre-registered prediction, written before running

> With a terminal half-life of 5.8–11.2 days on a 14-day interval, simulated trough occupancy
> never falls low enough for engaged fraction to drop materially between doses. Cumulative
> effector dwell time under the clinical schedule will therefore be within 10% of matched-AUC
> continuous infusion, and exhaustion at day 28 will be statistically indistinguishable between
> them.

**This dies cleanly.** If engaged fraction between doses falls by more than 25% from peak, the
hypothesis is false and gets reported false. No reinterpretation.

### B4. Sensitivity that must be run, not skipped

The conclusion depends on EC50, which is **not** established in vivo — the available anchor is an
in vitro co-culture EC50 of ~1 nM (PMC13082099), which is not the same quantity. Sweep EC50 across
three orders of magnitude and report the fraction of that space where the prediction holds, plus
the boundary where it flips. If the result only holds in a narrow EC50 window, say so and classify
it as parameter-dependent.

**Prior-art position, already audited:** population PK models are PK-only; QSP CD3-bispecific
models are non-spatial and not tarlatamab; the reference ABM used a stand-in half-life and flagged
it as a limitation. The gap is the combination, and only the combination is claimed.

---

## Later phases, gated

- **Phase 1 — exhaustion recalibration.** Only if Track A says it is needed. Method fixed in
  advance to avoid circularity: calibrate in an **in-vitro-replica configuration** (no influx,
  fixed E:T, well mixed) against Philipp's measured curve, then transplant the parameter unchanged
  into the in vivo configuration. Calibrate to the external curve, never to a schedule ranking.
- **Phase 2 — re-run L1** with correct schedule encoding and report every pre-registered endpoint.
- **Phase 4 — heritable per-cell antigen.** The highest-value structural change and the least
  data-constrained; treat with corresponding suspicion. Enables the antigen-low selection question.
- **Phase 5 — the attack battery** (`exp_attacks.py`, already written and pre-registered):
  exhaustion knockout, dose-matched controls, trafficking knockout.
- **Phase 6 — SCLC geometry.** Blocked. Requires nest dimensions from a retrievable source; must
  not be simulated on guessed geometry.

---

## Standing constraints, carried from the existing project

- Report every pre-registered endpoint regardless of outcome. No reinterpretation after seeing
  results. No parameter fishing toward a positive result.
- Preserve all code, seeds, raw outputs, pre-registrations, failed hypotheses and deviations.
- Final classification must be exactly one of: A hypothesis falsified, B model artefact, C known
  result, D robust computational prediction, E independently supported prediction.
- No cure claim, no clinical recommendation, no proposed change to tarlatamab dosing.
- A clean negative is a success. So is establishing that the L1 failure was a configuration bug.

## Immediate next actions

1. `exp_replicate.py` — replication rig, duty-cycle schedule generator, pre-registration in the
   docstring written before any run.
2. `pk.py` — two-compartment model, unit-tested against published half-life in isolation.
3. Run Phase 0 ablation, report the attribution table.
4. Correct `FINDINGS.md` once Phase 0 returns, with the ablation as evidence rather than assertion.
