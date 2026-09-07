# Phase 2 pre-registration

Written after calibration (Phases A–C) and **before any Phase 2 architecture result was
generated or viewed**. Committed before the positive control was run. Experiment L1 is preserved
untouched and remains bit-identical reproducible; nothing here reinterprets it.

Everything below is fixed in advance: parameter families, admissible ranges, the positive-control
pass criterion and its minimum effect size, endpoints, the attacks, and the stopping rules.

---

## 0. What changed since L1, and why

L1 could not test the hypothesis because its positive control failed. Two independent causes were
found, and both are fixed here. Neither fix was chosen by looking at a schedule ranking.

**Cause 1 — the exhaustion clock was mis-specified, not just mis-tuned.** L1 reached at most 12%
functional loss at day 28 against ~90% measured. Retrieving the full text of Philipp et al. (Blood
2022, PMID 35878001) yielded four quantitative readouts the original handoff did not have, two of
which force structural changes:

- The TFI arm does not return to its day-14 level after its second rest (58.7% vs 93.4% specific
  lysis). Recovery is **partial**, and the shortfall grows with cumulative exposure. A single
  fully-reversible exhaustion state cannot represent that.
- Granzyme B MFI is a per-cell readout and does not saturate the way % specific lysis does. Its
  CONT/TFI ratio pins the **scale** of remaining function, which the lysis percentages alone
  cannot. Without it the fit drives the assay-saturation constant to 10.0 and declares every arm
  functionally dead while still reproducing all five lysis points.

**Cause 2 — L1 tested the wrong schedule family.** The reference result being reproduced
(Obertopp/Basanta, bioRxiv 2025, PPR1121269) states: *"a 7-day TFI improved T-cell function over
continuous dosing during the initial 28-day treatment period. However, when simulations were
extended to a full 42-day cycle... this advantage was lost. In contrast, shorter TFIs consistently
outperformed both 7-day and continuous schedules... A translationally oriented Monday-through-Friday
(MO_FR) regimen also achieved comparable benefits."*

L1 only ever tested a **28-day cycle with the last k days off**. In that family a "2-day TFI" is a
single 2-day break in 28 days: a 7% dose reduction with one recovery opportunity. The reference
effect lives in **repeating** short breaks — the MO_FR regimen is 5 days on, 2 days off, every
week. L1 never tested that family at all.

Note that L1 was *not* wrong about the 7-day interval: it found TFI-7 loses by day 42, and so does
the reference. L1 simply never sampled the schedules where the effect exists.

---

## 1. Calibrated parameter families and accepted uncertainty

### 1.1 The exhaustion submodel (fixed structure, no further states)

Per T cell, while in antigen contact with engager present:

```
dEr/dt = (1 - rho) * k * h(C) * (1 - E)      reversible component
dEd/dt =      rho  * k * h(C) * (1 - E)      durable component
dC/dt  = 1                                    cumulative engaged exposure
```

while the engager is absent: `dEr/dt = -Er / tau_r`, and `Ed` is unchanged. `E = Er + Ed`,
remaining function `f = 1 - E`, kill hazard proportional to `f`. `h(C)` is an optional Hill lag on
the cell's own cumulative engaged exposure, Hill exponent fixed at 4.

Observation model for the assay only: predicted specific lysis `= 1 - exp(-a * f)`, one shared
constant `a`, first-order kill kinetics in a fixed-duration fixed-E:T assay.

**No further states will be added.** The pre-registered rule that the Hill lag is adopted only if
it more than halves SSE was applied and the lag was **not** adopted (414.0 → 414.0).

### 1.2 What the data identify, and what they do not

Grid search over 244,800 parameter sets, scored against seven measured points and nothing else.

| Quantity | Status | Value |
|---|---|---|
| `k`, accrual rate on remaining function | **identified** | 0.159 /day, at every point of both profiles |
| `rho`, durable share | **not identifiable** | SSE flat from 0.00 to 0.93 |
| `tau_r`, reversible recovery | **not identifiable** | ranges over three orders of magnitude |
| fraction of accrued exhaustion a 7-day break reverses | **identified as a combination** | admissible 4.3% – 55%, fit prefers the low end |
| `a`, assay saturation | weakly identified | 1.8 – 13.8, median 2.8 |

`rho` and `tau_r` are therefore **not** reported as a single best value. The identified combination
is carried forward as the sensitivity axis.

### 1.3 The three exhaustion representatives (sensitivity axis, none preferred)

All share `k = 0.159/day`. All are run in every architecture comparison.

| Label | 7-day break reverses | rho | tau_r | C50 | a | SSE |
|---|---|---|---|---|---|---|
| `rec_low` | 4.3% | 0.93 | 7.49 d | 0 | 5.85 | 420.7 |
| `rec_mid` | 30.6% | 0.69 | 1.60 d | 4 d | 2.57 | 543.6 |
| `rec_high` | 55.0% | 0.45 | 0.10 d | 6 d | 1.89 | 798.8 |

**Recorded deviation.** The pre-registered absolute acceptance region (RMSE ≤ 5 percentage points)
is empty: the best achievable is 7.69 pp. The reason is in the data, not the model — the three
independent day-28 readouts of the *same* CONT/TFI functional ratio span a 10-fold range
(proliferation 2.1, lysis-derived 14.7, granzyme B 23.2), and a model with one functional axis
predicts a single ratio for all three. The class is therefore defined relative to the best fit,
`SSE ≤ 2 × min`. This deviation is recorded rather than hidden, and the direction of the residual
is stated in §7.

### 1.4 Held-out check (not used in fitting)

CD2+ proliferation fold-change ratio, CONT/TFI at day 14: measured 26.8, model 28.6.

---

## 2. Admissible T-cell influx regimes (Phase C)

External anchor, independent of the in vitro calibration: Philipp Figure 1B, patients on
continuous blinatumomab infusion, ex vivo specific lysis of **peripheral** T cells fell from 73.1%
(day 0) to 17.4% (day 14) — to 0.238 of baseline — in a compartment where recruitment,
redistribution and engager-driven expansion are all occurring. Population-level function collapses
in vivo despite influx.

**Pre-registered admissibility criterion:** an influx regime is admissible only if, under
continuous dosing in the dispersed architecture, mean T-cell-pool functional capacity falls to
≤ 0.40 of its day-0 value by day 14. Threshold set at 0.40 rather than the measured 0.238 to be
deliberately generous to the model.

This criterion concerns reproduction of an external measurement. It is **not permitted to be
adjusted by any treatment-schedule ranking.**

### Result (Phase C, run before this document was committed; no architecture run existed)

Dispersed architecture, continuous dosing, 4 seeds, 28 days. Mean T-cell-pool functional capacity:

| influx | `rec_low` f(d14) | `rec_mid` f(d14) | `rec_high` f(d14) | T pool at d28 | verdict |
|---|---|---|---|---|---|
| 0 | 0.122 | 0.264 | 0.370 | ~25 | admissible for all three |
| 2.5e-5 | 0.343 | 0.580 | 0.665 | ~160 | admissible for `rec_low` only |
| **1.0e-4** (the value used in L1) | 0.533 | 0.828 | 0.901 | ~750 | **excluded for all three** |
| 4.0e-4 | 0.873 | 0.990 | 0.997 | — | excluded for all three |
| 1.6e-3 | — | — | — | — | excluded for all three |

**The influx rate used in experiment L1 is incompatible with the external patient data.** It holds
population functional capacity at 0.53–0.90 at day 14 against a measured 0.238. L1's own diagnosis —
that influx dilution was preventing population-level functional collapse — is confirmed and now
quantified against an external measurement rather than by judgement.

**Honest caveat, recorded now.** Zero influx passes the functional criterion but carries its own
unrealism: with no replenishment the T-cell pool falls from 200 to about 25 by day 28 through
background death, which is not what happens in patients. The single most physiologically defensible
admissible combination is therefore `rec_low` at influx 2.5e-5, which is the only one that both
reproduces the functional collapse and holds the pool roughly stable (160 at day 28).

**Admissible set carried forward** (per-representative union, not the intersection, for the reason
above): `rec_low`×0, `rec_low`×2.5e-5, `rec_mid`×0, `rec_high`×0. Four combinations.

**Structural limitation, recorded.** Within this model's structure, recruitment and population-level
functional collapse cannot both be represented at a realistic influx rate, because arrivals enter
with zero exhaustion. A variant in which arrivals are not naive is listed in §8 as a robustness
check. It is **not** permitted to be introduced as a rescue if the positive control fails.

---

## 3. Positive control (Phase E). Dispersed / leukaemia-like architecture only

The model must reproduce the qualitative fact that an appropriate short, **repeating**
treatment-free interval outperforms continuous exposure.

**PC-1 (primary, required).** In the dispersed architecture, at least one non-continuous schedule
gives lower malignant burden at the primary endpoint than continuous dosing, by a one-sided
Wilcoxon signed-rank test over paired seeds at p < 0.05, **and** with a median paired reduction of
at least **10%** of the continuous arm's burden.

The 10% floor is stated now so that a statistically significant but negligible difference cannot be
counted as a pass.

**PC-2 (reported, not required).** The winning schedule should be a short repeating interval
(off-fraction ≤ 3 days per cycle), matching the reference. A pass driven only by a long interval is
flagged as anomalous.

**Scope rule.** PC-1 must pass for at least one admissible (exhaustion representative × influx
regime) combination. Every combination is reported regardless of outcome. The architecture
comparison is then run in **every** combination where PC-1 passes, and any architecture conclusion
must hold across all of them to be reported as a result.

**If PC-1 fails in every admissible combination:** STOP. Do not run the architecture comparison.
Report `POSITIVE CONTROL FAILED — MODEL NOT VALIDATED FOR SCHEDULING QUESTION`, classification A.
Diagnose the cause. Any subsequent modification requires new external calibration and a new
pre-registration.

---

## 4. Schedules

**Family A** — as in L1, for continuity: 28-day cycle, last *k* days off, *k* ∈ {0, 2, 4, 7, 14}.

**Family B** — repeating cycles, the family in which the reference effect exists.
(on, off) days ∈ {(6,1), (5,2) = MO_FR, (4,3), (12,2), (7,7)}.

`k = 0` in Family A is continuous dosing and is the shared comparator for both families.

---

## 5. Endpoints

**Primary:** malignant burden `nB` at day 42, matching the reference's endpoint and L1's.

**Pre-registered contingency for the floor effect observed in L1:** if, in a given architecture and
parameter combination, more than half the seeds in the *continuous* arm reach `nB = 0` before day
42, the primary endpoint for that combination switches to **time to clearance**, with non-clearing
runs censored at day 42 and arms compared by log-rank. This is fixed now because L1 reached medians
of 4 and 7 cells out of 5,525 and flagged the loss of dynamic range.

**Secondary, all reported regardless of outcome:** `nB` at day 28; integrated burden ∫nB dt over
0–42 d; peak burden; cumulative kills; total engager exposure delivered; T-cell abundance;
instantaneous engaged fraction; cumulative engagement time per T cell; mean `Er`; mean `Ed`; mean
functional capacity; radial position of T cells; exhaustion of front-line vs reservoir T cells;
clearance fraction; full schedule ranking.

---

## 6. Architecture hypotheses (Phase F)

Architectures: `dispersed` (leukaemia-like), `multi` (four follicles), `follicle` (one dense
follicle). Total malignant and T-cell numbers identical across architectures; only geometry differs.
Paired by seed.

- **H1** Engaged fraction under continuous dosing is lower in `follicle` than `dispersed`, with
  `multi` intermediate. *(Held in L1 at 4.7×. Not novel; expected.)*
- **H2** Mean exhaustion at day 28 under continuous dosing is lower in `follicle` than `dispersed`.
  *(Held in L1 at 4.8×.)*
- **H3 — the claim.** The *value* of the best treatment-free interval, relative to continuous
  dosing, is smaller in `follicle` than in `dispersed`. Tested as an architecture × schedule
  interaction.
- **H4** The best interval *duration* differs between architectures.
- **H5** Across the geometry continuum (Phase H), the TFI benefit is predicted by a dimensionless
  quantity — candidates fixed now: mean engaged fraction; mean cumulative contact time per T cell
  per unit treatment time; the ratio of contact-driven exhaustion timescale to recovery timescale;
  surface-to-volume ratio of the malignant compartment.

**Falsification, stated before running.** If the TFI benefit is equal across architectures, H3 is
wrong. If it is *larger* in `follicle`, the sign is opposite to the mechanism and H3 is wrong. Both
readings are recorded here in advance.

Multiplicity: H3 is tested once per (representative × influx × schedule family) combination.
Reported p-values are accompanied by Holm correction across schedules within a combination. L1's
single surviving signal was one comparison of four with no correction, and was reported as such.

---

## 7. Attacks, all mandatory, all pre-registered

### 7.1 Trafficking (Phase G) — the decisive one

`lymphoid.py`'s `swap_prob` lets a T cell exchange places with a malignant B cell instead of
requiring a vacancy. `TRAFFICKING_PREREG.md` recorded, before any alternative was run, that
absolute volume exclusion is *probably wrong*: the 11 µm/min T-cell speed that calibrates the model
was itself measured in densely cellular lymph-node cortex. It is already established that at
`swap_prob = 0.5` the follicle engaged fraction rises from 0.048 to 0.576 — essentially the
dispersed value — with T cells reaching 93% of the way to the follicle centre **without any
killing**.

Run the full critical comparison at `swap_prob` ∈ {0.0, 0.5, 1.0}.

**Decision rule, unchanged from `TRAFFICKING_PREREG.md`:** if the architecture × schedule
interaction shrinks by more than half at `swap_prob = 0.5`, classify as **model artefact
(category B)** and do not present it as an FL prediction. Do not select the trafficking rule that
preserves the hypothesis.

### 7.2 Exhaustion knockout (Phase H)

Set accrual to zero, holding architecture and delivered exposure fixed. If the architecture ×
schedule interaction survives essentially unchanged, the proposed chain
architecture → contact → exhaustion → schedule is **false**, and the actual cause must be found and
reported.

### 7.3 Dose matching (Phase H)

Kill hazard is linear in `drug`, so a schedule with duty cycle *d* is matched by a continuous arm
at level *d*. If a treatment-free interval does not beat its own dose-matched twin, there is no
timing effect, only dose.

**This control has already been run experimentally.** Philipp's TFI arm at day 21 has had exactly
14 days of cumulative exposure, the same as the continuous arm at day 14; measured after its rest it
is at 58.7% specific lysis against 34.9%. A genuine timing effect exists in the real system. The
calibrated model reproduces it at about 1.2–1.5× in function, against ~1.4× in granzyme B — i.e.
the model **under-predicts** the matched-exposure timing benefit. That is the conservative
direction for this project's hypothesis and is recorded here so it cannot later be presented as a
strength.

### 7.4 Exhaustion distribution

Report exhaustion among engaged vs unengaged T cells, and by radius. If the mean is low while the
front is fully exhausted and continuously replaced from a fresh reservoir, that is a different
mechanism and the claim must be **restated, not preserved**.

---

## 8. Numerical robustness (Phase I)

Required before any surviving effect is believed: timestep (`dt = 5` vs `dt = 1`, on the specific
quantities carrying the final claim, not on generic ones), lattice resolution, domain size, seed
count and stochastic uncertainty. L1's note that "the discretisation is common-mode across arms and
cancels in paired comparisons" is **not** assumed here and must be demonstrated for the quantities
that carry the claim.

---

## 9. Stopping rules

1. Positive control fails in every admissible combination → **STOP**, classification A.
2. Interaction halves at `swap_prob = 0.5` → **STOP** pursuing it as biology, classification B.
3. Interaction survives exhaustion knockout → the proposed chain is false; report that, and
   identify the true cause before anything else.
4. No treatment-free interval beats its dose-matched twin → no timing effect; report dose, not
   timing.
5. Any claim not surviving §8 is not reported as a result.

No result is reported at a strength beyond what these rules permit.

---

## 10. Success and failure

**Success** is establishing one correct thing that was not previously established, **or**
rigorously demonstrating why the idea fails. A clean negative is a success. An artefact diagnosis is
a success. A failed positive control, honestly reported, is a success.

**Failure** is a result obtained by tuning toward a desired answer.

Final classification is exactly one of:

- **A** — model not validated (positive control failed)
- **B** — model artefact (depends on an implausible trafficking / numerical / structural assumption)
- **C** — hypothesis falsified (validated model, controls passed, architecture produces no
  scheduling consequence)
- **D** — known result
- **E** — robust computational prediction
- **F** — independently supported prediction (frozen before inspection, later matched by
  independent evidence)

---

## 11. Standing constraints

- No further general literature audit. Phase 8 is closed; see `PRIOR_ART.md`.
- No parameter fishing. No calibration to any schedule ranking, architecture difference, tumour
  burden, or desired conclusion.
- No new biological mechanism is claimed. Every individual link has prior art. A conjunction of
  known mechanisms is not a new mechanism.
- Every pre-registered endpoint is reported regardless of outcome. No hypothesis is reinterpreted
  after seeing results.
- All code, seeds, raw outputs, calibration fits, pre-registrations, deviations and failed
  hypotheses are preserved.
- No cure claims. No clinical recommendations. Simulated days under stated rate calibrations are
  not clinical predictions.
