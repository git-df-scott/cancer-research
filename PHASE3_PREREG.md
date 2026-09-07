# Phase 3 pre-registration

Written after the Phase 2 positive control failed (classification A) and **before any Phase 3
result was generated or viewed**. Committed before the sweep in §3 was run.

Experiments L1 and Phase 2 are preserved unaltered. Nothing here reinterprets either.

---

## 1. What failed, and what is being changed

The Phase 2 positive control failed on **dynamic range**, not on biology. Three of four admissible
combinations separated all ten schedule arms by under 1% of tumour burden, at 78–79% of lattice
occupancy. The continuous arm's burden flattened by day 21–28 because the domain filled up. See
`results/expE_diagnosis.txt`.

The arithmetic is unambiguous. At the initial burden of 5,541 cells with a 2-day cycle, the tumour
generates ~2,770 new cells/day. At the measured killing rate of ~3 targets/T-cell/day, holding that
steady requires **~924 T cells**. Experiment L1 had ~860 (with weak exhaustion, so a higher
effective rate) and **cleared everything — the floor**. Phase 2 has 10–142 and **loses to growth —
the ceiling**. Both experiments landed on opposite sides of a sharp control threshold.

**This is why enlarging the domain cannot fix it.** When killing is 5–100× below growth, any
starting burden runs away; a bigger lattice only delays saturation. The lever is the initial
effector-to-target ratio.

### The single change

**The initial T-cell number, `N_T`.** Nothing else. No biological rate, no exhaustion parameter, no
influx value, no schedule, no architecture definition, no endpoint.

`N_T` is an **initial condition** — an experimental design choice, like selecting an E:T ratio when
setting up an assay — not a biological parameter fitted to data.

### Why the existing value is externally wrong for a positive control

The positive control exists to reproduce a leukaemia-like result, and both reference systems it is
being compared against are far more effector-rich than 200 T cells against 5,541 targets (1:28):

| Reference system | effector:target |
|---|---|
| Philipp et al. 2022, the assay the exhaustion clock is calibrated to | **1:4** |
| Obertopp/Basanta, the model whose result is being reproduced | tumour seeded at 50% occupancy with **T cells already interspersed** |
| Phase 2 as run | **1:28** |

Matching Philipp's 1:4 would mean `N_T ≈ 1385`. That is close to the ~924 control threshold
computed above, which is *not* a coincidence — an assay is set up so that the effectors can
actually kill the targets, which is the same requirement as having dynamic range.

---

## 2. The safeguard that makes this not tuning

Selecting an operating point where an assay has dynamic range is standard practice. It becomes
tuning only if the selection can see the quantity under test. So:

> **The sweep in §3 uses the CONTINUOUS arm only. No treatment-free interval, no schedule
> comparison, no architecture comparison is computed, examined, or permitted to influence the
> choice of `N_T`.**

The window criterion below is a statement about whether the continuous arm sits between the floor
and the ceiling. It is structurally incapable of expressing a preference about schedule ranking,
because no non-continuous arm is simulated.

---

## 3. The sweep (calibration, not a hypothesis test)

Dispersed architecture, continuous dosing only, 42 days, 6 seeds.

- `N_T` ∈ {200, 400, 800, 1400, 2000, 2800} — 200 is Phase 2's value; 1400 is Philipp's 1:4;
  2800 is 1:2.
- influx ∈ {0, 2.5e-5} — the admissible set from Phase C, unchanged.
- exhaustion representative: all three (`rec_low`, `rec_mid`, `rec_high`), unchanged.

### Window criterion, fixed before running

`N_T` is **admissible** if, in the continuous arm:

1. median `nB42` **> 0** — the tumour is not eradicated (no floor), **and**
2. fraction of seeds reaching `nB = 0` **≤ 0.25**, **and**
3. median `nB42` **< n0** — the drug demonstrably does something, **and**
4. max burden at any point in the run **< 70% of lattice capacity** — never approaches the
   ceiling.

If several values qualify, take the one whose median `nB42` is closest to `n0 / 2`, i.e. the
centre of the available range. If none qualifies, **STOP and report that**; do not widen the grid
in search of a pass.

---

## 4. Phase 3 positive control

Identical to Phase 2's in every respect except `N_T`: same ten schedules across both families, same
eight paired seeds, same admissible exhaustion representatives and influx regimes, same primary
endpoint, same floor-effect contingency, same statistics.

**PC-1 (required), unchanged from `PHASE2_PREREG.md` §3:** at least one non-continuous schedule
reduces day-42 burden versus continuous at one-sided Wilcoxon p < 0.05 **and** by a median of at
least **10%**.

**Additional ceiling contingency, new and fixed now.** If the continuous arm's median `nB42`
exceeds 70% of lattice capacity in a given combination, that combination is reported as
**no dynamic range** and is excluded from the pass/fail determination rather than counted as a
failure. Phase 2 lacked this and so recorded a "failure" that was really an absence of measurement.

**If PC-1 fails in every combination that has dynamic range:** STOP. Report
`POSITIVE CONTROL FAILED` again, classification A. **Two failed positive controls under two
different pre-registrations is the end of this line of attack**, and the project is written up as a
negative rather than modified a third time.

---

## 5. What is NOT being changed, and one issue recorded for later

Unchanged: exhaustion submodel and its fitted parameters; the admissible influx set and the
external criterion that produced it; `p_kill`; `p_div`; `t_death`; the schedules; the endpoints;
the architecture definitions; lattice size; `dt`; the trafficking attack and its decision rule; the
stopping rules; the classification scheme.

**Recorded, not acted on.** `p_div = 1/2880` per minute — roughly 0.5 divisions per cell per day —
is cited in `README.md` as a generic "standard human tumour cell-cycle time". It is defensible for
the leukaemia-like dispersed arm, which is what the positive control uses. It is **too fast for
follicular lymphoma**, which is the archetypal indolent lymphoma, with grade 1–2 Ki-67 typically
well below the ≥50% that marks the adverse group even in grade 3A. If Phase 3's positive control
passes and the architecture comparison is reached, this must be addressed there, because the
follicle arms would be running at a leukaemia growth rate. It is **not** addressed here: changing
it now would alter a biological parameter, and the user's instruction for this phase was to fix
dynamic range only.

Fixing it later would also confound the design, since the project deliberately holds all biology
fixed and varies only geometry. Whether that is resolvable is itself an open question and is
recorded as such.

---

## 6. Standing constraints

All constraints from `PHASE2_PREREG.md` §11 remain in force: no calibration to any schedule
ranking, architecture difference or tumour burden; no new biological mechanism claimed; every
pre-registered endpoint reported regardless of outcome; all code, seeds, raw outputs and deviations
preserved; no cure claims and no clinical recommendations.
