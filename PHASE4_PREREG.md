# Phase 4 pre-registration

Written after three failed experiments (L1, Phase 2, Phase 3) and **before any Phase 4
tumour-dynamics result was generated or viewed**. Committed before the window sweep in §4 was run.

L1, Phase 2 and Phase 3 are preserved unaltered. Nothing here reinterprets them, and the
classification-A write-up of the previous line stands.

---

## 1. Why this is a new line and not a rescue

`PHASE3_PREREG.md` §4 closed the previous line after two failed positive controls, and forbade
modifying that model a third time. It is not being modified. This is the different model both
earlier pre-registrations explicitly contemplated:

> `PHASE2_PREREG.md` §2: *"within this model's structure, recruitment and population-level
> functional collapse cannot both be represented at a realistic influx rate, because arrivals enter
> with zero exhaustion."*

The distinction that matters: **the defect was identified by an external measurement, not by the
hypothesis needing to survive.** Blinatumomab is given by 28-day continuous intravenous infusion, so
the entire circulating T-cell pool is engager-exposed — not only cells inside the lesion. Philipp
Figure 1B assayed exactly that pool, peripheral T cells from r/r BCP-ALL patients on c.i.v.:

| | specific lysis | normalised to baseline |
|---|---|---|
| day 0, pre-treatment | 73.1% | 1.000 |
| day 14, on infusion | 17.4% | 0.238 |
| post-cessation | 48.5% | 0.663 |

**The cells available for recruitment are measured to be exhausted.** A model in which recruits
arrive naive is wrong about the recruits, as a matter of observation. That is the defect being
repaired, and it is the same defect that produced all three previous failures.

## 2. The change, and it is the only one

A scalar systemic compartment advancing on **drug exposure alone**, since a circulating T cell is
engager-exposed whether or not it is touching a blast:

```
drug present:  dEr_sys/dt = (1-rho_sys) k_sys (1-E_sys)      dEd_sys/dt = rho_sys k_sys (1-E_sys)
drug absent:   dEr_sys/dt = -Er_sys / tau_sys                 Ed_sys unchanged
```

Recruited T cells enter carrying `(Er_sys, Ed_sys)` and the compartment's accrued exposure, instead
of zeros. **Nothing else changes**: not the per-cell exhaustion clock, not `p_kill`, not `p_div`,
not `t_death`, not the schedules, endpoints, architectures, lattice, or `dt`. Both previous code
paths are retained and verified bit-identical, so L1, Phase 2 and Phase 3 all still reproduce.

## 3. Calibration (done, before this document was committed)

Fitted to the three Figure 1B values and **nothing else**. No schedule ranking, architecture
quantity or tumour burden was used as a target.

| Quantity | Status | Value |
|---|---|---|
| `k_sys` | **well identified** | **0.102 /day** (0.086–0.125), insensitive to post-timing |
| `rho_sys` | not separately identifiable | 0.00–0.43 |
| `tau_sys` | not separately identifiable | 0.5–16.6 d |
| **fraction of systemic exhaustion a 7-day break reverses** | **identified as a combination** | **0.34–0.72**, median ~0.60 |

Declared assumptions and uncertainties, stated rather than buried:

- The readout is treated as approximately **linear** in effector function over the observed 17–73%
  range, away from the assay ceiling. This differs from Phase B, where values sat at 88–93% and
  saturation had to be modelled explicitly.
- The **"post" sampling time is not stated in the paper**. Carried as a declared uncertainty over
  {day 35, day 42}. `k_sys` is unaffected (set by the day-14 point alone); `tau_sys` is not, so
  both are carried.
- **`rho_sys` is a separate parameter from the per-cell `frac_durable`, and this is forced.**
  Applying the per-cell value 0.93 to the systemic pool floors function at 0.123 against a measured
  0.663 — arithmetically impossible. It is also biologically expected: the circulating compartment
  recovers by **turnover** as well as by cells de-exhausting in place, a route no tissue-resident
  cell has.

`k_sys` at 0.64× the per-cell rate is a *consequence* of the fit, not an input. Circulating cells
receive tonic engager signal but less antigen-driven stimulation, so below 1.0 is the expected
direction, and it landing there is mild evidence the compartment is not mis-specified.

**Note against my own interest:** the systemic pool recovers ~60% of its exhaustion in a week
against ~4–10% for a tissue-resident cell. This gives treatment-free intervals a lever they did not
previously have, so the model is now *more* capable of showing a scheduling benefit. That is
precisely why §4 and §5 below are fixed in advance, and why the positive control still has to earn
its pass against tumour regrowth and lost drug time.

## 4. The window sweep (calibration, schedule-blind)

**Continuous dosing only.** No treatment-free interval, no schedule comparison, no architecture
comparison is computed here, so the selection cannot express a preference about schedule ranking.

Dispersed architecture, 42 days, 4 seeds. `N_T` ∈ {200, 400, 800, 1400} × influx ∈ {0, 2.5e-5,
1.0e-4, 4.0e-4}, per-cell representative `rec_low`, systemic variants for post ∈ {d35, d42}.

A regime is **admissible** only if it satisfies both, exactly as before:

1. **External check (unchanged from Phase C):** mean T-cell-pool function falls to ≤ 0.40 of its
   day-0 value by day 14, the criterion derived from Philipp Figure 1B.
2. **Dynamic range:** median `nB42` > 0; cleared ≤ 25% of seeds; median `nB42` < n0; peak burden
   < 70% of lattice capacity.

**The decisive question, and its answer is pre-committed either way:**

> Does any regime satisfy **both** tumour control (median `nB42` < n0) **and** binding exhaustion
> (mean pool function < 0.5)? In the previous model this was false in 24 of 24 cells and the window
> was empty.

**If the window is still empty: STOP.** Report that the systemic compartment does not open it, that
the empty-window finding is therefore robust to this repair, and classify **A**. Do not widen the
grid, do not adjust `k_sys`, do not try a fourth mechanism.

## 5. Positive control

Unchanged in every respect from `PHASE2_PREREG.md` §3 and `PHASE3_PREREG.md` §4 — same ten
schedules across both families, eight paired seeds, same primary endpoint, same floor and ceiling
contingencies, same statistics.

**PC-1 (required):** at least one non-continuous schedule reduces day-42 burden versus continuous
at one-sided Wilcoxon p < 0.05 **and** by a median of at least **10%**.

**PC-2 (reported, not required):** the winner should be a short repeating interval, off-fraction
≤ 3 days per cycle.

**This is the third and final model.** If PC-1 fails in every admissible regime, the project is
written up as a negative and no further mechanism is attempted. That is fixed now, before the
result exists, and it is not subject to revision by me after seeing it.

If PC-1 **passes**, the previously unreachable phases become available in the order already
pre-registered: architecture comparison (`PHASE2_PREREG.md` §6), then the trafficking attack with
its unchanged category-B decision rule (§7.1 and `TRAFFICKING_PREREG.md` T2), then exhaustion
knockout, dose matching, exhaustion distribution and the geometry continuum (§7.2–7.4), then
numerical robustness (§8). **A pass here does not license skipping any of them**, and the
trafficking attack remains the most likely thing to end the hypothesis.

## 6. Recorded, not acted on

- `p_div = 1/2880` per minute is a generic tumour cell-cycle time. Defensible for the leukaemia-like
  arm the positive control uses; **too fast for follicular lymphoma**, which is indolent. It must be
  faced before any architecture claim, and it confounds a design that holds all biology fixed while
  varying only geometry.
- The systemic compartment is **one-way**: it supplies the tissue but does not receive exhausted
  cells back from it. A limitation, not a claim.
- The systemic pool has no explicit size, so recruitment is never exhausted numerically. Also a
  limitation.

## 7. Standing constraints

All constraints from `PHASE2_PREREG.md` §11 remain in force: no calibration to any schedule ranking,
architecture difference or tumour burden; no new biological mechanism claimed; every pre-registered
endpoint reported regardless of outcome; all code, seeds, raw outputs and deviations preserved; no
cure claims and no clinical recommendations.
