# Task board — Astra

Five independent tasks. **None touches the calibration/replication critical path**, which is held
elsewhere: do not modify `calibrate_exhaustion.py`, `run_calib_scan.py`, `exp_replicate.py`,
`schedules.py`, or anything under `results/calib/` or `results/expR/`. Everything below is either a
new file or a document.

Read `docs/plan/GROUNDWORK.md` first — it states what is built, what was already established, and the three
corrections to the record. `docs/plan/PLAN.md` has the phase structure. Standing constraints at the bottom of
this file are not optional.

Ranked by value. A1 is the one that decides whether the lung claim survives at all.

---

## A1 — Find an in vivo potency anchor for tarlatamab (highest value, and it may kill the project)

**Why this matters more than anything else on the board.** The lung hypothesis is that a
half-life-extended engager on a Q2W interval never lets target occupancy fall, so exhaustion
accrues as if dosing were continuous. `pk.py` shows that claim is entirely EC50-dependent:

| EC50 | occupancy drop per inter-dose interval | verdict |
|---|---|---|
| 1.0 nM | 17–22% | holds, narrowly (threshold is 25%) |
| 3.0 nM | 37–45% | **fails** |
| 10 nM | 62–69% | fails |

The only anchor currently in hand is an **in vitro co-culture EC50 of ~1 nM** (PMC13082099). That
is a different measurement in a different system, and using it as in vivo potency is an assumption
we are currently carrying uncosted.

**Do:** find any of the following, in order of preference.

1. Receptor-occupancy or target-engagement data for tarlatamab in patients or in vivo models.
2. Exposure–response analysis from DeLLphi-300/301/303/304 — anything relating serum concentration
   to response, CRS incidence, or pharmacodynamic markers (cytokine induction, T-cell margination,
   lymphocyte redistribution).
3. Ex vivo cytotoxicity EC50 against patient-derived SCLC at physiological E:T.
4. The same quantities for **any** half-life-extended CD3 bispecific (AMG 562, AMG 757 preclinical,
   other HLE BiTEs), which would at least bound the class.

**Deliver:** `EC50_EVIDENCE.md` — every candidate source, what quantity it actually measures, why it
does or does not map to an in vivo EC50, and a defensible range with the reasoning shown. If the
honest answer is "no in vivo anchor exists", **say that plainly** and state what the project is
therefore allowed to claim. That is a real and useful result, not a failure.

**Do not** pick a value because it makes the hypothesis hold. If the evidence points above 3 nM,
the hypothesis is dead and the correct move is to write that down.

**Environment note:** nature.com, thelancet.com, wiley and ashpublications return 403 here. Europe
PMC REST works: `https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=TERMS&format=json&resultType=core&pageSize=10`,
full text at `https://pmc.ncbi.nlm.nih.gov/articles/PMC<id>/`. Note the FDA label and EMA assessment
report are often the best source for occupancy and exposure–response and are usually fetchable.

---

## A2 — Retrieve SCLC tumour-nest geometry (unblocks Lead 2)

**Why.** `docs/recon/LUNG_CANCER_RECON.md` identifies SCLC tumour-nest / CD3+ stroma compartmentalisation as a
better-evidenced analogue of the follicular geometry question. It is currently marked **not
runnable** because nest dimensions could not be retrieved. The follicular project's own history says
guessed geometry is exactly how you produce a model artefact, so this stays blocked until measured.

**Needed, with sources:** tumour-nest diameter or area distribution; CD3+/CD8+ density inside nests
versus in stroma; distance from nest boundary to nearest T cell; thickness of the intervening
fibroblast/matrix layer. Treatment-naive limited-stage tissue is preferred.

**Known starting points:** Cell Discovery 2024 (PMC11375181, 132 TMA cores from 44 treatment-naive
LS-SCLC, three spatial compartments) — 403 via nature.com, try PMC. Cell Rep Med 2026 on combined
SCLC, and the spatial transcriptomics papers cited in the recon.

**Deliver:** `SCLC_GEOMETRY.md` with a calibration table in the same style as `README.md`'s — value,
units, source, confidence. Mark anything from a conference abstract as **low confidence and not for
calibration**, following the precedent already set for the CD8 core/interfollicular percentages.

If the numbers are not obtainable, say so and Lead 2 stays blocked. Do not interpolate.

---

## A3 — Independently verify the reference-model reading

**Why.** A large part of the current position rests on one reading of Obertopp/Basanta
(bioRxiv 2025.11.17.688873 / PMC12667981): that their TFI arms are a recurring duty cycle and that
all arms rest for days 29–42. That reading drives the conclusion that L1's schedule encoding was
wrong. It should not go unchecked by a second reader.

**Confirm or refute, quoting the paper directly:**

1. Are TFI arms recurring cycles, or a single interruption? Quote the methods.
2. **How many days are ON between OFF-intervals?** `exp_replicate.py` currently sweeps this because
   it could not be found. If it is stated anywhere — methods, figure captions, supplementary, code
   or data availability — that removes a free parameter and R4 becomes unnecessary.
3. Do all arms including CONT rest days 29–42, or only the TFI arms?
4. Confirm E:T 1:4, no T-cell influx, and the 5 h division refractory.
5. Their exhaustion: PD-1 threshold value, accrual rate, and the ~40%-death-in-10-days figure.
6. Is code or data published? If so, the duty cycle can be read directly rather than inferred.

**Deliver:** `REFERENCE_AUDIT.md` — each point confirmed, refuted, or not determinable, with quotes.
**If any point is refuted, say so prominently.** Several downstream claims depend on this and it is
better to find an error here than after the runs.

---

## A4 — Build the antigen-heterogeneity module (Phase 4 groundwork)

**Why.** The dominant clinical escape route for tarlatamab is selection of a pre-existing DLL3-low
population — the reference case showed the resistant clone "did not evolve from the clone that
predominated at treatment initiation… it may have been present before therapy and become dominant
due to the selective pressure of treatment" (PMC13082099). A model with uniform antigen cannot
represent that at all.

**Deliver:** `antigen.py` — a subclass or mixin over `Lymphoid` adding a heritable per-cell antigen
level `A ∈ [0,1]`, where the kill hazard scales with antigen as well as with drug occupancy and
`(1 − E)`. Daughter cells inherit the parent's level, with optional switching to represent
lineage plasticity as distinct from selection.

Requirements:
- **Do not modify `lymphoid.py`.** Subclass or compose. The existing model stays runnable as-is.
- Seeding must support an explicit initial antigen distribution, including a small pre-existing
  antigen-low subclone at a settable frequency — that is the whole point.
- Report, per run: antigen distribution over time, the antigen-low fraction, and whether change is
  attributable to selection (subclone outgrowth) or switching. These are separable in-model and
  should be reported separately, because distinguishing them is exactly the open question the
  clinical literature says is unresolved.
- Write a pre-registration in the module docstring, **before** any run, in the style of
  `exp_replicate.py`. Do not run experiments yet — Phase 4 is gated behind the calibration.
- Anchor DLL3 heterogeneity to what is published: positivity is defined as expression on ≥25% of
  tumour cells; 85–94% of SCLC is positive on that rule but only ~46.9% is high by H-score.

---

## A5 — Build the analysis harness for experiments R and T

**Why.** `analyze_L1_full.py` exists because the project's rule is that every pre-registered
endpoint gets reported regardless of outcome. R and T need the same, and it should be written
**before** results land so it cannot be shaped by them.

**Deliver:** `analyze_expR.py` and `analyze_expT.py`, reading `results/expR/*.json` and
`results/expT/*.json`.

`analyze_expR.py` must report R1–R4 explicitly by name with HELD/FAILED, plus:
- median burden at days 16/28/42 per arm, with the duty-cycle sweep laid out as a grid;
- exhaustion trajectories, and specifically whether CONT clears 0.50 by day 16 (that is R2);
- **the four-factor ablation table** — the single most important output, since it is what L1 never
  produced. One row per toggled factor, showing movement in the CONT-minus-best-TFI contrast;
- paired-by-seed comparisons with Mann–Whitney and Cliff's delta, matching the existing convention.

`analyze_expT.py` must report T1–T5 by name, including the PK-only occupancy table and the EC50
boundary, and must state the EC50 conditionality in the headline rather than in a footnote.

Both must tolerate partial results and say how many runs are missing — runs are written one file at
a time and a batch may be incomplete.

---

## Standing constraints — these apply to every task

- **No parameter fishing.** Calibrate to external measured curves, never to a schedule ranking.
- Report every pre-registered endpoint regardless of outcome. No reinterpreting hypotheses after
  seeing results.
- Preserve all code, seeds, raw outputs, pre-registrations, failed hypotheses and deviations.
- Record deviations rather than hiding them, in the style already used in `README.md`.
- No cure claim, no clinical recommendation, no proposed change to tarlatamab dosing.
- Final classification is exactly one of: A hypothesis falsified, B model artefact, C known result,
  D robust computational prediction, E independently supported prediction. Use no stronger language
  than the evidence permits.
- **A clean negative is a success.** A1 concluding that no in vivo EC50 anchor exists, or A2
  concluding the geometry is unobtainable, are good outcomes honestly delivered — not failures.
