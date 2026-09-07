# Handoff prompt — continue the follicular lymphoma / T-cell-engager scheduling project

Continue an in-progress computational project. Clone and read before doing anything:

```
git clone -b lymphoma-tce-architecture https://github.com/git-df-scott/cancer-research.git
cd cancer-research
python -m venv .venv && .venv/bin/pip install numpy scipy matplotlib
```

Read in this order: `FINDINGS.md` (live status), `README.md`, `UNDERSTANDING_NHL.md` (the biology),
`PRIOR_ART.md` (novelty audit and claim matrix), `TRAFFICKING_PREREG.md`, then the docstrings of
`exp_schedule.py` and `exp_attacks.py`, which contain the pre-registrations.

## The question, unchanged

Does follicular lymphoma architecture genuinely invalidate or materially alter treatment-free-interval
conclusions derived from dispersed / leukaemia-like T-cell-engager models?

Candidate chain: follicular architecture → smaller fraction of T cells simultaneously contacting
malignant B cells → lower cumulative contact dwell time per effector → slower exhaustion → less
benefit from treatment-free recovery intervals → different optimal schedule.

## Where it stands, and the one thing that matters

Experiment L1 is complete: 150 runs, 10 seeds, 3 architectures × 5 schedules, all raw output in
`results/expL1.json`.

- **P1 and P2 HELD.** Follicular architecture gives a 4.7× lower engaged fraction (0.098 vs 0.463)
  and 4.8× lower exhaustion (0.024 vs 0.116) than dispersed. Links 1 and 2 of the chain work in-model.
  Neither is novel; geometry changing contact frequency is expected and has prior art.
- **P3, P4, P5 FAILED.** No treatment-free interval beat continuous dosing in any architecture.
  Continuous won everywhere and longer interruptions were monotonically worse.

**The experiment is uninformative, because the positive control failed.** The dispersed arm was
supposed to reproduce the reference result (Obertopp/Basanta bioRxiv 2025: short TFIs beat continuous
in well-mixed leukaemia-like disease) and produced the opposite ranking. A model that cannot
reproduce the effect where the effect is known cannot be used to argue the effect fails to transfer.

**Root cause, quantified.** Exhaustion reaches at most 12% functional loss at day 28. Philipp et al.
(Blood 2022, PMID 35878001) measured specific lysis falling 88.4% → 8.6% over 28 days of continuous
engager exposure, roughly 90% loss. So exhaustion never becomes the binding constraint: a 7-day
interval costs 25% of drug-time and buys back ~1-2% of function, and can never win in any geometry.
Most likely mechanical cause is T-cell influx diluting the exhaustion pool — the T-cell population
grows 200 → ~860 over a run, so most effectors at day 28 are recent arrivals.

## Your first task, and do not skip it

**Recalibrate exhaustion against Philipp's external curve, then verify the positive control passes.**

1. Make the model reproduce ~90% loss of killing function after 28 days of continuous engager
   exposure, and recovery on a treatment-free interval (their day-14 comparison: 93.4% with an
   interval vs 34.9% continuous). Candidate levers: exhaustion accrual rate, whether influx dilutes
   the pool, whether exhaustion should attach to the standing population rather than be averaged over
   arrivals. Calibrate to the measured curve, which is fixed and external. **Do not calibrate toward
   any schedule ranking.**
2. Then confirm the dispersed arm reproduces short-TFI-beats-continuous. **If it does not, stop and
   report that. The architecture comparison is meaningless until this passes.**
3. Only then re-run the architecture comparison.

## Then the attack battery, already written and pre-registered

`exp_attacks.py` is ready to run and contains, with decision rules fixed in advance:

- **Phase 3, exhaustion knockout.** If the architecture × schedule interaction survives with
  exhaustion disabled, the proposed causal chain is false; find what actually causes it.
- **Phase 4, dose-matched controls.** TFI arms deliver less drug. A continuous arm at level
  (28−k)/28 delivers matched exposure. If a TFI does not beat its own dose-matched twin, there is no
  timing effect, only dose.
- **Phase 5, trafficking knockout.** `lymphoid.py` has a `swap_prob` parameter letting T cells
  squeeze past malignant B cells instead of requiring a vacancy. **This is the biggest threat and it
  is already partly established:** at swap_prob 0.5 the follicle engaged fraction rises from 0.048 to
  0.576, essentially the dispersed value, with T cells reaching 93% of the way to the follicle centre
  *without any killing*. The 11 µm/min T-cell speed that calibrates the model was itself measured in
  densely cellular lymph-node cortex, so absolute volume exclusion is likely wrong. Pre-registered
  rule: if the interaction shrinks by more than half at swap_prob 0.5, classify as **model artefact**
  and do not present it as an FL prediction.
- **Exhaustion distribution.** Report exhaustion among engaged vs unengaged T cells and by radius. If
  the mean is low while the front is fully exhausted and continuously replaced from a fresh reservoir,
  that is a different mechanism and the claim must be restated, not preserved.

Then Phase 6 (parameter robustness — report the fraction of plausible space where the conclusion
survives, and find the boundary where ranking flips), Phase 7 (geometry continuum rather than binary
categories; find the dimensionless quantity that predicts any transition), Phase 9 (an independent
prediction written down *before* looking for data to test it).

## Hard constraints

- Do not restart the general literature audit. Phase 8 is done; the two closest works were read in
  full and neither duplicates the claim. See the claim matrix in `PRIOR_ART.md`.
- Do not search for a positive result by tuning parameters after seeing outcomes.
- Do not claim a new biological mechanism. Every individual link has prior art.
- Report every pre-registered endpoint regardless of outcome. Do not reinterpret hypotheses after
  seeing results.
- Preserve all code, seeds, raw outputs, pre-registrations, failed hypotheses and deviations.
- Final classification must be exactly one of: A hypothesis falsified, B model artefact, C known
  result, D robust computational prediction, E independently supported prediction. Use no stronger
  language than the evidence permits.
- Success means establishing one correct thing that was not previously established, **or** rigorously
  demonstrating why the idea fails. A clean negative is a success.

## Environment notes learned the hard way

- The local laptop has 6 cores and everything spawned inherits nice 5. **A previous L1 run deadlocked
  at 125/150 when system swap filled.** Cloud should be free of this, but keep per-result file writes
  (`rerun_missing.py` shows the pattern) so a stall costs one run, not everything.
- Runtime: ~54 s per 42-simulated-day run at dt=5, L=120.
- dt=5 with sub-stepped motility was validated against dt=1: agreement within 13%, and the
  discretisation is common-mode across arms. Smoke run matched to 9112 vs 9117 kills.
- Elicit MCP tools are blocked on this account. nature.com, thelancet.com, wiley, ashpublications
  return 403. Europe PMC REST works well:
  `https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=TERMS&format=json&resultType=core&pageSize=10`
  Full text at `https://pmc.ncbi.nlm.nih.gov/articles/PMC<id>/` (not `www.ncbi.nlm.nih.gov/pmc/`).

## Key sources already verified by direct fetch

| Fact | Source |
|---|---|
| Peri-follicular region is a barrier to immune infiltration; CD8 attack at the rim not the core | J Hematol Oncol 2022, PMC9396877 (13 paired FL/POD24 biopsies, 36-plex IMC) |
| Continuous engager exposure exhausts T cells; intervals reverse it | Philipp, Blood 2022, PMID 35878001 |
| CTL kills 2-16 targets/day in vivo, stays motile, cooperative above 2 contacts | Halle, Immunity 2016, PMID 26872694 |
| T-cell speed ~11 µm/min in lymph node cortex | Miller/Cahalan, Science 2002 |
| Closest prior model (leukaemia, well-mixed, random seeding) | Obertopp/Basanta, bioRxiv 2025, doi 10.1101/2025.11.17.688873 |
| Architecture+exhaustion, but emergent density, rate-parameter exhaustion, no scheduling | PhysiCell trajectory paper, PMC13060115 |
| Exhaustion from cumulative antigen integral, emergent density, no scheduling | TIED-ABM, PMC13541251 |

Note: precise spatial percentages (CD8 2.0% core vs 28.3% interfollicular) come from a **conference
abstract**, are low confidence, and are deliberately not used for calibration.

No hype. No cure claims. No parameter fishing. Find out whether the effect is real.
