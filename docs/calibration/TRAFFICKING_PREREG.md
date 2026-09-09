# Phase 5 pre-registration: attacking the central spatial assumption
Written 2026-09-07 BEFORE L1 results were examined and BEFORE any alternative was run.

## The assumption under attack

`lymphoid.py` currently lets a T cell move only into an **empty** lattice site. In a packed follicle
there are no empty sites, so T cells cannot enter until killing creates a vacancy. That rule is what
produces rim-to-core invasion, and it is the single load-bearing structural assumption of the whole
project.

## Why I now think that assumption is probably WRONG as stated

This is the honest reading of the data I already gathered, and it goes against my hypothesis.

- T cells migrate at ~11 um/min **in the lymph node cortex**, which is itself densely cellular
  packed lymphoid tissue, not open space (Miller/Cahalan, Science 2002). The measurement that
  calibrates T-cell speed was made in exactly the kind of crowded tissue my rule forbids traversing.
- Tfh cells move within intact germinal centres, which are among the most densely packed tissues in
  the body. Germinal-centre agent-based models let lymphocytes move through them.
- Lymphocytes are highly deformable and migrate amoeboid-style, squeezing between cells. Absolute
  volume exclusion at the scale of one cell diameter is a lattice convenience, not biology.

So "T cells cannot pass a malignant B cell" is very likely a model artefact. If the L1 effect exists
only under that rule, it is a **model artefact (category B)** and must not be presented as an FL
prediction.

## The alternative formulation, committed to now

**Neighbour exchange (swapping) with a density-independent exchange probability.** A moving T cell
that targets a site occupied by a malignant B cell may swap places with it, with probability
`swap_prob` per attempt. Volume exclusion is preserved (one cell per site, no overlap), but the
tissue becomes traversable. `swap_prob = 0` recovers the current absolute-exclusion model exactly.

This is the standard way motile cells traverse crowded tissue in lattice models, and it is the
minimal change that removes the suspect assumption without introducing new free structure.

### Values to test, chosen from data, before seeing any outcome

| swap_prob | Justification |
|---|---|
| 0.0 | current model, absolute exclusion; the assumption under attack |
| 0.5 | T cells traverse packed tissue at about half the rate of free movement. Conservative given that the 11 um/min calibration was itself measured in dense cortex |
| 1.0 | full traversal; occupancy does not impede motility at all. The extreme case |

I expect **0.5 to be the most biologically defensible** of the three, and I am recording that
expectation now so I cannot later prefer whichever value keeps the result.

## Pre-registered Phase 5 predictions

- **T1** Under swap_prob = 0.5, the engaged fraction in the follicle arm rises substantially
  (at least doubles relative to swap_prob = 0), because T cells reach interior targets without
  waiting for vacancies.
- **T2** If the architecture x schedule interaction seen in L1 (if any) shrinks by more than half at
  swap_prob = 0.5, the effect is classified **model-dependent** and will be reported as category B,
  not as an FL prediction.
- **T3** If the interaction survives essentially unchanged at swap_prob = 0.5 AND at 1.0, the effect
  does not depend on the vacancy rule and survives this attack.

## Decision rule, fixed in advance

Phase 5 is run at swap_prob 0.5 and 1.0 regardless of what L1 shows. If L1 shows no interaction,
Phase 5 still runs at swap_prob = 0.5 as a check that the *absence* of an interaction is not itself
an artefact of the exclusion rule.
