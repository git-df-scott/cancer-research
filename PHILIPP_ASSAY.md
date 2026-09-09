# Phase B — What actually generated Philipp's calibration numbers

Source: Philipp et al., *Blood* 2022, PMID 35878001, full methods read at
[PMC10652962](https://pmc.ncbi.nlm.nih.gov/articles/PMC10652962/). This supersedes every
summary-based description of the assay used earlier in this project, including the docstring of
`calibrate_exhaustion.py`.

## The single most important finding

**The molecule is AMG 562, a half-life-extended CD19xCD3 BiTE — not blinatumomab.**

This matters for the lung question. The exhaustion-versus-treatment-free-interval evidence this
project calibrates against was generated with a *half-life-extended* engager, the same class as
tarlatamab. It is a closer analogue than assumed, and it also means the "TFI literature descends
from a 2-hour-half-life molecule" framing in `LUNG_CANCER_RECON.md` is too broad: it is true of the
Obertopp/Basanta ABM's blinatumomab lineage, but not of Philipp, whose rest intervals were imposed
by physically removing drug at reculture rather than by letting concentration decay.

## Two different cultures, two different E:T ratios

Earlier work in this repository conflated them. They are separate experiments.

| | Chronic stimulation (where exhaustion accrues) | Cytotoxicity readout (where lysis is measured) |
|---|---|---|
| E:T ratio | **1:4** | **1:1** |
| Targets | irradiated OCI-Ly1 | fresh hCD19-Ba/F3 or OCI-Ly1 |
| Duration | 28 days, 4 cycles of 7 days | **72 hours** |
| Drug | AMG 562 present (continuous arm) | AMG 562 at 5 ng/mL |
| Replenishment | medium, targets and AMG 562 on day 3; T cells re-isolated and recultured with fresh targets and drug on day 7 | none — single 72 h incubation |
| Normalisation | n/a | vs a **control construct** (cBiTE) |

## The measured observable, stated exactly

Equation 1 of the paper:

```
% specific lysis = (1 - CD19+ target count WITH engager / CD19+ target count WITH control construct) x 100
```

Three consequences, each of which the previous calibration rig got wrong:

1. **The denominator is a no-engager control run with the same T cells** — not naive or
   unexhausted T cells. The old probe divided by "the same cells with E set to zero", which is a
   *naive-effector* reference. That is a different denominator and a different quantity.
2. **It is target depletion over 72 hours**, not an instantaneous kill rate. The old probe used a
   4-hour window and counted kill events.
3. **It is measured at E:T 1:1 against fresh targets**, on T cells harvested out of the 1:4 chronic
   culture. Exhaustion accrues in one condition and is read out in another.

## The 4-hour assay is a different endpoint entirely

The 4-hour incubation with GolgiStop/GolgiPlug followed by intracellular staining measures
**granzyme B content in T cells**, not target death. It must not be conflated with the 72 h
specific-lysis readout. Codex's finding 2 was correct on this point.

## Provenance of the four calibration values

| Value | Figure | Assay | n | Error |
|---|---|---|---|---|
| 88.4% (d7) vs 8.6% (d28), continuous | Fig 2E | AMG 562-mediated 72 h lysis | 6 | mean ± SEM, *P* = .0003 |
| 34.9% (continuous) vs 93.4% (TFI), d14 | Fig 3E | AMG 562-mediated 72 h lysis | 6 | mean ± SEM, *P* < .0001 |

All four are the same observable. `n = 6` donors, mean ± SEM.

## The treatment-free interval, as actually performed

Drug was absent during stimulation cycles 2 (days 7–14) and 4 (days 21–28). **Target cells were
still present** during those windows — it is a drug-free interval, not a rest from antigen. The
previous model rig imposed rest by setting `drug = 0` while targets remained, which is correct.

## Protocol → observable → required model observable

| Experimental protocol | Measured observable | Model observable required |
|---|---|---|
| 28 d chronic coculture, E:T 1:4, irradiated targets, drug ± per arm, replenished d3, recultured d7 | latent T-cell functional state | per-cell exhaustion `E` carried by cell identity, not by lattice site |
| harvest T cells → 72 h coculture, E:T 1:1, fresh targets, ± engager | % specific lysis vs control construct | `(1 - N_target(72h, drug on) / N_target(72h, drug off)) x 100`, same harvested cells in both arms |

## Marked unavailable rather than inferred

- Plate format and wash steps for the lysis assay.
- Medium composition (referenced to a supplemental table not retrieved).
- Flow-cytometry gating strategy.
- Absolute target seeding densities in either culture.
- Whether the day-3 replenishment adds targets to the existing well or replaces them.

The last item affects how the chronic culture is modelled and is currently **unavailable**. The
implementation adds targets back to the set point and records this as an assumption.
