# Phase 2: Philipp's population measurements, reconstructed from primary methods

Source: Philipp et al., *Blood* 2022, PMID 35878001, methods and results read at
[PMC10652962](https://pmc.ncbi.nlm.nih.gov/articles/PMC10652962/).

**This reconstruction retracts part of the structural counterexample recorded in
`docs/findings/STRUCTURAL_COUNTEREXAMPLE.md`.** See `docs/findings/COUNTEREXAMPLE_VERDICT.md`.

## The decisive definition

```
Fold change = CD2+ cell count day 3 / CD2+ cell count day 0
```

**This is a 3-day assay readout, not the chronic culture's population trajectory.** It measures the
proliferative *capacity* of T cells harvested from the chronic culture at a given timepoint, in the
same way the 72 h lysis assay measures their cytotoxic capacity. Both are functional endpoints
measured on harvested cells; neither reports the size of the chronic culture.

I previously read 4.1 vs 1.1 as "the real culture expands four-fold more in the rested arm over 28
days". That reading is wrong. The correct reading is "T cells harvested at chronic-culture day 14
expand 4.1-fold over 3 days if rested, versus 1.1-fold if continuously stimulated."

## Reconstructed table

| arm | chronic day | endpoint | value | assay | provenance |
|---|---|---|---|---|---|
| continuous | 7 | specific lysis % | 88.4 | 72 h, E:T 1:1 | Fig 2E, n=6, mean ± SEM |
| continuous | 14 | specific lysis % | 34.9 | 72 h, E:T 1:1 | Fig 3E, n=6 |
| continuous | 28 | specific lysis % | 8.6 | 72 h, E:T 1:1 | Fig 2E, n=6 |
| TFI | 14 | specific lysis % | 93.4 | 72 h, E:T 1:1 | Fig 3E, n=6 |
| TFI | 28 | specific lysis % | 58.7 | 72 h, E:T 1:1 | Fig 3E |
| continuous | 14 | CD2+ fold change | 1.1 | **3-day** proliferation assay | Fig 3, n=6, *P*=.002 |
| TFI | 14 | CD2+ fold change | 4.1 | **3-day** proliferation assay | Fig 3, n=6 |
| continuous | 28 | CD2+ fold change | 0.06 | **3-day** proliferation assay | Fig 3, *P*=.09 |
| TFI | 28 | CD2+ fold change | 2.8 | **3-day** proliferation assay | Fig 3 |

## Explicitly NOT available, and not inferred

- **Absolute viable T-cell counts in the chronic culture at any intermediate timepoint.** Not
  reported. So the chronic culture's population trajectory is *unmeasured*, not measured-and-
  contradicted.
- **Any division-history or proliferation marker.** No Ki-67, no CFSE, no CellTrace. CD2+ fold
  change is a net count ratio over 3 days and cannot separate division from survival.
- **An explicit statement that effectors were counted before the cytotoxicity assay.** The methods
  state `E:T = 1:1` but do not describe the counting step.
- Plate format, wash steps, medium composition, gating strategy.

## On whether equal numbers were plated

The methods specify **E:T = 1:1** for the cytotoxicity assay. A ratio cannot be set without
counting effectors, so equal-number plating is entailed by the stated design even though the
counting step is not described. Recorded as **strongly implied, not explicitly stated**.

This matters because it determines whether chronic-culture population size can enter the lysis
observable at all. If effectors are replated to a fixed ratio, it cannot.

The current model already standardises plating (`n_plate = 250`, sampled from the harvested
exhaustion distribution), which matches this reading of the protocol. That fix was made for a
different reason and turns out to be protocol-correct.
