# External validation candidates, assessed BEFORE computing any M1 prediction

Written while the clean reselection is still running and before any prediction is generated, so
the choice cannot be influenced by whether M1 happens to predict a candidate correctly.

## Admissibility criteria, fixed

A candidate is usable only if all hold:

1. same or comparable biological mechanism;
2. quantitatively extractable;
3. its source informed no current model default;
4. its source informed no fitted parameter;
5. its source did not drive M0-versus-M1 selection;
6. protocol reproducible in the model well enough to generate a prediction **without adding a
   mechanism**;
7. replicate/uncertainty information sufficient to set a defensible acceptance criterion.

## Candidates

| # | Source | Molecule / system | What it measures | Provenance | Usable? |
|---|---|---|---|---|---|
| 1 | Philipp, *Blood* 2022 (PMID 35878001) — d7/d14/d28 lysis | AMG 562, OCI-Ly1 | 72 h specific lysis | **Contaminated.** Fitted `exhaust_tonic`; `recover_tau` default cites it | **No** — it is the calibration set |
| 2 | Philipp, same paper — d14_tfi, d28_tfi | as above | 72 h specific lysis after rest | **Contaminated** via `recover_tau` (cites "TFI reinvigoration, Philipp 2022") | **No** |
| 3 | Philipp, same paper — CD2+ fold change | as above | 3-day expansion capacity | Same source; also M1 has no expansion observable | **No** |
| 4 | Zugmaier/Subklewe ASH 2018 abstract, "Treatment-Free Intervals Mitigate T-Cell Exhaustion…" | CD19xCD3 BiTE | specific lysis, long-term culture | Same group, same culture system, precursor of candidate 1. Reports d14 CONT 34±4.2, TFI 99±2.2 — the same experiment | **No** — same data lineage |
| 5 | Weber et al., *Science* 2021 | CAR-T, rest via dasatinib | function after rest | **Contaminated** — `recover_tau` default cites it explicitly | **No** |
| 6 | Halle et al., *Immunity* 2016 | in vivo CTL | killing capacity | **Contaminated** — `p_kill` provenance | **No** |
| 7 | Obertopp/Basanta PMC12667981 | ABM | simulated schedules | Not experimental data; also informed the kill-driven accrual question | **No** |
| 8 | **ASH 2023 abstract, "Accomplices in Cure: Blinatumomab + Dasatinib…"** | blinatumomab + dasatinib, long-term culture | specific lysis ~80%, PD-1+Tim-3+LAG-3+ 0.6% vs 31% at day 14 | **Clean on 3, 4, 5** — never used for any default, fit or model selection | **No — fails 6 and 7**, see below |

## Why candidate 8 fails, despite clean provenance

It is the strongest candidate and the only one that passes the provenance tests. It fails on two
other criteria, and the failures are not fixable by trying harder.

**Fails criterion 6 (no new mechanism).** The intervention is intermittent high-dose dasatinib
transiently switching off TCR signalling while the engager remains present. M1 has no TCR-signalling
term and no second drug. The only way to represent it is to assert that "dasatinib on" maps to
occupancy 0 — treating a signalling blockade as pharmacologically identical to engager withdrawal.
That mapping is an untested assumption of exactly the kind this project keeps being burned by, and
adopting it would make any pass uninterpretable: a success would confirm the assumed mapping, not
M1.

**Fails criterion 7 (uncertainty).** It is a conference abstract. No replicate structure, no error
bars, and the dasatinib duty cycle is not specified in the retrievable text — "intermittent
high-dose" without days on/off. A prediction cannot be computed without inventing the schedule, and
an acceptance interval cannot be set without dispersion.

This repository already treats conference-abstract numbers as low confidence and excludes them from
calibration (see `docs/recon/PRIOR_ART.md` on the CD8 core/interfollicular percentages). Admitting
one as a *validation* target would be a weaker standard than the one applied to calibration, which
is backwards.

## Conclusion

**NO PROVENANCE-CLEAN EXTERNAL VALIDATION IS AVAILABLE.**

Every quantitative observation of chronic engager stimulation with a functional readout that could
be located either comes from the paper that set M1's parameters, or cannot be predicted without
adding a mechanism and inventing a schedule.

This is stated as the result of the search, not as a reason to lower the bar. Classification A
remains unreachable. The correct next step is a prospective discriminating experiment, not a
weaker validation target.
