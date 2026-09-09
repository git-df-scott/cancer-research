# Codex — Groundwork review

Author: Codex.

Reviewed September 8, 2026, at commit `5d439e871693c3e5c3ab9f2c32ad30af5443ad9b` of `claude/lung-cancer-research-vl7oq3`. Scope: five executable modules, handoff claims, primary assay methods, reference source, and small independent diagnostics. No full calibration or tumour experiment was run; no repository code was changed or pushed.

**Verdict:** useful scaffolding and correct schedule arithmetic, but the calibration and downstream experiment are not yet a valid test of the biological hypothesis. The first task for Fable should be repairing the measurement and parameter handoff, not running the 32-cell fit.

Publication note: added to the branch after it advanced to `93f2363cb3be7d60458e31fb4ec03da2501a7c2d`. This review concerns the pinned `5d439e8` snapshot; subsequent calibration work has not been re-audited here. Statements about unchanged code or unrun experiments describe the review pass, not this later documentation commit.

## Confirmed

- Executing the schedule generators gives exposure fractions 2/3, 1 and 2/3 for reference CONT, L1 CONT and L1 TFI14, respectively. The reference paper's Figure 3/4 captions support a 28-day treatment window followed by a common 14-day rest. This establishes an exposure mismatch, not a causal attribution of the full ranking failure.
- All ten shipped PK self-tests pass. The quoted occupancy drops at the default volume ratio reproduce. These checks establish implementation consistency and matching of supplied quantities; they are not independent validation against measured concentration-time data.
- The reference's public code does increment PD-1 upon killing. That supports a model-structure difference, not proof that kill-driven exhaustion is the required biological repair.

## 1. Blocking: the calibration probe corrupts the experiment it measures

`calibrate_exhaustion.py:113–138` runs the exhausted and naive probes consecutively on the same live model. It saves only the exhaustion array and kill total, then restores that array after T cells have moved or died. Time, cell locations, target count, history and random-generator state remain advanced. Exhaustion is therefore restored to lattice positions, not the original cells.

Direct diagnostic using L=20, seed=7, dt=5, and all 40 initial T cells assigned E=0.6:

| Quantity | Before one probe | After |
|---|---:|---:|
| Simulation time, minutes | 0 | 480 |
| T-cell count | 40 | 38 |
| Mean E on living T cells | 0.6000 | 0.03158 |
| Non-T-cell sites carrying E | 0 | 38 |

This is a clean implementation counterexample, not an estimate of biological recovery. Day-7 probing alters the cells used for later continuous-arm endpoints. The TFI arm is not probed at day 7, introducing another asymmetry.

**Required repair:** use independent assay copies with explicit cell-state transfer into fresh assay conditions. Confirm the parent culture and RNG are unchanged by measurement. Never restore a position-indexed state array onto a moved population.

## 2. Blocking: the fitted quantity is not Philipp's measured quantity

The code computes `100 × kills_with_exhaustion / kills_after_zeroing_E`, using two four-hour replenished-target probes. Philipp's relevant target-lysis readout is a 72-hour assay at E:T=1:1; the four-hour experiment measures intracellular granzyme B. Chronic stimulation uses E:T=1:4 with specified replenishment/reculture, not random target replacement every simulation step. [Philipp et al., Blood 2022, Methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC10652962/), DOI [10.1182/blood.2022015956](https://doi.org/10.1182/blood.2022015956).

The harness's maintained target number also does not guarantee continuous cell contact. Its supposedly fixed effector population retains the default T-cell death rate. These can be useful modelling choices, but they are not a faithful assay replica without an observation model and sensitivity checks.

**Consequences:** the reported 87.1/62.5/96.4 values cannot establish “decay wrong, recovery right.” Matching one post-break point while the pre-break state and observation mapping are wrong does not identify `recover_tau`. Changing accrual changes the exhaustion present when recovery starts. Likewise, a ratio of measured lysis percentages does not directly equal remaining per-cell killing function.

**Required repair:** implement the target-count/control-normalized readout or justify a calibrated observation model; profile recovery jointly with accrual and observation uncertainty. Treat kill-driven versus dwell-driven effects as competing structural assumptions. A failure of one configuration does not uniquely establish which parameter must change.

## 3. Blocking: the fitted parameters and decision gates are not connected to R/T

`exp_replicate.py:126–129` and `exp_tarlatamab.py:142` construct `Lymphoid` without loading the calibration output or supplying fitted exhaustion parameters. Running the four commands in the handoff will therefore use the shipped defaults downstream unless code is manually changed.

`calibrate_exhaustion.py` computes SSE but defines no numerical acceptance threshold or accepted-fit artifact. Experiment R saves trajectories but does not enforce R1, and experiment T does not check for an accepted calibration or passing R verdict. R's docstring still contains the older decision tree, inconsistent with the new handoff.

**Required repair:** an immutable accepted-calibration artifact with parameter values, target-wise errors, uncertainty, acceptance rule, code/input hashes and seeds; explicit loading by R/T; a persisted R decision; refusal of dependent runs without required artifacts. Never accept the lowest grid loss merely because every other cell fits worse.

## 4. Blocking for T2–T4: exposure is not driving the claimed mechanism

The PK adapter is passed into the unchanged legacy model. Tonic exhaustion still depends on `drug > 0`, not its magnitude; recovery requires `drug == 0`. Thus every positive PK tail disables recovery. Also, `dwell_Tmin` is the final adjacency count of each day multiplied by 1440, summed across days (`exp_tarlatamab.py:149–152`). It neither integrates contact at each simulation step nor measures drug-mediated functional engagement.

**Required repair:** define and validate concentration-dependent activation and recovery, distinguishing adjacency from functional synapses. Accumulate the intended exposure/contact quantity at each timestep. Check zero, tiny positive and saturating exposure, no target, and withdrawal. Keep geometric contact as a separate diagnostic.

## 5. The comparators are not matched on the tested horizon

I integrated the analytic PK response over days 0–84 at the default 11.2-day half-life:

| Schedule | Listed total dose, mg | Dose strictly before day 84, mg | AUC(0–84), mg·day/L |
|---|---:|---:|---:|
| Q2W | 71 | 61 | 87.278 |
| “Matched” infusion | 71 | 71 | 95.522 |
| Weekly half-dose | 61 | 56 | 77.618 |

The Q2W history includes a 10-mg bolus exactly at day 84 that cannot affect the preceding simulated trajectory. Its listed total is nevertheless spread across the whole infusion comparator. More generally, equal eventual dose does not imply equal finite-window AUC. Weekly dosing also differs in loading and total dose.

**Required repair:** fix the matching definition and integration window before claiming timing effects; account explicitly for loading and terminal tails. Report actual dose and AUC. A different answer could otherwise be a dose/exposure effect.

## 6. Counterexample: T1 can fail at 1 nM without changing the half-life

The scripts use the default midpoint of the feasible peripheral-volume range; they do not sweep it despite the PK documentation saying they do. Retain the same CL, central volume, terminal half-life and 1 nM EC50, but set V2/V1 to 10% rather than 50% of its permitted maximum:

| Terminal half-life | Default V2/V1 | Default true steady-state drop | Alternative V2/V1 | Alternative drop |
|---|---:|---:|---:|---:|
| 5.8 days | 0.28933 | 22.30% | 0.05787 | **28.29% — fails T1** |
| 11.2 days | 1.02422 | 17.26% | 0.20484 | **26.92% — fails T1** |

These alternative parameters satisfy the code's mathematical constraints; they are not asserted to be clinically plausible. They show that the supplied constraints do not establish T1 at 1 nM. The missing hypothesis is a sufficiently constrained PK distribution shape, in addition to potency and the functional-response mapping.

For positive peak/trough concentrations P and T under the assumed `C/(C+K)` response, a 25% relative occupancy-drop boundary is

`K* = 0.25*T*P / (0.75*P - T)`

when `T < 0.75*P`. If that condition fails, the drop never exceeds 25% at finite positive K. This is an algebraic property of the assumed curve, not a biological discovery. At the default parameters the true steady-state thresholds are approximately 1.17 and 1.63 nM. They move with PK shape. The shipped table samples days 56–70, which closely approximates these default steady-state drops but is not an exact steady-state calculation.

The justified current label is **conditional calculation under an incompletely constrained PK/response model**. “Classification D” is not earned merely by appending “conditional on EC50.”

## Further issues to fix before scaling

- **Ablation is confounded.** `L1+ET` also changes grid size, tumour count and spatial density. It is a configuration bundle, not an isolated E:T intervention. Comparing “best TFI” across different candidate schedule sets adds selection effects. Report the bundle honestly or redesign controlled contrasts; one-at-a-time changes cannot generally identify interactions.
- **Reference source is available.** The paper explicitly links [ninaobertopp/TFI](https://github.com/ninaobertopp/TFI). I inspected commit `2858dd0534cce7e0443a3e7dfea76bd14b7d2740`, file `ALL`. It includes schedule arrays and commented alternatives, kill-linked PD-1 increments, a 12-hour killing refractory condition, and a probabilistic exhaustion transition. The checked snapshot alone does not establish every published schedule. Resolve these differences before calling the rig faithful.
- **Resume is unsafe across revisions.** R/T filenames omit calibration/configuration/code hashes, so existing results can silently be reused after parameters change; JSON writes are not atomic and existence alone does not detect incomplete files. The calibration scan writes only at the end and retains averages, not all individual seed results. Its default is one seed per cell. The “a stall costs one run” claim does not hold for this calibration implementation.
- **Smoke evidence was deleted.** GROUNDWORK explicitly says validation runs were removed. I did not independently reproduce the large-grid floor result. Preserve diagnostic results with a smoke-test label; deleting them reduces auditability rather than preventing overclaiming.
- **Coordinate with Phase 4.** PR #1 already contains a different calibration and systemic model. These new modules still use the older model. Keep their results and intended questions distinct rather than treating them as one continuous validated instrument.

## Recommended handoff to Fable

1. Repair and verify the assay readout and non-mutating probes.
2. Establish an externally justified acceptance rule and jointly assess recovery/decay uncertainty.
3. Connect accepted parameters to R/T, enforce gates, and make result reuse configuration-safe.
4. Correct the concentration-response coupling, dwell integration and exposure matching.
5. Constrain or sweep PK shape together with potency; narrow the T1 claim accordingly.
6. Run a small preserved pilot; only then commit to the full fit/replication/ablation/T batch.

The most valuable next test is a faithful assay measurement that leaves its parent culture unchanged. Until it passes, more calibration compute can produce a more precise fit to the wrong quantity.

## Reproduction and limits

Independent checks are in [CODEX_GROUNDWORK_REVIEW_CHECKS.json](CODEX_GROUNDWORK_REVIEW_CHECKS.json); [CODEX_verify_groundwork.py](CODEX_verify_groundwork.py) reproduces the small diagnostics using NumPy/SciPy and the pinned checkout. All ten PK self-tests were run. No 32-cell scan or 360-job campaign was run. The source review and synthetic mutation test establish implementation problems; they do not establish a therapeutic mechanism or quantify patient effects. Reference paper inspected via Europe PMC full-text XML; source access date September 8, 2026.
