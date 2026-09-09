# Disease-level candidate screen: which cancer should a compute-only effort attack?

Same discipline that has worked all session at the hypothesis level, applied at the disease level:
fix the criteria first, then kill candidates fast.

## What hand we are actually playing

We have public data, arithmetic, literature access, and a demonstrated willingness to discard our
own work. We do not have a wet lab, clinical data access, patients, or money. The screen selects
for what that hand can play, not for what is most important in the abstract.

## Criteria, fixed before scoring

| | |
|---|---|
| `cure_exists` | Is anyone cured today? If yes, the mechanism is proven and the question becomes "why not everyone" — far more tractable than "is cure possible" |
| `dependency` | Is there a defined molecular dependency? Curable cancers overwhelmingly have one |
| `gap_size` | How many fail current therapy? |
| `public_data` | Downloadable quantitative data — DepMap, TCGA, SEER, published PK |
| `actionable` | Does a drug or trial exist that a finding could redirect? |
| `uncrowded` | Can we add what 500 better-resourced labs cannot? |
| `compute_only` | Can the gap be attacked by analysis rather than experiment? |

Scored 0–3. **Zero on `cure_exists`, `public_data` or `actionable` kills the candidate outright** —
those are necessary, not merely desirable. Reproducible in `cancer_screen.py`.

## Result

**Killed on `cure_exists`:** glioblastoma, pancreatic, SCLC. Nobody is cured, so there is no proven
mechanism to extend. Not pessimism — the criterion working.

**Killed on `uncrowded` (scored 0):** childhood ALL, DLBCL, AML, myeloma, metastatic NSCLC. Real
gaps, but we will not out-resource those fields.

**Low `gap_size`:** testicular, papillary thyroid, Hodgkin. Already cured; little left to win.

**Two survivors tied at 18: high-risk neuroblastoma and Ewing sarcoma.**

## Tiebreak to Ewing, and why

Ewing has the property that defines curable cancers — **one defining lesion**, EWSR1-FLI1 fusion in
~85%. And Ewing lines are among the most PARP-inhibitor-sensitive in all of cancer cell line
screening: strong, reproducible, textbook.

Then the trial ran.

> Phase II olaparib in refractory Ewing sarcoma: **12 patients, 0 objective responses**, median time
> to progression 5.7 weeks. ([PMC4230717](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4230717/))

The strongest preclinical biomarker in oncology produced nothing. That gap has never been settled.
Published explanations: insufficient synthetic lethality; autophagy-mediated resistance;
**"failure to achieve in vitro levels of olaparib at the clinical dose"**; and the need for
combination.

**The third is a pharmacokinetics question that has been proposed and never computed.**

## The question we would answer

> Is the **free** plasma concentration of olaparib at the approved dose ever above the
> concentration that killed Ewing cells in vitro?

Likely culprit: protein binding. Olaparib is ~82% plasma protein bound; in vitro assays run in 10%
serum where binding is far lower. The free drug a tumour cell sees in a patient can sit several-fold
below the nominal in vitro IC50, and the correction is routinely skipped.

Why it is ours: public data on both sides (published olaparib popPK; Ewing IC50s in GDSC/CCLE and
primary papers), one day of work, no model requiring validation — the trap that consumed three
weeks — and the same analysis shape as the EC50-free occupancy work, where the method held up and
we caught our own errors twice.

## Why the answer matters either way

| If | Then |
|---|---|
| Free concentration never reached the in vitro IC50 | The biology was right and the *drug* was wrong. Combination strategies using *lower* PARP doses are the wrong direction; a more potent or better-distributed PARP inhibitor is |
| Free concentration did reach it | The preclinical model was wrong, and autophagy / insufficient synthetic lethality are where to look |

Nobody has drawn that line. The trial authors listed it as a possibility and moved on.

## Recorded caveats

- The trial is n=12. Small.
- The field has since moved to olaparib + temozolomide combinations.
- This does not cure Ewing. Best case it redirects a drug-development strategy in a cancer killing
  ~30% of localised and ~70% of metastatic patients.
- Neuroblastoma HR tied on score and is not eliminated — it is the fallback if the Ewing question
  turns out already settled somewhere we have not looked.
