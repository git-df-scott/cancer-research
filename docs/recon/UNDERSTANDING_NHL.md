# Understanding non-Hodgkin lymphoma before modelling it

Written 2026-09-07. This document exists because the previous session's experiment (`exp8_follicular.py`)
bolted a "lymphoma-shaped" geometry onto a generic solid-tumour model without first establishing what
lymphoma actually is. Four of its five pre-registered predictions failed, and in hindsight the more
serious problem is that the abstraction was wrong regardless of the result.

Sources are given inline. Claims marked **[unsourced]** are my background knowledge and have not been
checked in this session; treat them as provisional. A parallel 12-domain research sweep with adversarial
fact-checking is running and will be folded in.

---

## 1. "Curing non-Hodgkin lymphoma" is not a coherent target

NHL is not a disease. It is a family of roughly 60-plus entities in the WHO 5th edition (2022) and the
competing International Consensus Classification (2022), spanning B-cell, T-cell and NK-cell origins,
with wildly different natural histories. Some are already cured routinely; some are never cured; some
kill in weeks.

Any honest project has to name a specific entity and a specific failure mode. The two that matter most
by burden are diffuse large B-cell lymphoma (DLBCL, the most common, aggressive, curable in most
patients) and follicular lymphoma (FL, second most common, indolent, essentially never cured).

| Entity | Behaviour | Cure status |
|---|---|---|
| Burkitt | Extremely proliferative, Ki-67 approaching 100% | Highly curable with short intensive chemo |
| DLBCL | Aggressive | ~60-65% cured with frontline immunochemotherapy **[unsourced]** |
| Follicular | Indolent, median survival now ~20 years **[unsourced]** | Not cured; remission-relapse cycles |
| Mantle cell | Mixed | Generally not cured; TP53-mutant does badly |
| Peripheral T-cell | Aggressive | Outcomes lag far behind B-cell disease |

The paradox worth sitting with: **the fastest-growing lymphoma is the most curable.** Burkitt's
near-100% proliferation index is exactly what makes cycle-active chemotherapy able to sterilise it.
Indolence is not benign; it is what makes a lymphoma incurable.

---

## 2. The founder lesion is nearly ubiquitous and by itself harmless

Follicular lymphoma's signature is t(14;18), which places *BCL2* under the immunoglobulin heavy-chain
enhancer and blocks apoptosis.

The crucial fact is that this translocation is **not** the disease. t(14;18)-positive B cells are found
in roughly 10-50% of healthy people, with frequency rising with age (Oishi 2023, PMID 37380471).
What separates the healthy carrier from the future patient is the *second* class of hit: Schroers-Martin
et al. (2023, PMID 36939219) found *CREBBP* variants in 29% of pre-diagnostic FL cases versus 0% of
healthy BCL2-rearrangement carriers, years before clinical diagnosis.

FL's driver landscape is dominated by chromatin modifiers: *KMT2D*, *CREBBP*, *EP300*, gain-of-function
*EZH2*, *TNFRSF14*. It is best understood as an epigenetic disease of the germinal centre, not a
signalling-pathway disease. *CREBBP* loss additionally reduces MHC class II expression, which is an
immune-evasion event, not merely a growth event.

**Consequence for modelling:** a model whose resistance trait is a scalar "drug resistance" that arises
by mutation under drug pressure is telling the wrong story. The relevant heritable variation here
pre-exists treatment and is largely about immune visibility and differentiation state.

---

## 3. The germinal centre is the structure being corrupted

Most B-cell lymphomas arise from, or are arrested in, the germinal centre (GC) reaction. Understanding
the normal process is not optional background; it is the mechanism.

The GC is a transient structure in which B cells deliberately mutate their own immunoglobulin genes via
AID (activation-induced cytidine deaminase), then compete for survival signals. It has two zones:
the **dark zone**, where centroblasts proliferate and undergo somatic hypermutation, and the
**light zone**, where centrocytes are selected by T follicular helper cells and follicular dendritic
cells. Cells cycle between them.

Two features matter enormously and both are absent from the existing model:

1. **The GC is already a high-turnover mutate-select-die machine.** Enormous proliferation is matched by
   enormous apoptosis of non-selected cells. Lymphoma does not introduce this dynamic; it hijacks it.
2. **The light zone is physiologically hypoxic** (Pruitt & Abbott 2024, PMID 39559353; and PMID 39627218
   describing the LZ as "glycolysis-dominant, cell cycle arrest-inducing, hypoxic"). This hypoxia is a
   *regulated signalling environment in a small, well-vascularised structure*, not the diffusion-limited
   necrosis of an overgrown spheroid.

AID activity is why the GC is the origin of most B-cell lymphomas: a cell that is deliberately running a
mutator on its own genome will sometimes mis-target it.

**Consequence for modelling:** the existing model's hypoxia is generated by oxygen diffusing in from a
distant box boundary and being consumed until the core starves. That produces a necrotic core and a
"quiescence jail" that protects buried cells from cycle-specific drug. That mechanism is a solid-tumour
mechanism. Lymph nodes are lymphoid organs threaded with vasculature and sinuses; a neoplastic follicle
is on the order of a few hundred micrometres. The jail almost certainly does not exist in the form the
model builds it, and the entire quiescence-protection result does not transfer.

---

## 4. Why follicular lymphoma relapses: the reservoir, not the bulk

This is the central question and it has a well-supported answer that invalidates the framing I was using.

Relapse in FL frequently does **not** arise by linear descent from the treated tumour. It arises
divergently from a **common progenitor cell (CPC)**: a reservoir of t(14;18)-positive, lightly mutated,
memory-like B cells that carry the founder lesion and the earliest drivers but not the full driver set of
the clinical tumour. Diagnosis and relapse samples share an ancestor rather than one descending from the
other (Michaeli et al. 2022, PMID 36439408; Bai et al. 2024, PMID 39191762 tracking 94 longitudinal
biopsies from 44 patients and finding *CREBBP* and *KMT2D* mutations present at initiation).

**Consequence for modelling, and it is severe.** The existing project's entire narrative is: a resistant
subclone sits inside a tumour mass, competes for space with sensitive neighbours, and its expansion is
paced by how much drug is delivered nearby. That is a story about *evolution under drug pressure inside
the treated lesion*. For follicular lymphoma it appears to be the wrong story. The thing that causes
relapse is not inside the lesion competing for space. It is a separate, quiescent, differently-located,
drug-indifferent reservoir. You can win the competition inside the mass completely and still lose.

This is worth stating plainly because it is the single biggest thing the previous session got wrong.

---

## 5. The effective modern drugs are cells, not molecules

Lymphoma therapy has moved decisively to immune mechanisms, and this breaks the model's drug abstraction.

- **Rituximab / obinutuzumab** (anti-CD20) kill largely through *effector* mechanisms: antibody-dependent
  cellular cytotoxicity, complement, phagocytosis. They need a functioning immune system present.
- **Lenalidomide** works substantially through immune modulation rather than direct tumour killing.
- **Bispecific T-cell engagers** (mosunetuzumab, glofitamab, epcoritamab; CD20xCD3) physically bridge a
  T cell to a lymphoma cell.
- **CAR-T** is a living, dividing, exhausting drug.

The existing model represents a drug as a uniform scalar field that kills any cell attempting division.
For every agent listed above, that abstraction is wrong in a specific and consequential way: the killer
is a **discrete, motile, depletable, exhaustible agent that must physically find its target**.

### Quantitative anchors for a T-cell-based model

| Quantity | Value | Source |
|---|---|---|
| T-cell speed in lymph node cortex | ~11 µm/min, peaks >25 µm/min | Miller/Cahalan, Science 2002 |
| T-cell motility vs B cells | 5-6x higher motility coefficient | same |
| CTL killing capacity in vivo | 2-16 target cells per CTL per day | Halle et al., Immunity 2016, PMID 26872694 |
| Killing mode | Motile "kinapses", not static synapses | same |
| Cooperativity | Death probability rises when >2 CTLs contact a target | same |

That last row is the important one. **T-cell killing is cooperative and therefore genuinely spatial.**
A non-spatial model cannot represent "two T cells happened to converge on the same target."

---

## 6. Continuous exposure exhausts the effector: the real reason to pause treatment

This is the finding that reframes the whole project.

Philipp et al., *Blood* 2022 (PMID 35878001, PMCID PMC10652962) showed that continuous exposure to a
CD19xCD3 bispecific progressively destroys T-cell function, and that **treatment-free intervals reverse
it**:

- Continuous AMG 562 exposure: mean specific lysis fell from **88.4% at day 7 to 8.6% at day 28**
  (n=6, p=0.0003).
- With a treatment-free interval, day-14 specific lysis was **93.4% versus 34.9%** for continuous
  (n=6, p<0.0001).
- Exhaustion genes (*NR4A3*, *IRF4*, *PDCD1*, *LAG3*) were downregulated after the interval.
- In patients, 28-day continuous blinatumomab infusion likewise led to declining T-cell function.

The same principle holds for CAR-T: Weber et al., *Science* 2021 (PMID 33795428) showed transient rest,
via dasatinib or enforced CAR downregulation, restores function through epigenetic remodelling, and that
exhaustion "is not an epigenetically fixed state."

### The reframing

The existing project justified treatment holidays by **evolutionary** logic inherited from adaptive
therapy: keep sensitive cells alive so they suppress the resistant clone. In lymphoma treated with
T-cell engagers there is a completely different and better-evidenced reason to pause: **let the effector
cells recover.**

Both say "do not dose continuously," but they are different mechanisms, and critically they point in
**opposite directions on tumour burden**:

- Adaptive therapy deliberately *maintains* a large sensitive tumour population as a competitive
  suppressor.
- Antigen load is what drives T-cell exhaustion, so maintaining a large tumour population is exactly
  what would exhaust the effectors fastest.

If that tension is real, then adaptive therapy's core manoeuvre is actively harmful in T-cell-redirected
lymphoma therapy. That is a sharp, falsifiable claim and it is checkable in a model.

---

## 7. What has already been modelled, and the gap that is left

I searched specifically to avoid repeating existing work.

- **Germinal centre agent-based models are a mature field** (Meyer-Hermann and colleagues; e.g. PMIDs
  35911680, 36993964, 37818361, 38274834, 40108454). They model dark/light zones, cyclic re-entry, Tfh
  selection, somatic hypermutation, and are validated against intravital imaging. **None found extended
  to lymphoma or to therapy.**
- **Quantitative systems pharmacology models of T-cell engagers exist and are good** (Betts et al. 2019,
  PMID 31119428, trimer formation; Susilo et al. 2025, PMID 40284495, PBPK with a bell-shaped
  dose-synapse relationship). These are **non-spatial**. They cannot represent cooperativity, local
  effector depletion, or penetration into a dense aggregate.
- **Agent-based T-cell/tumour models exist for solid tumours** (Bergman et al. 2024, PMID 38515743,
  CD8 T cells in bladder cancer).
- **Closest of all, and directly on point:** Obertopp, Froid, Pilon-Thomas & Basanta, bioRxiv 2025
  (doi 10.1101/2025.11.17.688873, PMC12667981), "Rethinking the seven-day treatment-free interval in
  T-cell engager therapy using agent-based modeling." A 2D lattice ABM, tumour cells and T cells with
  inactive/active/exhausted states, PD-1 accumulation driving exhaustion, AMG 562 pharmacokinetics.
  It concludes the clinical 7-day treatment-free interval is **not** optimal and that shorter intervals
  of 2-3 days do better.

That last paper is essentially the experiment I was converging on, done first, by the Moffitt group that
invented adaptive therapy. Good. It means the well-mixed version of this question is answered and must
not be repeated.

**What it explicitly does not do**, from its own stated limitations and my reading:

1. It is **B-cell acute lymphoblastic leukaemia**, calibrated to an *in vitro petri-dish co-culture*.
2. It has **no tissue architecture**: fixed grid, 50% occupancy, **random rather than clustered seeding**
   (their words), T cells interspersed among tumour cells from t=0.
3. Tumour cells are **non-motile and pre-irradiated** (an in-vitro artefact).
4. **No hook effect** and **no cooperativity** (fixed 0.03 kill probability per contact, single T cell).

### The gap

Leukaemia is well-mixed; a T cell and a blast are already neighbours. **Lymphoma is not.** A follicular
lymphoma cell sits inside a dense follicle, and a T cell has to *get in*. The global effector-to-target
ratio is then close to meaningless: what matters is the local ratio at the invasion front, and the front
only advances as fast as killing opens space to move into.

That is precisely the dynamic this codebase already characterised for chemotherapy, where killing a
neighbour opened space for a resistant clone to expand into. The machinery transfers; the biology
inverts. Here, killing at the rim opens space for the **killer** to advance.

And it collides with exhaustion in a non-obvious way: a T cell that penetrates deep into a follicle is
surrounded on all sides by antigen, which is the fastest possible route to exhaustion. So there may be an
optimal penetration depth, and therefore an optimal schedule that differs from the well-mixed answer.

**Candidate contribution:** the optimal treatment-free interval is a function of tumour *architecture*.
Dispersed disease (leukaemia) and densely clustered disease (follicular lymphoma) should need different
schedules, for a reason that is mechanical rather than pharmacological.

There is a real-world test sitting right there. Blinatumomab in B-ALL is given by continuous infusion in
28-day cycles with a 7-day break. Bispecifics in follicular lymphoma are given on very different
schedules with long inter-dose gaps. If that divergence is mechanistically necessary rather than a
historical accident of drug development, a spatial model should be able to show why.

---

## 8. What this means for the existing code

| Component | Verdict |
|---|---|
| 2D lattice, one cell per site, Moore neighbourhood | **Keep.** Right abstraction for dense lymphoid tissue. |
| Paired design from identical pretreatment copies | **Keep.** This is the project's methodological strength. |
| Pre-registration of endpoints before running | **Keep.** Non-negotiable. |
| Front-fitting machinery (`frontfit.py`) | **Keep, invert.** Now fits an inward invasion front. |
| Multi-follicle builder (`follicular.py`) | **Keep.** Written yesterday for the wrong reason, useful for the right one. |
| Oxygen diffusion from box boundary | **Discard.** Lymph node is not a spheroid in medium. |
| Hypoxic quiescence jail | **Discard.** GC hypoxia is regulated signalling, not diffusion-limited starvation. |
| Drug as uniform scalar field | **Discard.** The killer is a discrete motile agent. |
| "Kills cells attempting division" | **Discard.** T-cell killing is not cycle-specific. |
| Resistance as an evolving scalar trait | **Discard for FL.** Relapse comes from a pre-existing reservoir, not in-lesion evolution. Antigen loss is the relevant escape mechanism. |

Roughly half the model survives. The half that survives is the half that took the longest to get right.

---

## 9. Verification pass (done by hand, not by the agent sweep)

The research sweep was stopped after its 12 domain surveys and before its fact-checking phase, so the
load-bearing claims were checked directly. Two corrections resulted.

| Claim | Status |
|---|---|
| Peri-follicular region is a barrier to immune infiltration; CD8 attack occurs at the rim, not the follicle core | **Verified.** Imaging mass cytometry, 13 paired FL/POD24 biopsies, 36-plex panel, J Hematol Oncol 2022, PMC9396877. Direct quotes confirmed. |
| Precise spatial figures: CD8 2.0% at core vs 28.3% interfollicular; TOX+PD-1+ 73.6% to 7.2%; 156 cores, 4.2M cells | **Downgraded.** Not found in the peer-reviewed literature; traces to a conference abstract. Treated as low confidence and NOT used for calibration. |
| Continuous engager exposure exhausts T cells; treatment-free intervals reverse it | **Verified.** Philipp et al., Blood 2022, PMID 35878001, PMCID PMC10652962. Specific lysis 88.4% (d7) to 8.6% (d28) continuous; 93.4% vs 34.9% at d14 with an interval. |
| CTL kills 2-16 targets/day in vivo, stays motile, cooperative above 2 contacts | **Verified.** Halle et al., Immunity 2016, PMID 26872694. |
| T-cell speed ~11 um/min in lymph node cortex | **Verified.** Miller, Wei, Parker & Cahalan, Science 2002. |
| Lymphoma cell loss factor ~93% (doubling 29 d, labelling index 29%, growth fraction 90%) | **Unverified.** Traces to 1970s kinetics literature I could not retrieve. Direction is corroborated by apoptotic-index data (0.59% low-grade vs 1.96% high-grade) but the exact figure is not used quantitatively. |
| Obertopp/Basanta agent-based model of engager treatment-free intervals | **Verified by direct fetch.** bioRxiv 2025, doi 10.1101/2025.11.17.688873, PMC12667981. Note it is a preprint, not peer-reviewed. |

An important consequence of the second row. The radial exhaustion gradient inside follicles is measured
in **untreated** biopsies, so it reflects baseline Tfh/Treg biology, not engager-induced killing. This
model therefore must **not** claim to retrodict that gradient. What it can address is what happens under
therapy given the barrier, which is a different and narrower claim.

---

## 10. The pilot measurement that defines the experiment

Same model, same cell numbers, same parameters, 21 days of continuous engager, only geometry differs.

| Architecture | Mean engaged fraction | Mean exhaustion at day 21 | Contact-minutes per T cell |
|---|---|---|---|
| One dense follicle | 0.085 | 0.032 | 2,118 |
| Dispersed, T cells interspersed | 0.610 | 0.139 | 10,679 |

A T cell in a follicle spends about a twelfth of its time in contact with a target, because only the
cells at the invasion front are touching anything. A T cell in a dispersed tumour is almost always in
contact. Since exhaustion accrues with dwell time under antigen, the same drug degrades the effector
pool roughly seven times faster in dispersed disease than in follicular disease, from geometry alone.

If treatment-free intervals exist to reverse exhaustion, and exhaustion barely accumulates in follicular
architecture, then intervals should buy much less in lymphoma than in leukaemia, and the cost of pausing
(letting the tumour regrow) should dominate sooner. That is experiment L1, pre-registered in
`exp_schedule.py`, and it is a direct test of whether the Obertopp/Basanta leukaemia result transfers.
