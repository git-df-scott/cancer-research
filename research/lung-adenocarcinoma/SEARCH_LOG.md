# Audit search and access record

Date: 2026-09-07. Focused audit, not a registered systematic review. Scope was selected to test the proposed EGFR-mutant lung adenocarcinoma persister project for prior art, credible therapeutic evidence, counterexamples and data feasibility. No model implementation or biological-data analysis was performed.

## Sources and methods

Used web search for discovery, primary journal articles and PubMed/Europe PMC records for scientific claims, ClinicalTrials.gov API v2 for five trial records, and NCBI GEO text manifests for three datasets. Reviews were used as leads to original studies, not as proof of efficacy or novelty. No paid database or exhaustive patent search was performed. No authors were contacted.

Full-text XML retrieval used `https://www.ebi.ac.uk/europepmc/webservices/rest/{PMCID}/fullTextXML`. Core metadata used the REST search endpoint with DOI, PMID or title queries. Some PMC pages served browser challenges; ordinary PMC HTML access and Europe PMC provided alternative public access. Several XML requests returned 404 despite a listed PMCID. Inspection depth is recorded per source in SOURCES.json; retrieval alone is not described as full appraisal.

## Search families and representative exact queries

- `EGFR mutant lung cancer osimertinib persister residual disease 2025 2026 clinical trial TEAD`
- `osimertinib persister eradication ferroptosis GPX4 lung cancer primary study`
- `osimertinib residual disease combination FLAURA2 MARIPOSA overall survival 2025 trial`
- `EGFR lung cancer persister senescence immune clearance 2025 2026`
- `"osimertinib" "persister" "sequential" BCL`
- `"osimertinib" "persister" "sequential" "ferroptosis"`
- `"EGFR" "persister" "mitochondrial transfer" fibroblasts`
- `"Therapy-Induced Evolution" "Lung Cancer" Maynard 2020`
- `osimertinib TEAD inhibitor clinical trial VT3989 IK-930 IAG933`
- `osimertinib navitoclax phase Ib trial 27 patients`
- `"TROP2" "TEAD" "persister"`
- `"persister" "osimertinib" "state" "optimal" schedule model`
- `"Transient IGF-1R" osimertinib`
- `"osimertinib" "APOBEC3" persister Nature 2023`
- `"osimertinib" "dacomitinib" "computational model"`
- `"TROP2" "YAP" "osimertinib"`
- `"TROP2" "CAR-T" "Baldacci" "toxicity"`
- `"persister" "osimertinib" "combination" "cross-resistance"`
- `"persister" "TEAD" "TROP2" combination timing`
- `"osimertinib" "PD-1" "VEGFR2" sequential`
- `"osimertinib" "senolytic" persister`
- `"osimertinib" "TEAD" "ferroptosis" combination`
- `"osimertinib" "persister" "feedback" treatment`
- `"osimertinib" "persister" "lineage" "combination" 2026`
- `"EGFR" "TP53" "RB1" transformation osimertinib study`
- `"osimertinib" "TROP2" "AND" gate persister`
- `"osimertinib" "adaptive" "persister" "model" schedule`
- `"GSE193259" "prediction" combination`
- `"GSE150949" "GSE193259"`

DOI/PMID follow-ups resolved bibliographic identity for the selected sources. Broad title searches occasionally returned unrelated records; these were excluded from the curated source ledger. Search-engine relative publication ages were not treated as authoritative dates.

## Specific prior-art findings that changed the recommendation

1. Criscione already compares upfront versus sequential drug combinations across multiple models.
2. Wang already reports transient IGF-1R treatment followed by non-regrowth in selected mouse models.
3. Baldacci and Brea already investigate TROP2 CAR-T approaches and discuss logic-gated safety strategies.
4. Liao 2026 adds another TROP2 ADC and preliminary clinical combination evidence; the result with one ADC must not be generalized to the class.
5. ResMap, June 2026, already benchmarks 94 compounds/57 targets across genotypes and oxygen conditions, and reports both beneficial and pro-persistence effects. A generic atlas would duplicate that objective.

## Located but not fully appraised

- Poels-associated osimertinib–dacomitinib pharmacokinetic/tumour-evolution modeling: primary full-text excerpts located, but complete bibliographic and methodological appraisal not completed. No novel scheduling-model claim is made.
- *Optimal Dosing Strategies in Non-Small Cell Lung Cancer: A Multi-Scale Modelling Approach to Combat Drug Tolerance*: discovery record located; primary publication status and full content not verified. Excluded as scientific evidence.
- *Differential Models of Time-Variant Tumor Growth Trajectories with Sensitive, Persister, Resistant Cell Population in Lung Tumors During Tyrosine Kinase Inhibitor Therapy*: publisher result located; full page access failed. Excluded as evidence for performance, but blocks any confident claim that such model structure is new.
- Chen et al. ASCO 2024 abstract, DOI 10.1200/JCO.2024.42.16_suppl.8565: n=11 EGFR/TP53/RB1 study reported that early platinum/etoposide did not prevent transformation. Meeting abstract only; not used as definitive clinical efficacy evidence.

## Data-access checks

Retrieved GEO series manifests for GSE193259 (112 sample records), GSE150949 (19), GSE302284 (6). Counts refer to archive records, not independent biological replicates or patients. Processed matrix/metadata and raw archive listings were checked; underlying expression matrices were not analyzed. Source-data supplements were identified but not analyzed.

ClinicalTrials.gov records: NCT04665206, NCT05228015, NCT02520778, NCT04780568, NCT05401110. Selected fields are preserved in TRIAL_RECORDS.json. Trial registration does not establish that a combination arm enrolled participants or produced results.

## Limits on negative search conclusions

No universal statement that an intervention has never been tried is justified. This audit does not cover every patent, conference abstract, preprint, unpublished study, drug program or clinical registry. Absence of a located exact match is a provisional gap, not proof of novelty. New searches and full protocol comparison would be required before beginning the proposed BRD4 discrepancy investigation.

Stopping instruction: audit only. Further data analysis, modeling, simulations and outreach await Scott's authorization.
