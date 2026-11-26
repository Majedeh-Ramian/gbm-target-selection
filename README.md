# GBM Target Selection (in silico, adult glioblastoma)

Independent computational project. **Status: Step 1 (target selection) complete; docking pipeline in progress.**
No drug candidates have been identified. This repository ranks targets only.

## Goal
Rank 7 candidate targets (EGFR, PDGFRA, CDK4, CDK6, MET, MDM2, PIK3CA) for adult glioblastoma, to choose targets for small-molecule docking.

## Data sources
| Evidence | Source | How obtained |
|---|---|---|
| Structures | RCSB PDB, human, resolution < 3 Å | RCSB Search API (counts change over time; run in Nov 2026) |
| Alteration frequency | cBioPortal, GBM TCGA PanCancer Atlas (n=378, samples with mutation and CNA data) | Read from the web interface, entered as constants |
| Overall survival | cBioPortal, same study, single-gene queries, log-rank | Read from the web interface, entered as constants |
| Dependency | DepMap Public 26Q3, CRISPR Chronos | **Computed** from CRISPRGeneEffect.csv, Model.csv |
| TP53 status | DepMap 26Q3, OmicsSomaticMutations.csv | **Computed** |

Result files are in the repository root (the notebook reads them from a `results/` folder when run locally).

DepMap files are not included (large). Download them from the DepMap portal (release 26Q3) into `data/`.

## Scoring (0-2 per criterion, then summed)
- Structures: >=50 = 2, >=10 = 1, else 0
- Alteration frequency: >=20% = 2, >=8% = 1, else 0
- Survival: log-rank p < 0.05 = 2, else 0
- Dependency: many = 2, some = 1, none = 0 (adult-type diffuse glioma lines, n=106)

## Dependency (computed)
Number of the 106 adult-type diffuse glioma lines with Chronos < -1:

| Gene | Lines < -1 | Fraction | Label |
|---|---|---|---|
| EGFR | 4 | 0.038 | none |
| PDGFRA | 4 | 0.038 | some |
| CDK4 | 17 | 0.160 | some |
| CDK6 | 38 | 0.358 | many |
| MET | 3 | 0.028 | some |
| MDM2 | 31 | 0.292 | many |
| PIK3CA | 6 | 0.057 | some |

## Results
Scores: MDM2 = 5; EGFR = 4; PIK3CA = 4; MET, PDGFRA, CDK6 = 3; CDK4 = 2.
The survival criterion scored 0 for all genes (not discriminating).

Sensitivity analysis (changing structure and alteration thresholds, 7 settings): MDM2 was never below the others but tied for first in 2 of 7 settings; EGFR scored 4 in all 7. **MDM2 and EGFR are the two leading candidates; the final choice is deferred to re-docking results.**

## MDM2 and TP53
Among the 106 lines, MDM2 dependency (Chronos < -1) was seen in 31 lines, 30 of them TP53 wild-type.

| | MDM2 dependent | Not dependent | Total |
|---|---|---|---|
| TP53 wild-type | 30 | 8 | 38 |
| TP53 damaging | 1 | 67 | 68 |

Median MDM2 effect: -1.60 (TP53 wild-type) vs -0.32 (TP53 damaging). "Damaging" = LikelyLoF, Hotspot, or HIGH VEP impact (rule fixed before viewing results). Using "any TP53 mutation" instead gave 29/37 vs 2/69, so the result does not depend on the definition.
**MDM2 dependency is largely restricted to the TP53 wild-type subset.** The MDM2 score should be read with this caveat.

## Limitations
- Thresholds and dependency labels were chosen after seeing data. Sensitivity analysis varied structure and alteration thresholds only, not criterion weights or the dependency rule.
- Dependency labels (none/some/many) were first assigned by eye and later checked against computed fractions. Small inconsistency: EGFR and PDGFRA both have 4 lines below -1 but different labels (none vs some). Relabeling by fraction would change some scores; not done.
- Cell lines are not patients. The 106 lines are about 83% IDH-wildtype glioblastoma, the rest other adult-type diffuse gliomas. TP53-mutant lines may be over-represented relative to tumors.
- Alteration frequency and survival values were read from the cBioPortal web interface and entered as constants; they were spot-checked (EGFR: 53%, n=378) but not all re-verified.
- PDB counts: 288 on the website vs 285 in code for EGFR (strict "<" vs "<=" and release timing).
- Progression-free and disease-specific survival results were deliberately not reported (multiple testing, ~21 tests).
- Adult data only; pediatric high-grade glioma is not covered.

## Next steps
Reference inhibitors (ChEMBL/PDB), re-docking for MDM2 and EGFR, screening and ADMET/BBB filtering on the better-validated target.

## Environment
Python 3.14, pandas, numpy, matplotlib, seaborn, lifelines, scipy, statsmodels, jupyter.
## Step 2: Reference inhibitors for MDM2 and EGFR (ChEMBL 37)

Notebook: `02_reference_inhibitors.ipynb`

### Data source
- ChEMBL release 37 (SQLite dump, `chembl_37_sqlite.tar.gz`, file dated 2026-05-29).
- The database file (~28 GB extracted) is not included. To reproduce, download it from the ChEMBL FTP site and change `db_path` in the notebook.
- Targets: MDM2 (`CHEMBL5023`) and EGFR (`CHEMBL203`), both SINGLE PROTEIN entries.

### Why a local database
The ChEMBL web API returned repeated HTTP 500 errors for the full activity query, so only a partial download (460 of 6,340 rows) was possible. That partial data was discarded. All results here come from the SQLite dump. The row count for MDM2 from the local database (6,340) matched the count reported by the API.

### Query filters
`standard_type` in (IC50, Ki, Kd), `standard_relation` = "=", `standard_units` = "nM", `assay_type` = "B" (binding assays only).

### Cleaning rules (applied identically to both targets)
1. Remove rows without a SMILES.
2. Remove rows with a ChEMBL `data_validity_comment` (e.g. "Outside typical range").
3. Remove rows flagged `potential_duplicate`.
4. Remove rows with an empty `standard_value` (one EGFR row).

| Step | MDM2 rows | EGFR rows |
|---|---|---|
| Raw query | 6,340 | 19,826 |
| After removing no-SMILES | 6,336 | 19,809 |
| After removing validity warnings | 6,300 | 19,528 |
| After removing potential duplicates | 5,878 | 17,299 |
| After removing empty value | 5,878 | 17,298 |
| Unique molecules | 4,785 | 10,691 |

pActivity = 9 − log10(value in nM). It agreed with ChEMBL's own `pchembl_value` to within 0.0055 for all rows in both targets.

### Reference sets
- **Full actives:** median pActivity per molecule ≥ 7 (MDM2: 3,240; EGFR: 5,499).
- **Strict actives:** full actives with at least 2 measurements and spread (max − min pActivity) ≤ 1 (MDM2: 441; EGFR: 1,054).
- **Known drugs** are kept separately as positive controls (see `*_known_drugs_check.csv`).

### Findings and limitations
- The thresholds (cleaning rules, pActivity ≥ 7, n ≥ 2, spread ≤ 1) were chosen after looking at the data, not before.
- Measurements disagree a lot between assays: 28% of MDM2 molecules and 42% of EGFR molecules with ≥ 2 measurements have a spread > 1 pActivity unit. For example, gefitinib spans 4.6 to 10.0 across ten different assays at its extremes.
- The strict set excludes all known drugs checked (nutlin-3, idasanutlin, navtemadlin, milademetan, erlotinib, gefitinib, afatinib, osimertinib), because these have wide spread or only one measurement. It is therefore not a good reference set on its own.
- Most molecules have a single measurement (MDM2: 82%; EGFR: 69%).
- IC50, Ki and Kd are pooled. IC50 depends on assay conditions.
- The EGFR data mixes wild-type and mutant assays (partially flagged by ChEMBL `variant_id`) and possibly other differences I did not investigate. Mutation status matters for GBM, so this needs care.
- Median pActivity is high (MDM2 7.79, EGFR 7.15), which reflects publication bias toward optimized compounds.
- 23 EGFR rows have pActivity ≥ 10.5 and were kept; I could not verify them.
- Some salt forms (milademetan tosylate, osimertinib mesylate) have no data under their own IDs, so I did not merge salts with parent compounds.
- Not done yet: PDB ligands, re-docking, screening, ADMET.
