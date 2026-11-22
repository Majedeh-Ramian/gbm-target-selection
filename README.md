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
