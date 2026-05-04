# Comparative oncogenomics of canine and feline mammary cancer reveals a conserved and tractable oncogenic core with human breast cancer concordance

![R](https://img.shields.io/badge/R-%E2%89%A54.0-blue)
![License](https://img.shields.io/badge/license-MIT-green)
![DOI](https://img.shields.io/badge/DOI-pending-lightgrey)

This repository contains a standalone, reproducible R pipeline for comparative oncogenomic analysis of canine and feline mammary cancer. The workflow integrates canine transcriptome analysis, feline mammary cancer gene signatures, ortholog mapping into a common human reference space, conserved Hallmark pathway analysis, Open Targets/ChEMBL translational annotation, drug-level/probe-linked indication expansion, target prioritisation, and optional validation analyses.

## Project aim

The pipeline was developed to identify conserved oncogenic programs and conserved actionable targets across canine and feline mammary cancer, with additional translational annotation using human druggability and clinical-candidate resources.

The primary analysis identifies:

- directionally conserved Hallmark programs across species;
- a strict conserved UP-regulated target core;
- druggability and clinical-candidate evidence for conserved targets;
- integrated drug-level and probe-linked indication evidence;
- a consensus-prioritised target shortlist;
- translational evidence tiers;
- optional network-context validation using STRING;
- optional human external concordance using TCGA-BRCA.

## Repository structure

```text
.
├── 00_run_pipeline.R
├── 00_install_packages.R
├── LICENSE
├── METHODS_LANGUAGE.md
├── README.md
├── R/
│   ├── 00_config.R
│   ├── 01_module_dog_gsea.R
│   ├── 02_module_cat_hallmarks.R
│   ├── 03_module_conserved_hallmarks.R
│   ├── 04_module_conserved_targets.R
│   ├── 05_module_drug_targets.R
│   ├── 05b_module_probe_linked_indications.R
│   ├── 06_module_target_prioritisation.R
│   ├── 06b_prioritisation_sensitivity_analysis.R
│   ├── 07_module_outputs.R
│   ├── 08_figure1_conserved_hallmarks.R
│   ├── 09_figure2_overlap_support.R
│   ├── 10_figure3_target_prioritisation.R
│   ├── 11a_build_background_tractability_universe_light.R
│   ├── 11b_permutation_validation_light.R
│   ├── 13_network_validation_prioritised_targets.R
│   ├── 14_module_tcga_brca_external_concordance.R
│   └── 15_figure4_tcga_brca_concordance.R
├── data/
├── cache/
├── results/
└── figures/
```

## Data availability and required input files

Raw input files are not included in this repository because of file size and/or data-licensing constraints. The canine transcriptome dataset is derived from GEO accession `GSE119810`. The feline upregulated and downregulated gene lists were curated from published feline mammary cancer study cited in the manuscript. These literature-derived gene signatures constitute the feline input for the comparative analysis.

The primary pipeline expects the following files under `data/`:

```text
data/GSE119810_CMT_222S_FPKM.csv
data/CAT_up.txt
data/CAT_down.txt
```

For optional TCGA-BRCA validation, local GDC-derived TCGA-BRCA STAR-count files must be available in the directory specified by `PATH_CONFIG$tcga_dir` in `R/00_config.R`:

```text
data/TCGA_BRCA_STAR_Counts/
```

On Windows, use a short project path for TCGA files, because very long paths can cause file-reading failures even when the files are present.

## R dependencies

Core CRAN packages:

```r
dplyr
readr
tidyr
stringr
ggplot2
httr
jsonlite
purrr
tibble
```

Core Bioconductor packages:

```r
limma
biomaRt
fgsea
msigdbr
clusterProfiler
BiocParallel
```

Optional CRAN packages:

```r
igraph
ggrepel
```

Optional Bioconductor packages:

```r
STRINGdb
edgeR
TCGAbiolinks
```

To install the main pipeline dependencies, run:

```r
source("00_install_packages.R")
```

To install optional network and TCGA validation dependencies as well, run:

```r
install_optional_packages()
```

To install all core and optional dependencies at once, run:

```r
source("00_install_packages.R")
install_all_packages()
```

## How to run the pipeline

From the project root directory, run:

```r
source("00_run_pipeline.R")
```

Execution is controlled from:

```text
R/00_config.R
```

Main execution flags are defined in `EXECUTION_CONFIG`. Core analysis steps are enabled by default. Optional validation modules are disabled by default and can be enabled in `EXECUTION_CONFIG`.

## Main pipeline steps

### Step 1 — Dog transcriptome analysis and Hallmark GSEA

The canine RNA-seq expression matrix is mapped to human Ensembl gene space using high-confidence 1:1 dog-human orthologs from Ensembl release 104. Differential expression is performed using limma with block-aware modelling of paired tumor/normal samples where available.

The default dog analysis mode is:

```r
dog_analysis_mode = "all_blockaware"
```

An alternative sensitivity mode is available:

```r
dog_analysis_mode = "complete_pairs_only"
```

Hallmark enrichment is performed using `fgsea`. The pipeline uses seeded, serial execution for reproducibility.

Key outputs:

```text
results/DOG_sample_structure_QC.csv
results/DOG_mapping_1to1_highconf.rds
results/DOG_expr_humanENSG_1to1.rds
results/DOG_DE_limma_all_blockaware_full.csv
results/DOG_DE_limma_all_blockaware.rds
results/DOG_DE_limma_all_blockaware_sig_FDR0.05_logFC1.csv
results/DOG_Hallmark_fgsea.rds
results/DOG_Hallmark_fgsea_simple_plus.csv
results/DOG_up_pathways_padj0.05.txt
results/DOG_down_pathways_padj0.05.txt
results/DOG_up_leadingEdge_ENSG.txt
results/DOG_down_leadingEdge_ENSG.txt
```

### Step 2 — Cat gene-list mapping and Hallmark ORA

Feline upregulated and downregulated gene lists are mapped to human gene symbols. The ORA universe is fixed as the intersection between the canine-derived human ortholog universe and MSigDB Hallmark genes.

Hallmark ORA is performed using `clusterProfiler::enricher()`.

Important reproducibility note: `pvalueCutoff = 1`, `qvalueCutoff = 1`, and `pAdjustMethod = "BH"` are set explicitly so that the full tested Hallmark result table is returned. Filtering is then performed downstream using the configured FDR threshold.

Key outputs:

```text
results/CAT_up_mapped_humanSymbols.txt
results/CAT_down_mapped_humanSymbols.txt
results/CAT_up_universeFixed.txt
results/CAT_down_universeFixed.txt
results/CAT_Hallmark_ORA_UP_universeFixed.csv
results/CAT_Hallmark_ORA_DOWN_universeFixed.csv
```

### Step 3 — Directionally conserved Hallmarks

Dog and cat Hallmark results are intersected by direction. A Hallmark is considered conserved when it is significant in both species and shows the same direction of regulation.

Key output:

```text
results/CONSERVED_HALLMARKS_summary.csv
```

### Step 4 — Strict conserved core targets

For each conserved UP Hallmark, dog leading-edge genes are intersected with the universe-filtered cat UP genes. The strict conserved UP core is defined as the union of these Hallmark-wise cross-species intersections.

Key outputs:

```text
results/Conserved_Core_UP_genes_strict.txt
results/Conserved_Core_DOWN_genes_strict.txt
results/core_conserved_targets.txt
results/Conserved_UP_core_gene_frequency.csv
results/CONSERVED_CORE_TARGETS_by_hallmark.csv
results/DOG_LE_vs_CAT_UP_Fisher_byHallmark.csv
results/DOG_LE_vs_CAT_DOWN_Fisher_byHallmark.csv
results/CONSERVED_HALLMARKS_strict_summary.csv
```

### Step 5 — Druggability and clinical-candidate annotation

Conserved upregulated core targets are queried against the Open Targets Platform GraphQL API. The module retrieves target-level drug/clinical-candidate evidence and tractability annotations, then queries ChEMBL for ATC classifications when drug identifiers are available.

The pipeline uses the current Open Targets GraphQL `target(ensemblId: ...)` query with the `drugAndClinicalCandidates` and `tractability` fields. Because Open Targets is a live resource and its API schema may change over time, raw JSON responses are cached locally.

Open Targets cache:

```text
cache/OT_cache_v3_drugAndClinicalCandidates/
```

The Open Targets access date is not hardcoded in this README. The final run timestamp is written automatically to:

```text
results/REPRODUCIBILITY_RECORD.txt
```

Key outputs:

```text
results/OT_query_QC.csv
results/OT_query_error_log.csv
results/drug_target_table_raw_OpenTargets.csv
results/OT_phase_parsing_diagnostic.csv
results/drug_target_table_with_ATC_and_fallback.csv
results/tractability_table_raw.csv
results/TRACTABILITY_summary_by_target_fixed.csv
results/conservation_score_targets.csv
results/target_drug_summary_with_indication_counts.csv
results/target_master.csv
results/target_master_primary_targetlevel.csv
```

### Step 5B — Drug-level and probe-linked indication expansion

Step 5B is part of the main workflow. It adds an additional Open Targets evidence layer using drug-level indication expansion and high-quality chemical probe-linked drug evidence.

The integrated evidence is built from the union of:

```text
primary target-level drug/candidate rows
primary drug-level indication expansion
high-quality probe-linked drug-level indication expansion
```

By default, Step 5B overwrites `target_master.csv` with integrated evidence columns so that Step 6 and Figure 3 use the updated translational annotation layer.

This behaviour is controlled in `R/00_config.R`:

```r
run_step5b_probe_indications = TRUE
probe_overwrite_target_master = TRUE
probe_min_phase = 1
```

Key outputs:

```text
results/STEP5B_method_parameters.csv
results/OT_chemical_probe_links_all_unfiltered.csv
results/OT_chemical_probe_links_high_quality_only.csv
results/OT_probe_drug_phase_check.csv
results/OT_probe_drugs_excluded_by_phase_filter.csv
results/OT_chemical_probe_drug_links.csv
results/OT_target_drug_links_for_deep_indication_query.csv
results/OT_drug_level_indications_raw.csv
results/OT_target_drug_deep_indication_rows.csv
results/OT_primary_target_level_rows_for_integration.csv
results/OT_integrated_drug_indication_rows.csv
results/target_deep_indication_summary.csv
results/target_master_integrated_deep_indications.csv
results/QC_targets_with_phase_gain_from_deep_indications.csv
results/QC_targets_with_breast_gain_from_deep_indications.csv
results/QC_targets_with_oncology_gain_from_deep_indications.csv
results/target_master.csv
```

### Step 6 — Consensus target prioritisation

Target prioritisation uses a consensus point-based framework combining rank-derived evidence layers and binary modality evidence.

Quantitative/count-based features are converted into rank-based points using a top-N tie-aware scheme:

```text
conservation_score
max_drug_phase
n_drugs
n_indications
n_oncology_drug_indications
n_breast_drug_indications
```

Small-molecule tractability is included as a binary support flag rather than ranked, because it represents the availability of a therapeutically relevant modality class rather than an ordinal quantitative variable.

The final consensus priority score is calculated as the sum of these feature-specific points. The exact top-N shortlist is defined by `priority_rank`; by default, `top_n_targets = 15`.

Key outputs:

```text
results/target_prioritisation_scored.csv
results/shortlist_overall_top15_exact.csv
results/shortlist_onco_or_approved_all.csv
results/shortlist_overall_inclusive_top15_cutoff.csv
results/shortlist_onco_or_approved_inclusive_top15_cutoff.csv
```

Note: legacy filenames containing `inclusive_top15_cutoff` are retained as compatibility aliases. The current prioritisation and network validation use the exact top-N targets by `priority_rank`.

### Step 6B — Prioritisation sensitivity analysis

This optional module compares primary target-level evidence from Step 5 with the integrated Step 5B evidence layer.

Enable it in `R/00_config.R`:

```r
run_sensitivity_prioritisation = TRUE
```

Relevant module:

```text
R/06b_prioritisation_sensitivity_analysis.R
```

### Step 7 — Summary outputs and translational evidence tiers

The output module creates manuscript-ready target tables and translational evidence tier summaries.

Translational evidence tiers:

```text
Approved/phase >=2 evidence
Early oncology-linked evidence
Tractable, no retained drug evidence
Biology only
```

Key outputs:

```text
results/PIPELINE_METADATA.csv
results/TOP_TARGETS_FORMATTED.csv
results/TRANSLATIONAL_GAP_summary.csv
results/TRANSLATIONAL_GAP_target_tiers.csv
results/TRANSLATIONAL_GAP_tier_counts.csv
results/TRANSLATIONAL_GAP_clinical_leverage_counts.csv
results/SHORTLIST_drug_evidence.csv
results/PIPELINE_SUMMARY.txt
```

### Step 8 — Figure 1

Generates the conserved Hallmark program figure panels.

Key outputs:

```text
results/FIG1_panels/
results/FIG1_tables/
```

### Step 9 — Figure 2 and Table 1

Generates gene-level overlap support panels and the conserved Hallmark summary table.

Key outputs:

```text
results/TABLE1_conserved_hallmark_programs.csv
results/FIG2_panels/
results/FIG2_tables/
```

### Step 10 — Figure 3

Generates target prioritisation and translational landscape panels.

Key outputs:

```text
results/FIG3_panels/
results/FIG3_tables/
```

## Optional validation modules

Optional validation steps are controlled from `EXECUTION_CONFIG` in:

```text
R/00_config.R
```

### Background tractability validation

The background validation module constructs a biologically constrained comparison universe from Hallmark-linked genes.

Relevant modules:

```text
R/11a_build_background_tractability_universe_light.R
R/11b_permutation_validation_light.R
```

### Network validation using STRINGdb

Network validation uses STRING v12.0 protein interaction data through the `STRINGdb` workflow and local STRING cache. The first run requires internet access to download STRING mapping and protein-link files; subsequent runs reuse cached files.

The network comparison uses the exact top-N prioritised targets, defined by `priority_rank` order from:

```text
results/target_prioritisation_scored.csv
```

By default, N = 15.

The network validation tests whether prioritised targets show higher network centrality than expected by chance. Random target sets of equal size are sampled without replacement from the mapped conserved-core network. Empirical one-sided p-values are calculated as the proportion of random sets with mean centrality greater than or equal to the observed mean centrality of the prioritised set. P-values across tested network metrics are adjusted using the Benjamini-Hochberg procedure.

Primary network metrics:

```text
degree
weighted_degree
betweenness
harmonic_centrality
eigenvector_centrality
```

Relevant module:

```text
R/13_network_validation_prioritised_targets.R
```

Key outputs:

```text
results/NETWORK_validation/tables/NETWORK_summary.csv
results/NETWORK_validation/tables/NETWORK_group_summary.csv
results/NETWORK_validation/tables/NETWORK_permutation_summary_with_BH.csv
results/NETWORK_validation/tables/NETWORK_permutation_summary_manuscript_ready.csv
results/NETWORK_validation/tables/NETWORK_method_parameters.csv
results/NETWORK_validation/plots/
results/NETWORK_validation/NETWORK_STRING_core_graph.rds
```

### TCGA-BRCA external concordance

The TCGA validation module uses locally downloaded GDC/TCGA-BRCA STAR-count files. The TCGA directory must be specified in `PATH_CONFIG$tcga_dir`.

Relevant modules:

```text
R/14_module_tcga_brca_external_concordance.R
R/15_figure4_tcga_brca_concordance.R
```

Key outputs:

```text
results/TCGA_BRCA_core_gene_results.csv
results/TCGA_BRCA_concordance_summary.csv
results/FIG4_TCGA_A_core_concordance_dotplot.png
results/FIG4_TCGA_B_sigconc_barplot.png
results/FIG4_TCGA_summary.txt
```

## Reproducibility notes

The pipeline is designed to improve reproducibility through:

- fixed Ensembl release for ortholog mapping;
- cached biomaRt mappings;
- cached Open Targets JSON responses;
- cached ChEMBL ATC queries;
- explicit ORA cutoffs in `clusterProfiler::enricher()`;
- seeded serial execution for `fgsea`;
- recorded Open Targets query diagnostics;
- cached STRINGdb mapping/interactions;
- written pipeline metadata and session information;
- exported intermediate tables for each major analysis stage.

Important live resources:

```text
Open Targets Platform GraphQL API
ChEMBL API
STRINGdb
Ensembl BioMart
MSigDB Hallmark gene sets through msigdbr
GDC/TCGA
```

Because these resources may update over time, cache files and final output tables should be retained with the submitted code archive.

## Expected key final outputs

Main manuscript-relevant outputs include:

```text
results/CONSERVED_HALLMARKS_summary.csv
results/CONSERVED_HALLMARKS_strict_summary.csv
results/CONSERVED_CORE_TARGETS_by_hallmark.csv
results/Conserved_Core_UP_genes_strict.txt
results/target_master_primary_targetlevel.csv
results/target_master_integrated_deep_indications.csv
results/target_master.csv
results/target_prioritisation_scored.csv
results/TOP_TARGETS_FORMATTED.csv
results/TRANSLATIONAL_GAP_summary.csv
results/TRANSLATIONAL_GAP_target_tiers.csv
results/TABLE1_conserved_hallmark_programs.csv
results/FIG1_panels/
results/FIG2_panels/
results/FIG3_panels/
results/REPRODUCIBILITY_RECORD.txt
```

Optional validation outputs include:

```text
results/NETWORK_validation/
results/TCGA_BRCA_core_gene_results.csv
results/TCGA_BRCA_concordance_summary.csv
results/FIG4_TCGA_A_core_concordance_dotplot.png
results/FIG4_TCGA_B_sigconc_barplot.png
```

## Troubleshooting

### Open Targets returns no drug rows or API errors

Check:

```text
results/OT_query_QC.csv
results/OT_query_error_log.csv
```

If Open Targets has changed its GraphQL schema, update field names in:

```text
R/05_module_drug_targets.R
R/05b_module_probe_linked_indications.R
R/11a_build_background_tractability_universe_light.R
```

### Network validation fails because optional packages are missing

Install optional packages:

```r
source("00_install_packages.R")
install_optional_packages()
```

### Network validation download is interrupted

Delete the incomplete STRING protein-links file from the local STRING cache and rerun the network module. The first run requires a stable internet connection; subsequent runs use the cache.

### PDF export fails for network plots

PNG output is saved by default. PDF export can fail on Windows when the path is very long or when a PDF is already open in a viewer. Use `save_pdf = FALSE` or shorten the project path.

### TCGA validation cannot find files

Check that `PATH_CONFIG$tcga_dir` points to the local folder containing downloaded GDC/TCGA-BRCA STAR-count files. On Windows, avoid deeply nested project paths.

## License

This repository is released under the MIT License. See `LICENSE` for details.

The license applies to the code in this repository. It does not apply to third-party datasets or downloaded database files.

## Code availability

The source code for this pipeline is available at:

https://github.com/Jasonpat/comparative-oncogenomics-canine-feline-mammary-cancer


## Citation

Patergiannakis I.-S., Pappas IS. Comparative oncogenomics of canine and feline mammary cancer reveals a conserved and tractable oncogenic core with human breast cancer concordance. [Journal], 2026. DOI: pending.
