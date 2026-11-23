# Longitudinal gut microbiome stability in type 2 asthma

R code for the analyses in:

> *Longitudinal gut microbiome stability in type 2 asthma: a prospective cohort study* (manuscript in preparation)

Sequencing data: NCBI Sequence Read Archive, BioProject **PRJNA1565519**.
Supplementary Table S9 maps each SRA sample to its de-identified subject, visit and sequencing run.

## Study design

- 44 adults with type 2 asthma, stool sampled at two visits (Time 1, Time 2)
- 41 non-asthma controls, sampled once
- 16S rRNA gene V3–V4 amplicon sequencing (Illumina MiSeq, 2 × 300 bp), three sequencing runs
- ASVs inferred with DADA2; genus, family and phylum tables obtained by summing ASV counts per taxon

## Analysis steps

Run the scripts in order from this folder. Each script reads the outputs of the previous ones.

| Script | Analysis | Manuscript output |
|---|---|---|
| `R/00_setup.R` | Paths, packages, shared settings | – |
| `R/01_prepare_data.R` | Taxon filtering (mean relative abundance ≥ 0.1% or ≥ 0.5%, prevalence ≥ 10%); CLR (pseudo-count 0.5) | input for all steps |
| `R/02_diversity.R` | Alpha diversity; Bray–Curtis PCoA, PERMANOVA; within- vs between-patient dissimilarity; distance to controls | Figure 1A–B, Figure 2A–D |
| `R/03_stability.R` | ICC per genus; CLR change per genus; community types (PAM, Aitchison distance) | Figure 2E–F, Figure 3 |
| `R/04_differential_abundance.R` | Patients vs controls, MaAsLin2 (age + sex; + BMI) | Figure 1D, Figure E4B–D, Table S2 |
| `R/05_clinical_associations.R` | Genus–clinical associations: pooled mixed-effects model and baseline-adjusted change model; 0.5% threshold; treatment subgroups | Figure 4, Figures E7–E8, Tables S6–S7 |
| `R/06_networks.R` | SPIEC-EASI co-abundance networks, GeneNet co-change networks, modules, module–clinical associations | Figure 5, Figure E9, Tables S3–S5 |
| `R/07_sensitivity.R` | Excluding 9 previously reported samples; sample map; run-adjusted differential abundance | Tables S8–S10 |

To run everything at once:

```bash
Rscript run_all.R
```

All outputs are written to `results/` (tables, figures, MaAsLin2 output, intermediate objects).

Not yet included: Figure 1C (LEfSe, run on Galaxy; `07_sensitivity.R` writes the LEfSe input format and reproduces the Kruskal–Wallis step in R), and the supplementary panels in Figures E1–E3, E4A, E5 and E6 (batch effect, treatment step, Aitchison PCoA, phylum composition, clustering validation).

## Input data

Individual-level clinical data cannot be shared publicly because of participant privacy; they are available from the corresponding author on reasonable request with institutional review board approval. The scripts expect the following files in `../01_data_objects/`:

| File | Content |
|---|---|
| `ASV_table_verified_UPDATED_20260617.rds` | ASV count matrix (samples × ASVs) after DADA2 quality control |
| `taxonomy_verified_UPDATED_20260617.rds` | ASV taxonomy |
| `metadata_verified_UPDATED_20260617.rds` | Sample metadata: `sample_id`, `subject_name_cohort`, `analysis_group` (Control / Asthma Time1 / Asthma Time2), `age`, `sex`, `bmi`, `sequencing_time` |
| `CRF_clinical_data_clean_UPDATED_20260617.rds` | Clinical data per patient at both visits: ACT, FEV1 % predicted, FEV1/FVC, FeNO, blood eosinophils, GINA treatment step |
| `genus/family/phylum_count_all_UNFILTERED_20260922.rds` | Taxon count tables (ASV counts summed per taxon) |

The ASV table can be regenerated from the SRA reads with DADA2.

## Key methods

- **CLR:** pseudo-count 0.5, applied after filtering, separately for each taxonomic level and threshold.
- **MaAsLin2:** pre-computed CLR values as input (`normalization = "NONE"`, `transform = "NONE"`). FDR (Benjamini–Hochberg) is applied over taxa within each outcome or group term, not across all model terms.
- **SPIEC-EASI:** filtered genus **counts** are supplied, because `spiec.easi()` applies its own CLR transformation (`clr(x + 1)`).
- **Significance symbols in figures:** \* P < 0.05, \*\* P < 0.01, \*\*\* P < 0.001 (unadjusted); † FDR q < 0.05.
- **Random seed:** 42 for all permutation, clustering and network steps.

## Reproducibility

Every script ends with a check against the results reported in the manuscript. In our runs, all tables and statistics were reproduced exactly (differences ≤ 1e-14, i.e. floating-point rounding). Final figure layouts were assembled from the panels in `results/figures/`.

## Software

R 4.x with: Maaslin2, vegan, permute, compositions, irr, cluster, SpiecEasi, GeneNet, corpcor, igraph, ggraph, tidygraph, ggforce, ggalluvial, ggplot2, patchwork, ggtext, dplyr, tidyr, readr. Exact versions are in `sessionInfo.txt`.

## License

MIT