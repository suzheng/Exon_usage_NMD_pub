# Exon Usage in NMD Genes

# Pipeline and Analysis Scripts

This repository contains the scripts and notebooks used for data processing, analysis, and visualization for our study in NMD gene exon usage. The code processes RNA-seq data to compute Percent Spliced In (PSI) values, performs quality control, and conducts downstream analyses, including exon usage studies and correlations with variant data.

## Directory Structure

```
pipeline/
├── config/
│   └── config.yaml
├── workflow/
│   ├── Snakefile.psi_quantification_bam_input.SE.-RP0001.py
│   ├── Snakefile.psi_quantification_bam_input-RP0001.py
│   ├── Snakefile.rnaseqc.py
│   ├── Snakefile.fastqc.py
│   └── scripts/
│       ├── generate_junction_bed.sh
│       └── PSI.sh
notebooks/
├── shared.R
├── 1_get_age_dependent_exon_genes-generate_psi_meta_dfs.ipynb
├── 2_get_age_dependent_exon_genes-PCAs.ipynb
├── 3_get_age_dependent_exon_genes-correct_batch_effect.ipynb
├── 4_get_age_dependent_exon_genes-PCAs-IR_corrected.ipynb
├── 5_get_age_dependent_exon_genes-assoc_analysis.ipynb
├── 6_get_age_tissue_subtype_dependent_exons.ipynb
├── 7_visualize_age_tissue_subtype_dependent_exons.ipynb
├── 8_exon_usage_cor_clinvar_target_genes-cal_var_count_every_exon.ipynb
├── 9_exon_usage_cor_clinvar_target_genes-cor_analysis.ipynb
├── combine_plots.ipynb
```

## Configuration

### `config/config.yaml`

- **Purpose**: Stores configuration settings for the pipeline.
- **Key Settings**:
  - `ref_genome`: Path to the reference genome index for the STAR aligner.
  - `SI`: Directory path for Singularity images.

## Workflows

### `workflow/Snakefile.psi_quantification_bam_input.SE.-RP0001.py`

- **Purpose**: Defines the workflow for PSI quantification using single-end BAM inputs.
- **Main Components**:
  - Converts BAM files to FASTQ.
  - Aligns reads using STAR.
  - Calculates PSI values.
  - Performs quality control checks.

### `workflow/Snakefile.psi_quantification_bam_input-RP0001.py`

- **Purpose**: Similar to the above but tailored for paired-end BAM inputs.
- **Main Components**:
  - Converts paired-end BAM files to FASTQ.
  - Aligns paired-end reads using STAR.
  - Calculates PSI values.
  - Performs quality control checks.

### `workflow/Snakefile.rnaseqc.py`

- **Purpose**: Runs RNA-SeQC for quality assessment of aligned BAM files.
- **Main Components**:
  - Processes sorted BAM files.
  - Generates RNA-SeQC metrics reports.

### `workflow/Snakefile.fastqc.py`

- **Purpose**: Executes FastQC for quality control of FASTQ files.
- **Main Components**:
  - Runs FastQC on FASTQ files.
  - Generates FastQC HTML reports.

## Scripts

### `workflow/scripts/generate_junction_bed.sh`

- **Purpose**: Generates BED files for junctions from STAR output.
- **Functionality**:
  - Processes STAR's `SJ.out.tab` to create BED-formatted junction files.

### `workflow/scripts/PSI.sh`

- **Purpose**: Calculates the Percent Spliced In (PSI) metric.
- **Functionality**:
  - Counts exon inclusion and exclusion.
  - Filters junctions.
  - Computes PSI values based on inclusion and exclusion counts.

## Notebooks

### `notebooks/shared.R`

- **Purpose**: Sets up shared configurations and loads necessary R libraries for data analysis and visualization. Defines plotting themes and styles for consistent appearance across figures.

### `notebooks/1_get_age_dependent_exon_genes-generate_psi_meta_dfs.ipynb`

- **Purpose**: Processes PSI data to generate metadata and PSI matrices. Filters out low-quality samples and exons with excessive missing values.

### `notebooks/2_get_age_dependent_exon_genes-PCAs.ipynb`

- **Purpose**: Performs Principal Component Analysis (PCA) on the PSI data to identify patterns and potential batch effects. Visualizes differences in exon usage across samples.

### `notebooks/3_get_age_dependent_exon_genes-correct_batch_effect.ipynb`

- **Purpose**: Corrects for batch effects in the inclusion and exclusion counts of exons using the `limma` package. Adjusts the data for covariates like intronic rates.

### `notebooks/4_get_age_dependent_exon_genes-PCAs-IR_corrected.ipynb`

- **Purpose**: Re-runs PCA on the intronic rate-corrected PSI data to assess the effectiveness of batch effect correction and explore biological variance.

### `notebooks/5_get_age_dependent_exon_genes-assoc_analysis.ipynb`

- **Purpose**: Conducts association analysis to identify exons with significant differential usage between age groups or tissue types. Applies statistical tests to determine significant correlations.

### `notebooks/6_get_age_tissue_subtype_dependent_exons.ipynb`

- **Purpose**: Analyzes exon usage dependencies based on age groups and muscle tissue subtypes (skeletal and cardiac). Identifies significant changes in exon usage across developmental stages and tissue types.

### `notebooks/7_visualize_age_tissue_subtype_dependent_exons.ipynb`

- **Purpose**: Visualizes relationships and overlaps in delta PSI values between original and intronic rate-adjusted data. Generates Venn diagrams and combined plots for age and tissue subtype-dependent exon usage.

### `notebooks/8_exon_usage_cor_clinvar_target_genes-cal_var_count_every_exon.ipynb`

- **Purpose**: Correlates exon usage with clinical variant data from ClinVar and population variant data from gnomAD. Calculates the density of various variant types within each exon.

### `notebooks/9_exon_usage_cor_clinvar_target_genes-cor_analysis.ipynb`

- **Purpose**: Performs correlation analyses between exon usage metrics and variant densities from ClinVar and gnomAD. Generates plots to visualize these correlations and assesses statistical significance.

### `notebooks/combine_plots.ipynb`

- **Purpose**: Combines and arranges all generated plots into comprehensive figures for reporting. Enhances plot aesthetics and ensures consistency across all visualizations.

## Usage

To execute the pipeline, navigate to the `workflow` directory and run the appropriate Snakemake command, specifying the desired Snakefile. Ensure that the `config.yaml` file is properly configured with the correct paths.

Example command for running a Snakefile:

```bash
snakemake --cores all -F --use-singularity --singularity-args "-B DIR_OUTSIDER:DIR_DOCKER" -C out_dir=$PWD meta_dir=meta_ERP003613_SRP028336 -s Snakefile.psi_quantification_bam_input.SE.-RP0001.py TARGET_FILE
```

Replace `TARGET_FILE` with the target file you wish to generate.

## Dependencies

- **Software**:
  - [Snakemake](https://snakemake.readthedocs.io/) for workflow management.
  - [STAR](https://github.com/alexdobin/STAR) aligner for RNA-seq alignment.
  - [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) for FASTQ quality control.
  - [RNA-SeQC](https://github.com/getzlab/rnaseqc) for RNA-Seq quality assessment.
  - R and various R packages (e.g., `ggplot2`, `dplyr`, `limma`) for data analysis and visualization in the notebooks.

- **Data**:
  - Reference genome indexed for STAR aligner.
  - Junction BED files generated from STAR outputs for PSI calculation.

## Notes

- Ensure all paths in `config/config.yaml` are correctly set to point to your data and reference files.
- The notebooks assume that the PSI data has been generated and is accessible in the expected directories.
- The `combine_plots.ipynb` notebook requires all the plot objects saved from previous analyses to generate the final figures.

---