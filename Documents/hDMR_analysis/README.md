<div align="center">
    <img src="../../images/hDMRs.jpg" alt="hDMR analysis" width="800"/>
</div>


### Haplotype-based Differentially Methylated Region (hDMR) Analysis

This module focuses on identifying and characterizing **haplotype-based differentially methylated regions (hDMRs)** across individuals and populations.  
The analysis aims to capture allele-specific methylation differences that may underlie imprinting control and haplotype-biased regulation.

- **hDMR Detection**:  
  For each population, **10 individuals** were randomly selected to perform hDMR detection using **Metilene (v0.2.8)** with the same parameters applied in pDMR identification — a minimum number of CpGs ≥ 5, a maximum CpG distance ≤ 1000 bp, and detection mode = 1 (*de-novo*).  
  The analysis was repeated **five times**, followed by a full-cohort run using all individuals.  
  We further required **segment length ≥ 100 bp** to improve robustness.  
  The complete set of hDMRs was defined as regions with mean methylation difference (α) > 0 and q-value < 1, and filtered using the official **Metilene** post-processing script (`metilene.output.pl`) to retain those with **q-value ≤ 0.05**, **α ≥ 0.1**, and **length ≥ 50 bp**.

- **Background Distribution**:  
  The genome was partitioned into **800 bp intervals** (corresponding to the average length of hDMRs, 784 bp).  
  For each interval, methylation differences between haplotype 1 and haplotype 2 were calculated to generate a background distribution, representing baseline haplotype variation across the genome.

- **Union hDMR and ICR Definition**:  
  A union set of hDMRs across all individuals was created by merging overlapping regions based on their maximum genomic span.  
  Regions overlapping with gene bodies or their **2 kb upstream** promoter sequences were defined as potential **imprinting control regions (ICRs)**.  
  For each merged ICR, Δ<sub>hDMR</sub> was defined as the arithmetic mean of absolute methylation differences between haplotypes across all constituent hDMRs:

<p align="center">
  Δ<sub>hDMR</sub> = mean(|Meth<sub>hap1</sub> − Meth<sub>hap2</sub>|)
</p>

  The **LenPro** parameter (hDMR density per gene) was calculated as the total hDMR length within a gene divided by the gene length.

- **Identification of Known and Candidate Imprinting Genes**:  
  ICRs were classified into two groups based on overlap with public datasets and population-specific criteria:
  - **Known imprinting genes** — Δ<sub>hDMR</sub> ≥ 0.3 and overlap with the **GENEIMPRINT** database ([https://www.geneimprint.com](https://www.geneimprint.com)) by ≥ 50%.  
  - **Candidate imprinting genes** — population-specific ICRs with **LenPro ≥ 0.1**, **Δ<sub>hDMR</sub> ≥ 0.5**, and **no overlap** with GENEIMPRINT.  

  These candidate ICRs were visualized in *Supplementary Fig. S11* and provide insight into population-specific allele-biased methylation and imprinting potential.

---

#### Summary of Scripts

| Script Name | Function | Description |
|--------------|-----------|-------------|
| `haplotype_extract.sh` | Prepare haplotype-level methylation input | Separates phased methylation by haplotype from modBAM files |
| `hDMR_metilene.sh` | Identify haplotype-based DMRs using Metilene | Applies pDMR parameters with segment length ≥ 100 bp |
| `merge_hDMR.sh` | Merge hDMRs across individuals | Generates a union set of hDMRs by merging overlapping regions |
| `calculate_Delta_hDMR.py` | Compute Δ<sub>hDMR</sub> and LenPro | Calculates per-region methylation difference and density metrics |
| `ICR_annotation.R` | Annotate and classify ICRs | Defines known and candidate imprinting genes based on Δ<sub>hDMR</sub> and LenPro |
| `hDMR_summary_plot.R` | Visualization of hDMR distribution | Plots hDMR density, Δ<sub>hDMR</sub> spectrum, and ICR overlaps |

> **Note:**  
> Users should configure file paths, haplotype identifiers, and reference annotation files (GENCODE/GENEIMPRINT) before execution.  
> Recommended to use consistent phasing output from **NanoMethPhase** to ensure accurate haplotype assignment.
