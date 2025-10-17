### CpG Annotation and Statistics

This module focuses on the genome-wide characterization and comparative analysis of CpG methylation profiles across individuals and populations.  
It includes sequencing quality assessment, comparative evaluation with WGBS datasets, principal component analysis, and functional annotation of CpG sites.  
All relevant scripts are located in the `scripts/cpg_annotation/` directory, and detailed workflows are provided in [**Document/CpG/cpg_annotation.md**](Document/CpG/cpg_annotation.md).

---

#### Sequencing Quality Assessment

Raw signal files (**POD5**) were basecalled to FASTQ using **Dorado (v0.3.1)** with the `dna_r9.4.1_e8_hac` model.  
Sequencing metrics such as **read-length N50**, **mean Q-score**, and **predicted error rate** (calculated as *10^(-Q/10)*) were summarized using **NanoStat (v1.6.0)**.  
Per-sample coverage was computed with **Samtools (v1.9)** using the `depth` function, and cohort-level metrics were represented as arithmetic means of individual values.  

To assess data consistency, Pearson correlation coefficients were calculated between methylation levels at overlapping CpG sites across individuals.  
A **LOESS** (locally weighted scatterplot smoothing) regression was applied to model the relationship between sequencing coverage and the number of detected CpG sites per sample, revealing saturation trends in CpG detection with increasing coverage.

---

#### Comparative Analysis with WGBS Datasets

We compared our **ONT-derived methylation profiles** with publicly available **WGBS datasets** to evaluate CpG coverage and consistency.  
CpG sites located on **autosomes (chr1–22)** were retained for comparison.  
- **GSE181854:** two healthy Chinese males (NM1, NM2) and three healthy Chinese females (NF1–NF3).  
- **GSE186458:** twenty-two white blood cell (WBC) samples (GSM6810026–GSM6810048).  
- **GSE80911:** two stem cell samples (GSM2138820, GSM2138821), with CpGs filtered by coverage ≥ 3×.  

The mean value was used to summarize CpG counts across studies.  
Differences in CpG detection between ONT and WGBS were quantified using the following delta value:


<p align="center">
    <img src="https://latex.codecogs.com/svg.image?\Delta=\frac{CpG_{ONT}-CpG_{WGBS}}{CpG_{WGBS}}" title="\Delta=(CpG_{ONT}-CpG_{WGBS})/CpG_{WGBS}" />
</p>


where **CpG<sub>ONT</sub>** and **CpG<sub>WGBS</sub>** represent the number of CpGs detected by ONT and WGBS, respectively.  

---

#### Principal Component Analysis (PCA)

Shared CpG sites across individuals were merged using **Bedtools (v2.29)** with the `unionbedg` option, excluding sites with missing values (`NA`).  
CpGs exhibiting low population differentiation (≤ 0.05) were removed, where differentiation was defined as:

<p align="center">
    <img src="https://latex.codecogs.com/svg.image?diff(pop_1,pop_2)=\left|\frac{\sum\beta_1}{n}-\frac{\sum\beta_2}{m}\right|" title="diff(pop_1,pop_2)=|Σβ₁/n−Σβ₂/m|" />
</p>

Here, **n** and **m** denote the number of samples in populations 1 and 2, and **β** represents the methylation level of a CpG site.  

**PCAtools (v2.6.0)** was used to perform principal component analysis after excluding the top 10% of highly similar CpGs.  
Population differentiation was assessed using **t-tests** (95% confidence interval) based on the first four principal components.  
To further validate population clustering, **ANOSIM** (Analysis of Similarities) was conducted with the **Vegan (v2.6.4)** R package using the **Bray–Curtis** distance metric (`vegdist` function).

---
#### Functional Annotation

CpG methylation levels were annotated to genomic and regulatory features to investigate functional enrichment patterns.  
The overall methylation level (ML) of a functional element was calculated as:

<p align="center">
    <img src="https://latex.codecogs.com/png.image?ML(i,j)=\frac{\sum \beta}{\sum N_{CpG}}" alt="ML(i,j)=Σβ/ΣN_CpG" />
</p>

where **N<sub>CpG</sub>** is the number of CpGs within the region, and **β** represents the methylation level.
Methylation density (MD) was computed as:

<p align="center">
    <img src="https://latex.codecogs.com/png.image?MD(i,j)=\frac{\sum N_{CpG}}{j-i}" alt="MD(i,j)=ΣN_CpG/(j−i)" />
</p>

Annotation files were downloaded from [UCSC Genome Browser](http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database).  
For features available only in GRCh37, **liftOver** was used to convert coordinates to GRCh38.

**Annotation categories:**
- **Gene-related regions:** gene body, intergenic, exon, intron, CDS, UTR, TSS200, TSS1500, and promoter (−2 kb upstream of TSS).  
- **Histone marks:** DNase, H3K27ac, H3K27me3, H3K36me3, H3K4me1/2/3, H3K9ac, H3K9me3, and H4K20me1.  
- **Repeat elements:** SINE, LINE, LTR, satellites, retroposons, low-complexity regions (LCRs), and simple sequence repeats (SSRs).  
- **CpG islands:** CpG islands (CGIs), shores (±2 kb), shelves (±2 kb beyond shores), and open sea regions.  
- **Telomeric regions:** downstream of chromosome start and upstream of chromosome end (10 kb, 100 kb, and 1 Mb windows).

Statistical comparisons between two groups were conducted using the **Wilcoxon rank-sum test** (p < 0.05)

For multiple groups, the **Kruskal–Wallis test** was applied to evaluate differences among independent populations.

---

### Summary of Scripts

| Script Name | Purpose / Function | Notes |
|--------------|-------------------|-------|
| `cpg_merge.sh` | Merges CpG sites across individuals using Bedtools `unionbedg`. | Input: methylation files from multiple samples. |
| `pca_analysis.R` | Performs PCA and visualizes principal components using PCAtools. | Output: PCA plots and variance explained by components. |
| `anosim_test.R` | Conducts ANOSIM clustering validation using Vegan with Bray–Curtis distance. | Evaluates population-level methylation similarity. |
| `filter_lowdiff_cpg.py` | Filters CpG sites with low differentiation (≤ 0.05) between populations. | Ensures CpGs with stable methylation are excluded. |

---

#### Notes
- All parameters can be customized in the corresponding shell or R scripts under the `scripts/cpg_annotation/` directory.  
- Users should modify input file paths and reference genome versions according to their datasets before running the analysis.
