### Population-specific Differentially Methylated Region (pDMR) Analysis

This module focuses on the identification and interpretation of **population-specific differentially methylated regions (pDMRs)** and their associated **differentially methylated genes (DMGs)** across the North, South, and Xizang populations.  
It provides a framework for inter-population comparison, gene-level annotation, and altitude-associated methylation modeling.

---

### 1. pDMR Detection

Pairwise population comparisons were performed for **North vs. South**, **Xizang vs. North**, and **Xizang vs. South** using **Metilene (v0.2.8)**.  
Each analysis included both **random subsampling (10 individuals per population)** and a **full-cohort analysis** using all available samples.  
The complete set of DMRs was generated and filtered using the official **Metilene** scripts (`metilene.input.pl` and `metilene.output.pl`), retaining significant regions based on q-value and methylation difference thresholds.  
Detected pDMRs were categorized as **hyper** (mean methylation difference > 0) or **hypo** (mean methylation difference < 0), and their genomic distributions were visualized using **TBtools (v1.115)** in circos format.

---

### 2. Differentially Methylated Gene (DMG) Identification

Differentially methylated genes were defined by integrating pDMRs and differential methylated CpG sites overlapping **gene bodies** or **±2 kb regulatory regions**.  
Genes were retained if DMRs contained ≥ 20 CpGs and covered more than 50% of the gene region.  
Population-level comparisons were performed to identify genes showing distinct methylation profiles among high- and low-altitude groups.  
Clustering and visualization were implemented using the **ComplexHeatmap (v2.16.0)** package in R.

---

### 3. Altitude-associated Methylation Modeling

Genes were classified into **high-altitude adaptation (HA)** and **novel high-altitude adaptation (nHA)** categories based on prior literature validation.  
Statistical evaluations of methylation differences between populations employed **Wilcoxon rank-sum tests**, **Cliff’s Delta**, and **mean difference** metrics.  
A **logistic regression model** was used to fit methylation–altitude relationships:

<p align="center">
  y = y<sub>min</sub> + (y<sub>max</sub> − y<sub>min</sub>) / (1 + e<sup>−k(x − m₀)</sup>)
</p>

where *m₀* represents the methylation inflection point, and *k* reflects the rate of methylation change with altitude.  
Model fitting significance was assessed using the **F-statistic p-value**, and R² values were used to evaluate fit performance.

---

### 4. Enrichment Analysis and Discovery of Novel Altitude-adaptive Genes

To interpret the functional implications of DMGs, we performed **GO** and **KEGG** enrichment analyses using **Metascape (v3.5.20230501)**, applying the following thresholds:  
minimum intersection ≥ 3, p-value ≤ 0.01, and enrichment factor ≥ 1.5.  
The top 20 clusters were used to construct an **enrichment network**, where terms with similarity > 0.3 were connected by edges.  
Outlier clusters containing fewer terms were removed using **Cytoscape (v3.10.1)**.

For identifying **novel high-altitude adaptation (HA) genes**, we applied a multi-step filtering strategy:  
- Removed nodes with degree < 15 and −log(p-value) < 15 in the enrichment network.  
- Selected one representative node per cluster based on the highest degree.  
- Retained the **top 5 nodes** with the highest degrees in the largest subnetwork, prioritizing higher −log(p-value) in case of ties.  

Finally, differences in methylation profiles for representative genes were visualized using the **NanoMethViz (v2.6.0)** R package (`plot_gene` function).

---


#### Summary of Scripts

| Script Name | Function | Description |
|--------------|-----------|-------------|
| `prepare_input.sh` | Generate input files | Merges population BED files for Metilene input |
| `pDMR_metilene.sh` | Detect population-specific DMRs | Runs Metilene to identify DMRs between population pairs |
| `filter_pDMR.sh` | Filter significant DMRs | Applies thresholds for q-value, α, and minimum DMR length |
| `DMG_annotation.py` | Identify DMGs | Maps filtered pDMRs to gene regions |
| `altitude_model.R` | Fit logistic models | Evaluates altitude–methylation associations |
| `pDMR_visualization.R` | Plot results | Generates circos and heatmap visualizations |

> **Note:**
> Ensure consistent reference genome (GRCh38) and coordinate alignment across datasets before execution.  
> Recommended to maintain balanced sample sizes for each population when performing subsampling analyses.
