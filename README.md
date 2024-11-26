# ChinaMeth: Comprehensive DNA Methylation Analysis for Diverse Population

![Version](https://img.shields.io/badge/version-1.0.0-blue)
![Language](https://img.shields.io/badge/language-R-blue)
![License](https://img.shields.io/badge/license-MIT-green)
![Platform](https://img.shields.io/badge/platform-linux%20|%20macOS-brightgreen)

<div align="center">
    <img src="images/ChinaMeth.png" alt="ChinaMeth" width="300"/>
</div>

### A comprehensive DNA methylation atlas for the Chinese population through nanopore long-read sequencing of 106 individuals
Welcome to **ChinaMeth**, a project dedicated to exploring the complex patterns of DNA methylation across diverse Chinese populations using advanced long-read sequencing technologies. This repository provides insights into our comprehensive study of methylation dynamics and their role in genetic regulation and environmental adaptation. It also includes scripts for data processing and visulization to facilitate methylation studies.

<div align="center">
    <img src="images/chinameth.png" alt="chinameth" width="900"/>
</div>

## Key Features and Advantages
- **Population-Specific Insights**: By examining distinct populations (e.g., northern, southern, and xizang groups), ChinaMeth provides a deeper understanding of the methylation landscape and the potential influence of regional adaptations.
- **Advanced Technology**: Fully leveraging the advantages of long-read sequencing, ChinaMeth addresses limitations found in second-generation studies, providing high-resolution methylation data that uncovers patterns often overlooked by traditional methods.
- **Interactive Data Exploration**: ChinaMeth provides an online, interactive data visualization platform to facilitate research on methylation patterns across populations. This open-access resource allows researchers to explore and analyze methylation data seamlessly, supporting diverse epigenetic studies

## Novel Findings
ChinaMeth has unveiled several unique findings, such as:
- Elevated methylation levels in specific structural variations, with significant implications for regional adaptation.
- Identification of unique haplotype-based methylation patterns, especially among genes related to environmental resilience.
- Insights into hypo-methylation patterns distinctive to certain populations, suggesting underlying epigenetic regulation mechanisms influenced by local environments.

## Benefits for Future Researchers
- **Comprehensive ONT Analysis Pipeline**: ChinaMeth offers a complete ONT analysis workflow, including essential scripts for key processing steps, streamlining data analysis for new users and experienced researchers alike.
- **Open Data Access with Interactive Visualization**: The platform provides open access to methylation data with an interactive visualization interface, enabling researchers to explore and interpret methylation patterns across populations with ease.
- **Novel Findings for Reference**: ChinaMeth presents a range of unique insights and findings, providing valuable references and hypotheses for future research in epigenetics and population studies.

## Analysis Workflow
ChinaMeth provides a comprehensive workflow for DNA methylation data analysis. This workflow encompasses essential analysis scripts for processing methylation data, conducting SV-methylation correlation analyses, and generating visualizations.

### Methylation workflow
<div align="center">
    <img src="images/Meth_workflow.png" alt="Methylation workflow" width="500"/>
</div>

The Methylation analysis pipeline includes the following scripts:

1. **methylation_calling.sh**: Performs DNA basecalling, alignment, and methylation calling.
2. **methylation_phasing.sh**: Performs phasing and calculates haplotype-specific methylation frequencies.
3. **hDMR_calculate.sh**: Calculates and filters DMRs and DMCs.

### Additional Enrichment (Guppy + Nanopolish Traditional Workflow)

For users interested in the traditional Guppy + Nanopolish workflow, we provide the following scripts:

- **nanopy.sh**: Executes the traditional Guppy + Nanopolish workflow for basecalling and methylation calling.
- **calculate_methylation_frequency.py**: Computes methylation frequency based on Nanopolish results.

### SV workflow
<div align="center">
    <img src="images/SV_workflow.png" alt="SV worflow" width="400"/>
</div>

**SV Workflow** is a comprehensive workflow designed for analyzing methylation patterns in structural variants (SVs). The workflow consists of the following three main steps:
1. **SV Data Preprocessing**: Integrate individual data into a population-level dataset and perform data filtering and quality control. This step ensures data integrity and reliability for downstream analysis.
2. **SV Methylation Level Analysis**: Perform detailed analysis of methylation patterns in SV regions, with a focus on:
- **DEL Workflow**: Analyze methylation patterns in deletion (DEL) regions to explore their distribution and potential functional roles within the genome.
- **INS Workflow**: Investigate methylation patterns in insertion (INS) regions to uncover changes in methylation levels during the insertion process.
3. **Identifying Transpoable Elements**:Examine changes in methylation patterns of insertion (INS) regions as transposable elements (TEs). Using the **MEG Workflow**, this step explores the mechanisms behind methylation dynamics in TEs and their impact on genome stability.

### DEL workflow
<div align="center">
    <img src="images/del_pipeline.png" alt="DEL Pipeline" width="600"/>
</div>

The DEL analysis pipeline includes the following scripts:

1. **sv_sampleFilter.sh**: Filters and standardizes SV data for individual samples. 
2. **merge_pop.sh**: Merges SV data across populations, then filters and standardizes the merged data. 
3. **DEL_pop.sh**: Calculates sDMR (significant Differentially Methylated Region) methylation levels for DELs within populations. 
4. **DEL_plot.R**: Generates scatter and density plots for DEL methylation levels.

### INS workflow
<div align="center">
    <img src="images/ins_pipeline.png" alt="INS Pipeline" width="600"/>
</div>

The INS analysis pipeline includes the following scripts:

1. **extractReadFromINS.py**: Extracts methylation signals and sequences around INS (Insertion) variants. 
2. **compareSide2kbINS.sh**: Compares methylation levels between INS regions and their upstream/downstream 2kb regions. 
3. **ins_pop_merge.sh**: Merges individual methylation data files into a population-level file. 
4. **INS_plot.R**: Generates scatter and density plots for INS methylation levels.

## Feature
- Data preprocessing scripts
- Data analysis methods
- Visulization tools for results

## Directory Structure
- `data/`: Contains raw and processed data files
    - `DEL/`: Scripts specific to DEL analysis, including sample filtering, population merging, methylation level calculation, and plotting.
- `scripts/`: Contains scripts for data processing and analysis
- `plots/`: Contains generated plots and visualizations
- `docs/`: Documentation and Usage instructions

## Website
Explore CpG and three types of DMR distributions, including sDMR, hDMR, and pDMR, on our interctive [ChinaMeth](http://bioinformatics.hit.edu.cn/methylation).

## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for more details.

## Citation
If you use ChinaMethAtlas in your research, please cite the following paper: **A Comprehensive Tool for DNA Methylation Analysis in the Chinese population.**

## Contact
For any questions, please contact [email](yli21b@hit.edu.cn)
