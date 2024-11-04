# ChinaMethAtlas

### A comprehensive DNA methylation atlas for the Chinese population through nanopore long-read sequencing of 106 individuals
ChinaMethAtlas is a project focused on the analysis of DNA methylation data for the Chinese population. The project aims to provide tools and methods for the preprocessing, analysis, and visualization of methylation data. This repository contains scripts for data processing and visulization to facilitate methylation studies.

## Analysis Workflow
### Methylation workflow
![Methylation Workflow](images/Meth_workflow.png)
### SV workflow
![SV Workflow](images/SV_workflow.png)

The DEL analysis pipeline includes the following scripts:
1. **sv_sampleFilter.sh**: Filters and standardizes SV data for individual samples. 
2. **merge_pop.sh**: Merges SV data across populations, then filters and standardizes the merged data. 
3. **DEL_pop.sh**: Calculates sDMR (significant Differentially Methylated Region) methylation levels for DELs within populations. 
4. **DEL_plot.R**: Generates scatter and density plots for DEL methylation levels.

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

## Installation and Usage

### Clone the Repository

```bash
git clone https://github.com/YLeeHIT/ChinaMethAtlas.git
cd ChinaMethAtlas
```

### Usage Example
1. Data Preprocessing:
    ```bash
    python scripts/data_preprocessing.py
    ```

2. Data Analysis:
    ```bash
    bash scripts/data_analysis.sh
    ```

3. Results Visulization:
    ```bash
    Rscript scripts/plot_results.R
    ```

## Website
Explore CpG and three types of DMR distributions, including sDMR, hDMR, and pDMR, on our interctive [ChinaMeth](http://bioinformatics.hit.edu.cn/methylation).

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for more details.

## Citation

If you use ChinaMethAtlas in your research, please cite the following paper: **A Comprehensive Tool for DNA Methylation Analysis in the Chinese population.**

## Contact

For any questions, please contact [email](yli21b@hit.edu.cn)
