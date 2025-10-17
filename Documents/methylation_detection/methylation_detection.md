### Methylation workflow
<div align="center">
    <img src="../../images/Meth_workflow.jpg" alt="Methylation workflow" width="500"/>
</div>

The Methylation analysis pipeline includes the following scripts:

1. **methylation_calling.sh**: Performs DNA basecalling, alignment, and methylation calling.
2. **methylation_phasing.sh**: Performs phasing and calculates haplotype-specific methylation frequencies.
3. **hDMR_calculate.sh**: Calculates and filters DMRs and DMCs.

### Additional Enrichment

For users interested in the traditional Guppy + Nanopolish workflow, we provide the following scripts:

- **nanopy.sh**: Executes the traditional Guppy + Nanopolish workflow for basecalling and methylation calling.
- **calculate_methylation_frequency.py**: Computes methylation frequency based on Nanopolish results.
