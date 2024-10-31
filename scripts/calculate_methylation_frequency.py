#!/usr/bin/env python
# coding=utf-8
"""
Nanopolish Methylation Frequency Calculator

This script calculates the methylation frequency at genomic CpG sites from the output
of `nanopolish call-methylation`. It computes the methylation percentage based on 
the log likelihood ratios (LLRs) of methylation calls across reads at each CpG site 
or group of CpG sites, optionally allowing for groups to be split into individual 
CpG sites. 

Original Source:
    This code was inspired by the `calculate_methylation_frequency.py` script from the 
    [nanopolish GitHub repository](https://github.com/jts/nanopolish).

    Usage:
    ------
    Run the script with input TSV files containing methylation calls. 
    Use the `-c` option to specify a call threshold for LLRs, and `-s` to split groups.

    Example:
        python calculate_methylation_frequency.py -c 2.0 -s sample_methylation_calls.tsv > methylation_frequency.tsv

        Arguments:
        ----------
        - `-c` or `--call-threshold` (float): Specifies the LLR threshold for considering a site methylated.
          Default is 2.0.
          - `-s` or `--split-groups` (flag): If set, splits multi-CpG groups into individual CpG sites.

          Output:
          -------
          The script prints a tab-separated table with the following columns:
          - chromosome: Chromosome of the CpG site.
          - start: Start position of the CpG site.
          - end: End position of the CpG site.
          - num_motifs_in_group: Number of CpG motifs in the group.
          - called_sites: Number of called CpG sites.
          - called_sites_methylated: Number of methylated CpG sites.
          - methylated_frequency: Methylation frequency (methylated/called).
          - group_sequence: The sequence of the CpG group.

"""
import sys
import csv
import argparse
import gzip

class SiteStats:
    def __init__(self, group_size, sequence):
        self.num_reads = 0
        self.called_sites = 0
        self.called_sites_methylated = 0
        self.group_size = group_size
        self.sequence = sequence

def update_call_stats(key, num_called_sites, is_methylated, sequence):
    """Updates methylation statistics for each site or site group."""
    if key not in sites:
        sites[key] = SiteStats(num_called_sites, sequence)
    sites[key].num_reads += 1
    sites[key].called_sites += num_called_sites
    if is_methylated:
        sites[key].called_sites_methylated += num_called_sites

# Set up argument parser
parser = argparse.ArgumentParser(description="Calculate methylation frequency at genomic CpG sites")
parser.add_argument('-c', '--call-threshold', type=float, default=2.0, help="Log likelihood ratio threshold for methylation call")
parser.add_argument('-s', '--split-groups', action='store_true', help="Split multi-CpG groups into individual sites")
args, input_files = parser.parse_known_args()

sites = {}

# Process input files and gather per-site statistics
for f in input_files:
    with (gzip.open(f, 'rt') if f.endswith(".gz") else open(f)) as in_fh:
        csv_reader = csv.DictReader(in_fh, delimiter='\t')
        for record in csv_reader:
            num_sites = int(record['num_motifs'])
            llr = float(record['log_lik_ratio'])

            # Skip ambiguous calls below threshold
            if abs(llr) < args.call_threshold * num_sites:
                continue

            is_methylated = llr > 0
            sequence = record['sequence']
            chrom = record['chromosome']
            start = int(record['start'])
            end = int(record['end'])

            # Split multi-CpG groups if required
            if args.split_groups and num_sites > 1:
                cg_pos = sequence.find("CG")
                first_cg_pos = cg_pos
                while cg_pos != -1:
                    key = (chrom, start + cg_pos - first_cg_pos, start + cg_pos - first_cg_pos)
                    update_call_stats(key, 1, is_methylated, "split-group")
                    cg_pos = sequence.find("CG", cg_pos + 1)
            else:
                key = (chrom, start, end)
                update_call_stats(key, num_sites, is_methylated, sequence)

# Output results
print("\t".join(["chromosome", "start", "end", "num_motifs_in_group", "called_sites", 
                 "called_sites_methylated", "methylated_frequency", "group_sequence"]))

for key in sorted(sites.keys(), key=lambda x: (x[0], x[1], x[2])):
    stats = sites[key]
    if stats.called_sites > 0:
        frequency = stats.called_sites_methylated / stats.called_sites
        print(f"{key[0]}\t{key[1]}\t{key[2]}\t{stats.group_size}\t{stats.called_sites}\t"
              f"{stats.called_sites_methylated}\t{frequency:.3f}\t{stats.sequence}")

