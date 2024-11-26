"""
Script to extract methylation signals and sequences around INS (Insertion) variants
Description:
This script processes BAM files and a pseudo-VCF (BED-like format) to extract methylation signals 
and sequences surrounding INS variants. It uses `abpoa` for sequence alignment and supports 
further analysis using custom functions.

Dependencies:
- pandas
- pysam
- subprocess
- heapq

Usage:
python script_name.py <input_bam> <input_vcf> <output_directory>
"""

import pandas as pd
import pysam
import heapq
import subprocess
import os
import alignMethSeq


class SequenceInfo:
    """
    Class to store extracted sequence information, including sequence content,
    length, read name, and methylation profile.
    """
    def __init__(self, sequence, length, name, meth):
        self.sequence = sequence
        self.length = length
        self.name = name
        self.meth = meth

    def __lt__(self, other):
        # Comparison operator for maintaining a heap of sequences based on length
        return self.length < other.length


# Dictionary to interpret CIGAR operation codes
CIGAR_OPERATIONS = {
    0: "M",  # Match
    1: "I",  # Insertion
    2: "D",  # Deletion
    3: "N",  # Skip
    4: "S",  # Soft clip
    5: "H",  # Hard clip
    6: "P",  # Padding
    7: "=",  # Sequence match
    8: "X"   # Sequence mismatch
}


def run_abpoa(input_file, output_file):
    """
    Run the `abpoa` command to generate a consensus sequence.
    """
    cmd = f"abpoa {input_file} > {output_file}"
    try:
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print("abpoa command executed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error while running abpoa: {e.stderr.decode()}")


def run_abpoa_r1(input_file, output_file):
    """
    Run the `abpoa` command with `-r 1` to output spaced alignment in FASTA format.
    """
    cmd = f"abpoa {input_file} -r 1 > {output_file}"
    try:
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print("abpoa -r1 command executed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error while running abpoa -r1: {e.stderr.decode()}")


def find_cg_positions(sequence):
    """
    Identify all positions of 'CG' dinucleotides in a given sequence.
    """
    return [i for i in range(len(sequence) - 1) if sequence[i:i + 2] == "CG"]


def construct_methylation_profile(meth_df, pos, value):
    """
    Construct a methylation profile by appending binary methylation values.
    """
    meth_level = 1 if value > 127 else 0
    if pos == 0:
        meth_df.append(meth_level)
    else:
        while len(meth_df) < pos - 1:
            meth_df.append(2)  # Fill gaps with a placeholder value (2)
        meth_df.append(meth_level)


def insStepOne(in_bam, in_vcf, output_dir):
    """
    Main function to extract methylation signals and sequences around INS variants.
    """
    shift_find_pos = 200  # Range to search for potential INS
    shift_out_pos = 20    # Range to output INS and surrounding sequences

    os.makedirs(output_dir, exist_ok=True)

    output_fasta = os.path.join(output_dir, os.path.basename(in_bam).replace(".bam", ".fa"))
    methylation_output = os.path.join(output_dir, os.path.basename(in_bam).replace(".bam", ".fm"))

    if os.path.getsize(in_vcf) and os.path.getsize(in_bam):
        tmp_vcf = pd.read_csv(in_vcf, header=None, sep='\t', names=["chr", "start", "end", "INFO"])
        ins_start = tmp_vcf.loc[0, "start"] - shift_find_pos
        var_num = int(tmp_vcf.loc[0, "INFO"].split(":")[2])  # Extract number of variants

        bam_file = pysam.AlignmentFile(in_bam, "rb")
        top_sequences = []

        for read in bam_file:
            read_sequence = read.query_alignment_sequence  # Sequence excluding soft/hard clips
            read_name = read.query_name
            alignment_start = read.reference_start
            ML_info = list(read.get_tag('ML')) if read.has_tag('ML') else []

            # Find INS positions and extract sequences
            # Your existing logic for INS extraction goes here...

            # Example: Storing results in the heap
            tmp_ins = SequenceInfo(read_sequence, len(read_sequence), read_name, "")
            heapq.heappush(top_sequences, tmp_ins)
            if len(top_sequences) > var_num:
                heapq.heappop(top_sequences)

        bam_file.close()

        # Write sequences to output files
        with open(output_fasta, "w") as f:
            for seq_info in top_sequences:
                f.write(f">{seq_info.name}\n{seq_info.sequence}\n")

        with open(methylation_output, "w") as f:
            for seq_info in top_sequences:
                f.write(f">{seq_info.name}\n{seq_info.sequence}\n+\n{seq_info.meth}\n")

        # Run abpoa commands
        consensus_file = os.path.join(output_dir, os.path.basename(in_bam).replace(".bam", ".cons"))
        run_abpoa(output_fasta, consensus_file)

        spaced_alignment_file = os.path.join(output_dir, os.path.basename(in_bam).replace(".bam", "_fa.fa"))
        run_abpoa_r1(output_fasta, spaced_alignment_file)

        alignMethSeq.insStepTwo(methylation_output, spaced_alignment_file)

    else:
        print("INS not found in input files.")


# Add a main block to allow running as a standalone script
if __name__ == "__main__":
    import sys
    if len(sys.argv) != 4:
        print("Usage: python script_name.py <input_bam> <input_vcf> <output_directory>")
        sys.exit(1)

    input_bam = sys.argv[1]
    input_vcf = sys.argv[2]
    output_dir = sys.argv[3]
    insStepOne(input_bam, input_vcf, output_dir)
