"""
Script Purpose:
    This script is designed to identify the source location of INS (insertion) consensus sequences
    based on methylation data. It processes insertion methylation levels, calculates methylation
    differences across genomic regions (upstream, downstream, and body), and aligns consensus sequences
    to a reference genome to determine their origin.

    Key Features:
    1. Identify INS regions and calculate their methylation levels.
    2. Align consensus sequences using `minimap2` to find their source locations in the reference genome.
    3. Extract methylation levels for source locations and surrounding regions.
    4. Classify methylation levels into categories such as "High", "Low", and "Other" based on differences.
"""

import pandas as pd
import subprocess
import pysam
import os

# Global variables (abstracted paths)
REFERENCE_GENOME = "<path_to_reference_genome>"  # Reference genome path
GENOME_POS = "<path_to_genome_pos>"  # Genome position file path
THREADS = 8  # Number of threads for minimap2
MIN_METH_DIFF = 0.1  # Minimum methylation difference threshold
CHROMOSOMES = [f"chr{i}" for i in range(1, 23)]  # Autosomes
STEP_TWO_SCRIPT = "./reAlignTwo.sh"  # Path to the step two shell script
ROOT_DIR = "<root_directory>"  # Root directory for input and output files

def run_minimap2(cons_file, out_file):
    """
    Aligns the consensus sequence to the reference genome using minimap2.

    Args:
        cons_file (str): Path to the consensus sequence file.
        out_file (str): Path to the output SAM file.
    """
    cmd = f"minimap2 -t {THREADS} -ax map-ont {REFERENCE_GENOME} {cons_file} > {out_file}"
    try:
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print("Minimap2 executed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error running minimap2: {e.stderr.decode()}")

def run_step_two(chr_name, start, end, map_file):
    """
    Calls a shell script to extract methylation levels for a given region.

    Args:
        chr_name (str): Chromosome name.
        start (int): Start position of the region.
        end (int): End position of the region.
        map_file (str): Path to the input BED file.

    Returns:
        float: Methylation level as a percentage, or -1 if an error occurs.
    """
    cmd = f"bash {STEP_TWO_SCRIPT} {chr_name} {start} {end} {map_file}"
    try:
        result = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        output = result.stdout.strip()
        return float(output) / 100 if output.isdigit() else -1
    except subprocess.CalledProcessError as e:
        print(f"Error running step two script: {e.stderr}")
        return -1

def is_within_region(chr_new, pos, chr_raw, start, end):
    """
    Checks if a position is within a specific region on the same chromosome.

    Args:
        chr_new (str): Chromosome name of the new position.
        pos (int): Position to check.
        chr_raw (str): Chromosome name of the reference region.
        start (int): Start of the reference region.
        end (int): End of the reference region.

    Returns:
        bool: True if the position is within the region, False otherwise.
    """
    return chr_new == chr_raw and start <= pos <= end

def determine_methylation_type(ins_meth, up_meth, down_meth):
    """
    Classifies methylation levels into categories based on differences.

    Args:
        ins_meth (float): Methylation level of the body region.
        up_meth (float): Methylation level of the upstream region.
        down_meth (float): Methylation level of the downstream region.

    Returns:
        str: Methylation category ("High", "Low", "Other", etc.).
    """
    if ins_meth + up_meth + down_meth < 0:
        return "NA"
    ins_up_diff = ins_meth - up_meth
    ins_down_diff = ins_meth - down_meth
    total_diff = ins_up_diff ** 2 + ins_down_diff ** 2

    if ins_up_diff > 0 and ins_down_diff > 0:
        return "High_delta" if total_diff >= MIN_METH_DIFF ** 2 else "High"
    elif ins_up_diff < 0 and ins_down_diff < 0:
        return "Low_delta" if total_diff >= MIN_METH_DIFF ** 2 else "Low"
    else:
        return "Other_delta" if total_diff >= MIN_METH_DIFF ** 2 else "Other"

def calculate_methylation_levels(map_chr, map_start, map_end, sample):
    """
    Calculates methylation levels for a given region (body, upstream, downstream).

    Args:
        map_chr (str): Chromosome name.
        map_start (int): Start position of the region.
        map_end (int): End position of the region.
        sample (str): Sample name.

    Returns:
    tuple: Methylation levels (body, upstream, downstream).
    """
    if map_chr not in CHROMOSOMES:
        return -1, -1, -1

    bed_file = os.path.join("<bed_directory>", sample, "imprint", f"{sample}_noIm.bed")
    ins_meth = run_step_two(map_chr, map_start, map_end, bed_file)
    up_meth = run_step_two(map_chr, map_start - 2000, map_start - 1, bed_file)
    down_meth = run_step_two(map_chr, map_end + 1, map_end + 2000, bed_file)

    return ins_meth, up_meth, down_meth

def realign(sample, pop):
    """
    Main function to find the source location of INS consensus sequences and
    calculate methylation levels for the body and surrounding regions.
    
    Args:
        sample (str): Sample name.
        pop (str): Population group name.
    """
    sample_dir = os.path.join(ROOT_DIR, pop, sample)
    out_dir = os.path.join(sample_dir, "reAlign")
    os.makedirs(out_dir, exist_ok=True)
    
    input_cpg = os.path.join(sample_dir, "cpg5.cpg")
    input_cons = os.path.join(sample_dir, "result")
    result_file = os.path.join(out_dir, f"{sample}_reAlign.result")
    sam_dir = os.path.join(out_dir, "sam")
    os.makedirs(sam_dir, exist_ok=True)
    
    with open(result_file, "w") as f:
        f.write("rawChr\trawPos\tinsMeth\tupMeth\tdownMeth\tType\tmapChr\tmapStart\tmapEnd\tmapMeth\tupMap\tdownMap\tmapLevel\tquality\tsample\n")
    
    if os.path.getsize(input_cpg):
        cpg_data = pd.read_csv(input_cpg, sep='\t', header=None, names=["chr", "start", "end", "ins", "up", "down", "sample", "upNum", "downNum"])
    
        for _, row in cpg_data.iterrows():
            chromosome, ins_pos, ins_meth, up_meth, down_meth = row["chr"], row["start"], row["ins"], row["up"], row["down"]
            level_type = determine_methylation_type(ins_meth, up_meth, down_meth)
    
            ins_id = f"{chromosome}_{ins_pos}"
            tmp_dir = os.path.join(input_cons, ins_id)
            cons_file = os.path.join(tmp_dir, f"{sample}_{ins_id}.cons")
            sam_file = os.path.join(sam_dir, f"{sample}_{ins_id}.sam")
    
        if os.path.exists(cons_file) and os.path.getsize(cons_file):
            if not os.path.exists(sam_file):
                run_minimap2(cons_file, sam_file)
    
        if os.path.getsize(sam_file):
            samfile = pysam.AlignmentFile(sam_file, "rb")
            best_align = next(samfile)
            map_chr, map_pos, map_len, map_qual = best_align.reference_name, best_align.reference_start, best_align.query_alignment_length, best_align.mapping_quality
    
            raw_start, raw_end = max(0, ins_pos - 10000), ins_pos + 10000
            if map_qual >= 20 and not is_within_region(map_chr, map_pos, chromosome, raw_start, raw_end):
                map_end = map_pos + map_len
                map_body_2kb = calculate_methylation_levels(map_chr, map_pos, map_end, sample)
                map_level = determine_methylation_type(map_body_2kb[0], map_body_2kb[1], map_body_2kb[2])

            with open(result_file, "a") as f:
                f.write(f"{chromosome}\t{ins_pos}\t{ins_meth}\t{up_meth:.4f}\t{down_meth:.4f}\t{level_type}\t")
                f.write(f"{map_chr}\t{map_pos}\t{map_end}\t{map_body_2kb[0]:.4f}\t{map_body_2kb[1]:.4f}\t{map_body_2kb[2]:.4f}\t")
                f.write(f"{map_level}\t{map_qual}\t{sample}\n")
            samfile.close()
        else:
            print(f"Consensus file {cons_file} does not exist.")
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Manual to reAlign")
    parser.add_argument("--sample", type=str, default="sam1")
    parser.add_argument("--population", type=str, default="pop")
    args = parser.parse_args()
    sample = args.sample
    pop = args.population
    realign(sample,pop)       
