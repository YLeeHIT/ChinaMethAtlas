"""
Script Purpose:
This script is designed to merge individual INS analysis results into a single population-level result.
It processes individual result files, calculates various methylation statistics, and outputs 
merged population data for further analysis.

Main Features:
1. Merges individual result files into a single, sorted population-level file.
2. Annotates INS source locations and calculates overlap regions.
3. Computes population-level statistics, including methylation levels and INS classifications.
4. Filters INS data based on support thresholds to retain high-confidence results.

Workflow:
- Step 1: Merge individual files into a single population-level file, sorted by genomic position.
- Step 2: Annotate INS source regions using external shell scripts.
- Step 3: Calculate methylation statistics (mean, overlap center, etc.) at the population level.
- Step 4: Classify INS and mapping types based on methylation patterns.
- Step 5: Filter results based on support thresholds and save final population-level results.

Usage:
- Update input paths and parameters as required.
- Run the `run_popSta(pop)` function with the population name as an argument.

Dependencies:
- Python libraries: pandas, numpy, subprocess, os
- External tools: Shell scripts for annotation, data processing.

"""

import pandas as pd
import numpy as np
import subprocess
import os

# Step 1: Merge individual files into a single population file
def merge_file(inDir, pop):
    """
    Merges individual INS result files into a single population-level file.

    Args:
        inDir (str): Input directory containing individual result files.
        pop (str): Population name.

    Returns:
        str: Path to the merged output file.
    """
    outDir = os.path.join(inDir, "out")
    outfile = os.path.join(outDir, f"{pop}_merge.result")
    os.makedirs(outDir, exist_ok=True)

    cmd = f"cat {inDir}/*result | grep -v '^rawChr' | sort -k1,1V -k2,2n -k7,7V > {outfile}"
    try:
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print("File merge executed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred during file merge: {e.stderr.decode()}")
    return outfile


def bash_anno(pop, outdir):
    """
    Runs a shell script to annotate the INS source regions and perform statistical analysis.

    Args:
        pop (str): Population name.
        outdir (str): Directory containing the input and output files for annotation.
    """
    annoShell = "<path_to_annotation_script>"
    cmd = f"bash {annoShell} {pop} {outdir}"
    try:
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print("Annotation shell script executed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred during annotation: {e.stderr.decode()}")


def safe_nanmean(x):
    """
    Safely calculates the mean of a list while ignoring NaN values.

    Args:
        x (list): A list of numerical values.

    Returns:
        float: Mean of the list, or NaN if the list is empty or contains only NaN values.
    """
    return np.nan if len(x) == 0 or np.isnan(x).all() else np.nanmean(x)


def calculate_sampleNum(df):
    """
    Calculates the unique number of samples in the 'sample' column.

    Args:
        df (pd.DataFrame): Input DataFrame.

    Returns:
        float: Half of the unique sample count, or 0 if invalid.
    """
    sampleNum = df['sample'].nunique()
    try:
        return sampleNum / 2 if isinstance(sampleNum, (int, float)) else 0
    except Exception as e:
        print(f"Error calculating sample number: {e}")
        return 0


def find_max_overlap_center(starts, ends):
    """
    Finds the region with the maximum overlap and calculates its center.

    Args:
        starts (list): List of start positions.
        ends (list): List of end positions.

    Returns:
        tuple: Center of the maximum overlap region, overlap count, and overlap length.
    """
    intervals = sorted(zip(starts, ends))
    overlap_intervals = []
    current_start, current_end = intervals[0]
    overlap_count = 1

    for i in range(1, len(intervals)):
        if intervals[i][0] <= current_end:
            current_start = max(current_start, intervals[i][0])
            current_end = min(current_end, intervals[i][1])
            overlap_count += 1
        else:
            overlap_intervals.append((current_start, current_end, overlap_count))
            current_start, current_end = intervals[i]
            overlap_count = 1

    overlap_intervals.append((current_start, current_end, overlap_count))
    max_overlap = max(overlap_intervals, key=lambda x: x[2])
    overlap_center = int((max_overlap[0] + max_overlap[1]) / 2)
    overlap_length = max_overlap[1] - max_overlap[0]

    return overlap_center, max_overlap[2], overlap_length


def determine_map_type(row):
    """
    Determines the mapping type (Low, High, Other) based on methylation levels.

    Args:
        row (pd.Series): A row of the DataFrame.

    Returns:
        str: Mapping type.
    """
    if row['mapMethMean'] < min(row['mapUpMean'], row['mapDownMean']):
        return 'mapLow'
    elif row['mapMethMean'] > max(row['mapUpMean'], row['mapDownMean']):
        return 'mapHigh'
    else:
        return 'mapOther'


def determine_ins_type(row):
    """
    Determines the INS type (Low, High, Other) based on methylation levels.

    Args:
        row (pd.Series): A row of the DataFrame.

    Returns:
        str: INS type.
    """
    if row['insMeth'] < min(row['insUp'], row['insDown']):
        return 'insLow'
    elif row['insMeth'] > max(row['insUp'], row['insDown']):
        return 'insHigh'
    else:
        return 'insOther'


def statistic_result(infile):
    """
    Generates statistics for transType, mapType, and insType.

    Args:
        infile (pd.DataFrame): Input DataFrame.

    Returns:
        tuple: Statistics for transType, mapType, and insType.
    """
    counts_trans = infile['transType'].value_counts()
    counts_map = infile['mapType'].value_counts()
    counts_ins = infile['insType'].value_counts()

    trans_stats = f"{counts_trans.get('Low', 0)}:{counts_trans.get('High', 0)}"
    map_stats = f"{counts_map.get('mapLow', 0)}:{counts_map.get('mapHigh', 0)}:{counts_map.get('mapOther', 0)}"
    ins_stats = f"{counts_ins.get('insLow', 0)}:{counts_ins.get('insHigh', 0)}:{counts_ins.get('insOther', 0)}"

    return trans_stats, map_stats, ins_stats


def write_sta_result(statistics, pop, outDir, popMerge, resultFile):
    """
    Writes the statistics to an output file.

    Args:
        statistics (tuple): Statistics for transType, mapType, and insType.
        pop (str): Population name.
        outDir (str): Output directory.
        popMerge (pd.DataFrame): Merged DataFrame.
        resultFile (pd.DataFrame): Filtered DataFrame.
    """
    outStaFile = f"{pop}_result.sta"
    outStaPath = os.path.join(outDir, outStaFile)

    with open(outStaPath, "w") as f:
        f.write(f"Population:\t{pop}\n")
        f.write(f"Raw_reFile:\t{popMerge.shape[0]}\n")
        f.write(f"Filter_file:\t{resultFile.shape[0]}\n")
        f.write(f"Statistics:\t{statistics}\n")


def run_popSta(pop):
    """
    Main function to process population-level data, perform merging, filtering, and annotation.

    Args:
        pop (str): Population name.
    """
    rootDir = "<path_to_population_directory>"
    popDir = os.path.join(rootDir, pop)
    outDir = os.path.join(popDir, "out")

    infile = merge_file(popDir, pop)
    if os.path.exists(infile):
        popFile = pd.read_csv(infile, sep='\t', names=["insChr", "insPos", "insMeth", "insUp", "insDown", "insType",
                                                       "mapChr", "mapStart", "mapEnd", "mapMeth", "mapUp", "mapDown",
                                                       "mapType", "quality", "sample"])

        # Clean and process data
        popFile.replace(-1, np.nan, inplace=True)
        popTheshold = calculate_sampleNum(popFile)

        popMerge = popFile.groupby(['insChr', 'insPos', 'mapChr'], as_index=False).agg({
            'sample': list,
            'insMeth': 'mean',
            'insUp': 'mean',
            'insDown': 'mean',
            'mapMeth': list,
            'mapUp': list,
            'mapDown': list,
            'mapStart': list,
            'mapEnd': list,
            'quality': list
        })

        popMerge['Support'] = popMerge['sample'].apply(len)
        popMerge['mapMethMean'] = popMerge['mapMeth'].apply(safe_nanmean)
        popMerge['mapUpMean'] = popMerge['mapUp'].apply(safe_nanmean)
        popMerge['mapDownMean'] = popMerge['mapDown'].apply(safe_nanmean)

        popMerge[['overlapCenter', 'overlapNum', 'overlapLength']] = popMerge.apply(
            lambda row: pd.Series(find_max_overlap_center(row['mapStart'], row['mapEnd'])), axis=1)

        popMerge['transType'] = popMerge.apply(
            lambda row: 'High' if row['insMeth'] > row['mapMethMean'] else 'Low', axis=1)
        popMerge['mapType'] = popMerge.apply(determine_map_type, axis=1)
        popMerge['insType'] = popMerge.apply(determine_ins_type, axis=1)

        resultFile = popMerge[popMerge['Support'] >= popTheshold]
        resultFile.to_csv(os.path.join(outDir, f"{pop}_filtered.result"), sep='\t', index=False)
