import heapq
import subprocess
import os
import pysam
import pandas as pd
import argparse
import alignMethSeq


class SequenceInfo:
    """
    Data structure to store sequence information and support heap comparison.
    """
    def __init__(self, sequence, length, name, meth):
        self.sequence = sequence
        self.length = length
        self.name = name
        self.meth = meth

    def __lt__(self, other):
        return self.length < other.length

cigar_operations = {
            0: "M", 1: "I", 2: "D", 3: "N", 4: "S", 5: "H", 6: "P", 7: "=", 8: "X"
            }


def run_abpoa(input_file, output_file):
    cmd = f"abpoa {input_file} > {output_file}"
    try:
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print("abpoa command executed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error running abpoa: {e.stderr.decode()}")


def run_abpoa_r1(input_file, output_file):
    cmd = f"abpoa {input_file} -r 1 > {output_file}"
    try:
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print("abpoa -r 1 command executed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error running abpoa -r1: {e.stderr.decode()}")


def find_cg_positions(sequence):
    return [i for i in range(len(sequence) - 1) if sequence[i:i + 2] == "CG"]


def construct_meth_level(meth_df, pos, value):
    level = 1 if value > 127 else 0
    if pos == 0:
        meth_df.append(level)
    else:
        while len(meth_df) < pos - 1:
            meth_df.append(2)
        meth_df.append(level)

def ins_step_one(in_bam, in_vcf, out_dir):
    shift_find_pos = 200
    shift_out_pos = 20
    os.makedirs(out_dir, exist_ok=True)

    out_fa = os.path.join(out_dir, os.path.basename(in_bam).replace(".bam", ".fa"))
    out_fm = os.path.join(out_dir, os.path.basename(in_bam).replace(".bam", ".fm"))

    if not (os.path.getsize(in_bam) and os.path.getsize(in_vcf)):
        print("Input BAM or VCF file is empty.")
        return

    tmp_vcf = pd.read_csv(in_vcf, header=None, sep='\t', names=["chr", "start", "end", "INFO"])
    ins_start = tmp_vcf.loc[0, "start"] - shift_find_pos
    var_num = int(tmp_vcf.loc[0, "INFO"].split(":")[2])

    bamfile = pysam.AlignmentFile(in_bam, "rb")
    top_sequences = []

    for read in bamfile:
        read_seq = read.query_alignment_sequence
        read_full = read.query_sequence
        read_name = read.query_name
        cigar_tuples = read.cigartuples
        alignment_start = read.reference_start
        ml_info = list(read.get_tag('ML')) if read.has_tag('ML') else []

        cg_positions = find_cg_positions(read_full)
        ins_pos_range = (ins_start - 1, ins_start + 2 * shift_find_pos)

        if not cigar_tuples:
            continue

        cigar_ref, cigar_read, clip_len, clip_flag = 0, 0, 0, 1
        ins_dict = {}

        for op, length in cigar_tuples:
            op_code = cigar_operations.get(op, '?')
            if op_code in {"S", "H"}:
                if clip_flag == 1:
                    clip_len = length
                continue
            elif op_code == "D":
                cigar_ref += length
            elif op_code == "I":
                index_pos = alignment_start + cigar_ref
                if ins_pos_range[0] < index_pos < ins_pos_range[1] and length >= 30:
                    ins_dict[cigar_read] = length
                cigar_read += length
            else:
                cigar_ref += length
                cigar_read += length
            clip_flag += 1


        if ins_dict:
            max_pos = max(ins_dict, key=ins_dict.get)
            max_len = ins_dict[max_pos]
            seq_start = max(0, max_pos - shift_out_pos)
            seq_end = min(len(read_seq), max_pos + max_len + shift_out_pos)
            sequence = read_seq[seq_start:seq_end]

            meth_profile = []
            cg_start = seq_start + clip_len
            cg_end = seq_end + clip_len
            for i, cg_index in enumerate(cg_positions):
                if cg_start <= cg_index <= cg_end:
                    cg_rel_start = cg_index - cg_start + 1
                    construct_meth_level(meth_profile, cg_rel_start, ml_info[i])
                elif cg_index > cg_end:
                    break
            meth_profile.extend([2] * (len(sequence) - len(meth_profile)))
            meth_str = "".join(map(str, meth_profile))

            heapq.heappush(top_sequences, SequenceInfo(sequence, max_len, read_name, meth_str))
            if len(top_sequences) > var_num:
                heapq.heappop(top_sequences)

    bamfile.close()

    with open(out_fa, "w") as f:
        for seq_info in top_sequences:
            f.write(f">{seq_info.name}\n{seq_info.sequence}\n")

    run_abpoa(out_fa, os.path.join(out_dir, os.path.basename(in_bam).replace(".bam", ".cons")))
    
    with open(out_fm, "w") as f:
        for seq_info in top_sequences:
            f.write(f">{seq_info.name}\n{seq_info.sequence}\n+\n{seq_info.meth}\n")

    fa_out = os.path.join(out_dir, os.path.basename(in_bam).replace(".bam", "_fa.fa"))
    run_abpoa_r1(out_fa, fa_out)

    alignMethSeq.insStepTwo(out_fm, fa_out)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Extract insertion sequence and methylation signal from BAM/VCF.")
    parser.add_argument('--bam', required=True, help='Input BAM file path')
    parser.add_argument('--vcf', required=True, help='Input pseudo-VCF (BED-like) file path')
    parser.add_argument('--outdir', required=True, help='Output directory path')
    args = parser.parse_args()

    ins_step_one(args.bam, args.vcf, args.outdir)
