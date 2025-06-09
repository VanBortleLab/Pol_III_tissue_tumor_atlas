#!/bin/bash
# ----------------SLURM Parameters----------------
#SBATCH -p normal
#SBATCH -n 1
#SBATCH --mem=200g
#SBATCH -N 1
#SBATCH -J ENCSR204SMO.heart_right_ventricle.1_2_4
#SBATCH -D /home/personal/ATAC/sh_files
#SBATCH -e output.8.ENCSR204SMO.heart_right_ventricle.1_2_4.txt
#SBATCH -o error.8.ENCSR204SMO.heart_right_ventricle.1_2_4.txt

# ----------------Load Modules--------------------
module load FastQC
module load Trim_Galore
module load Bowtie2
module load SAMtools
module load deepTools
module load BEDTools
# ----------------Commands------------------------

# ──────────────────────────────────────────────────────────────
# 0) Trimming and QC
# ──────────────────────────────────────────────────────────────


fastq_dir="/home/personal/ATAC/fastq_dir"
trim_dir="/home/personal/ATAC/trim_dir"

trim_galore --dont_gzip -o "${trim_dir}" \
	 "${out_dir4}/ENCFF520LZN.fastq.gz"  \
	--fastqc

# ──────────────────────────────────────────────────────────────
# 1) Alignment
# ──────────────────────────────────────────────────────────────	
index="/home/personal/bowtie_index/GRCh38_noalt_as/GRCh38_noalt_as"

align_dir="/home/personal/ATAC/align_dir"
	
bowtie2 -x "${index}" \
	 -U "${out_dir4}/ENCFF520LZN_trimmed.fq" \
	 -S "${align_dir}/ENCSR802GEV.cerebellum.1_1_1.sam"

	 
# ──────────────────────────────────────────────────────────────
# 2) Samtools processing
# ──────────────────────────────────────────────────────────────	 

# sam to bam	 
samtools view -S -b "${align_dir}/ENCSR802GEV.cerebellum.1_1_1.sam" \
	> "${align_dir}/ENCSR802GEV.cerebellum.1_1_1.bam"

# bam sorting
samtools sort "${align_dir}/ENCSR802GEV.cerebellum.1_1_1.bam" \
	-o "${align_dir}/ENCSR802GEV.cerebellum.1_1_1.sorted.bam"

# bam indexing
samtools index "${align_dir}/ENCSR802GEV.cerebellum.1_1_1.sorted.bam"

# ──────────────────────────────────────────────────────────────
# 3) Counts extraction
# ──────────────────────────────────────────────────────────────
bed_file_annotations="/home/personal/Tools/RNACentral_bins1.10.100_condense_pm150.bed"

counts_dir="/home/personal/ATAC/count_dir"
	 
bedtools coverage -a "${bed_file_annotations}" \
	-b "${align_dir}/ENCSR802GEV.cerebellum.1_1_1.sorted.bam" \
	> "${counts_dir}/ENCSR802GEV.cerebellum.1_1_1.bed" -counts

	 
