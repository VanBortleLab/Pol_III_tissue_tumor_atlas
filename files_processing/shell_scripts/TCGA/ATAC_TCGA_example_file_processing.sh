#!/bin/bash
# ----------------SLURM Parameters----------------
#SBATCH -p normal
#SBATCH -n 1
#SBATCH --mem=1200g
#SBATCH -N 1
#SBATCH -J TCGA-LIHC-124
#SBATCH -D /home/labs/kvbortle_lab/simonl3/ATAC/TCGA/sh_fil/
#SBATCH -e output.2.TCGA-LIHC-124.txt
#SBATCH -o error.2.TCGA-LIHC-124.txt

# ----------------Load Modules--------------------
module purge
module load BEDTools/2.28.0-IGB-gcc-8.2.0
# ----------------Commands------------------------
#
# TCGA files were processed from bam files. Refer to Metadata.
# ──────────────────────────────────────────────────────────────
# 1) Counts extraction
# ──────────────────────────────────────────────────────────────
bed_file_annotations="/home/personal/Tools/RNACentral_bins1.10.100_condense_pm150.bed"

align_dir="/home/personal/ATAC/align_dir"
counts_dir="/home/personal/ATAC/count_dir"

	 
bedtools coverage -a "${bed_file_annotations}" \
	-b "${align_dir}/8ff9ecd0-eb23-431a-b4cf-6eb0c19daa0f_atacseq_gdc_realn.bam" \
	> "${counts_dir}/TCGA-LIHC-124_RNAcentral.txt" -counts
