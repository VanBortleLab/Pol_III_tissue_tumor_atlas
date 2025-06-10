#!/bin/bash
# ----------------SLURM Parameters----------------
#SBATCH -p normal
#SBATCH -n 1
#SBATCH --mem=40g
#SBATCH -N 1
#SBATCH --mail-user=simonl3@illinois.edu
#SBATCH --mail-type=FAIL
# ----------------Load Modules--------------------
module purge
module load SAMtools
module load deepTools

### Example for generating bw files for Tissue - Heart

# ──────────────────────────────────────────────────────────────
# PART I - For every unique experiment
# ──────────────────────────────────────────────────────────────

# ──────────────────────────────────────────────────────────────
# 1) Changing names to the original bam file
# ──────────────────────────────────────────────────────────────

align_dir='/home/personal/align_dir' ## All the bam files for each experiment and its replicates is here

bam_dir='/home/personal/bam_dir' ## Output for merged bam file per experiment

cp  "${align_dir}/99ffdefa-d1b5-4cca-b71d-f7063065857e_atacseq_gdc_realn.bam" "${bam_dir}/summary_UCEC_8.bam"
			   
# ──────────────────────────────────────────────────────────────
# 2) Calculating sequencing depth and downscaling factor
# ──────────────────────────────────────────────────────────────			   

### Calculating sequencing depth for further downscaling
depth_dir="/home/personal/depth_dir"
   
samtools coverage "${bam_dir}/summary_UCEC_8.bam" -d 0 -o \
				  "${depth_dir}/summary_UCEC_8.txt"

#### Downscaling factor value
sum=$(awk 'NR >= 1 && NR <= 26 {sum += $4} END {print 250000000 / sum}' "${depth_dir}/summary_UCEC_8.txt")

### Indexing bam file			  
samtools index "${bam_dir}/summary_UCEC_8.bam"

echo "downscaling factor "$sum"" 

# ──────────────────────────────────────────────────────────────
# 3) Transforming bam file into BigWig file downscaling. 
# ──────────────────────────────────────────────────────────────

#### Downscaling experiment by its respective factor and making it a BigWig file
bw_dir="/home/personal/bw_dir"

bamCoverage -b "${bam_dir}/summary_UCEC_8.bam" \
            -o "${bw_dir}/summary_scaled_normal_UCEC_8.bw" \
            --binSize 20 --scaleFactor $sum --smoothLength 60 

bamCoverage -b "${bam_dir}/summary_UCEC_8.bam" \
            -o "${bw_dir}/summary_scaled_rpkm_UCEC_8.bw" \
            --binSize 20 --scaleFactor $sum --smoothLength 60 --normalizeUsing RPKM


# ...
# Perform the previous for each experiment that belongs to the Tissue - Heart                 

# ──────────────────────────────────────────────────────────────
# PART II - For every experiment that belongs to the respective Tissue
# ──────────────────────────────────────────────────────────────

# ──────────────────────────────────────────────────────────────
# 1) Merging downscaled BigWig files
# ──────────────────────────────────────────────────────────────

module load deepTools
final_bw_dir="/home/personal/final_bw"

bigwigAverage -b "${bw_dir}/summary_scaled_normal_UCEC_1.bw"
                 "${bw_dir}/summary_scaled_normal_UCEC_2.bw"
                 "${bw_dir}/summary_scaled_normal_UCEC_3.bw"
                 "${bw_dir}/summary_scaled_normal_UCEC_4.bw"
                 "${bw_dir}/summary_scaled_normal_UCEC_5.bw"
                 "${bw_dir}/summary_scaled_normal_UCEC_6.bw"
                 "${bw_dir}/summary_scaled_normal_UCEC_7.bw"
                 "${bw_dir}/summary_scaled_normal_UCEC_8.bw"   ## <-- This is our file example 
                  --binSize 20 --outFileFormat bigwig -o "${final_bw_dir}/summary_normal_UCEC.bw"
                  
bigwigAverage -b "${bw_dir}/summary_scaled_rpkm_UCEC_1.bw"
                 "${bw_dir}/summary_scaled_rpkm_UCEC_2.bw"
                 "${bw_dir}/summary_scaled_rpkm_UCEC_3.bw"
                 "${bw_dir}/summary_scaled_rpkm_UCEC_4.bw"
                 "${bw_dir}/summary_scaled_rpkm_UCEC_5.bw"
                 "${bw_dir}/summary_scaled_rpkm_UCEC_6.bw"
                 "${bw_dir}/summary_scaled_rpkm_UCEC_7.bw"
                 "${bw_dir}/summary_scaled_rpkm_UCEC_8.bw"   ## <-- This is our file example 
                  --binSize 20 --outFileFormat bigwig -o "${final_bw_dir}/summary_rpkm_UCEC.bw"
                  
                  