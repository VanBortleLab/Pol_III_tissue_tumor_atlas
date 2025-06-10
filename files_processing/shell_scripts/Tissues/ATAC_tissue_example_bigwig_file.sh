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
# 1) Merging replicates per experiment
# ──────────────────────────────────────────────────────────────

align_dir='/home/personal/align_dir' ## All the bam files for each experiment and its replicates is here

bam_dir='/home/personal/bam_dir' ## Output for merged bam file per experiment

samtools merge "${bam_dir}/summary_Heart_ENCSR925LGW.bam" \
			   "${align_dir}/ENCSR925LGW.heart_left_ventricle.1_1_1.sorted.bam" \
			   "${align_dir}/ENCSR925LGW.heart_left_ventricle.1_1_2.sorted.bam" \
			   "${align_dir}/ENCSR925LGW.heart_left_ventricle.1_1_3.sorted.bam" \
			   "${align_dir}/ENCSR925LGW.heart_left_ventricle.1_2_1.sorted.bam" \
			   "${align_dir}/ENCSR925LGW.heart_left_ventricle.1_2_2.sorted.bam" \
			   "${align_dir}/ENCSR925LGW.heart_left_ventricle.1_2_3.sorted.bam" \
			   "${align_dir}/ENCSR925LGW.heart_left_ventricle.1_2_4.sorted.bam" 
			   
# ──────────────────────────────────────────────────────────────
# 2) Calculating sequencing depth and downscaling factor
# ──────────────────────────────────────────────────────────────			   

### Calculating sequencing depth for further downscaling
depth_dir="/home/personal/depth_dir"
			   
samtools coverage "${bam_dir}/summary_Heart_ENCSR925LGW.bam" -d 0 -o \
				  "${depth_dir}/summary_Heart_ENCSR925LGW.txt"

#### Downscaling factor value
sum=$(awk 'NR >= 1 && NR <= 26 {sum += $4} END {print 250000000 / sum}' "${depth_dir}/summary_Heart_ENCSR925LGW.txt")

### Indexing bam file			  
samtools index "${bam_dir}/summary_Heart_ENCSR925LGW.bam"

echo "downscaling factor "$sum"" 

# ──────────────────────────────────────────────────────────────
# 3) Transforming bam file into BigWig file downscaling. 
# ──────────────────────────────────────────────────────────────

#### Downscaling experiment by its respective factor and making it a BigWig file
bw_dir="/home/personal/bw_dir"

bamCoverage -b "${bam_dir}/summary_Heart_ENCSR925LGW.bam" \
            -o "${bw_dir}/summary_scaled_Heart_ENCSR925LGW.bw" \
            --binSize 20 --scaleFactor $sum --smoothLength 60 

bamCoverage -b "${bam_dir}/summary_Heart_ENCSR925LGW.bam" \
            -o "${bw_dir}/summary_scaled_rpkm_Heart_ENCSR925LGW.bw" \
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

bigwigAverage -b "${bw_dir}/summary_scaled_Heart_ENCSR074WMH.bw" \
                 "${bw_dir}/summary_scaled_Heart_ENCSR133CMC.bw" \
                 "${bw_dir}/summary_scaled_Heart_ENCSR167SYF.bw" \
                 "${bw_dir}/summary_scaled_Heart_ENCSR204PZT.bw" \
                 "${bw_dir}/summary_scaled_Heart_ENCSR204SMO.bw" \
                 "${bw_dir}/summary_scaled_Heart_ENCSR925LGW.bw" \
                 "${bw_dir}/summary_scaled_Heart_ENCSR925LGW.bw" \
				 "${bw_dir}/summary_scaled_Heart_ENCSR212LYK.bw" \
				 "${bw_dir}/summary_scaled_Heart_ENCSR310RJN.bw" \
				 "${bw_dir}/summary_scaled_Heart_ENCSR310UDW.bw" \
				 "${bw_dir}/summary_scaled_Heart_ENCSR327DCG.bw" \
				 "${bw_dir}/summary_scaled_Heart_ENCSR330WUK.bw" \
 				 "${bw_dir}/summary_scaled_Heart_ENCSR374ZRW.bw" \
				 "${bw_dir}/summary_scaled_Heart_ENCSR390SLL.bw" \
				 "${bw_dir}/summary_scaled_Heart_ENCSR437OOJ.bw" \
				 "${bw_dir}/summary_scaled_Heart_ENCSR439TZT.bw" \
				 "${bw_dir}/summary_scaled_Heart_ENCSR451JSB.bw" \
				 "${bw_dir}/summary_scaled_Heart_ENCSR452OSK.bw" \
				 "${bw_dir}/summary_scaled_Heart_ENCSR522FGI.bw" \
				 "${bw_dir}/summary_scaled_Heart_ENCSR710SMN.bw" \
				 "${bw_dir}/summary_scaled_Heart_ENCSR731ODJ.bw" \
				 "${bw_dir}/summary_scaled_Heart_ENCSR736CIN.bw" \
 				 "${bw_dir}/summary_scaled_Heart_ENCSR769DGC.bw" \
				 "${bw_dir}/summary_scaled_Heart_ENCSR822BAD.bw" \
				 "${bw_dir}/summary_scaled_Heart_ENCSR848TMJ.bw" \
				 "${bw_dir}/summary_scaled_Heart_ENCSR851EBF.bw" \
				 "${bw_dir}/summary_scaled_Heart_ENCSR855BMI.bw" \
				 "${bw_dir}/summary_scaled_Heart_ENCSR925LGW.bw" \ ## <-- This is our file example 
                  --binSize 20 --outFileFormat bigwig -o "${final_bw_dir}/summary_normal_Heart.bw"
                  
bigwigAverage -b "${bw_dir}/summary_scaled_rpkm_Heart_ENCSR074WMH.bw" \
                 "${bw_dir}/summary_scaled_rpkm_Heart_ENCSR133CMC.bw" \
                 "${bw_dir}/summary_scaled_rpkm_Heart_ENCSR167SYF.bw" \
                 "${bw_dir}/summary_scaled_rpkm_Heart_ENCSR204PZT.bw" \
                 "${bw_dir}/summary_scaled_rpkm_Heart_ENCSR204SMO.bw" \
                 "${bw_dir}/summary_scaled_rpkm_Heart_ENCSR925LGW.bw" \
                 "${bw_dir}/summary_scaled_rpkm_Heart_ENCSR925LGW.bw" \
				 "${bw_dir}/summary_scaled_rpkm_Heart_ENCSR212LYK.bw" \
				 "${bw_dir}/summary_scaled_rpkm_Heart_ENCSR310RJN.bw" \
				 "${bw_dir}/summary_scaled_rpkm_Heart_ENCSR310UDW.bw" \
				 "${bw_dir}/summary_scaled_rpkm_Heart_ENCSR327DCG.bw" \
				 "${bw_dir}/summary_scaled_rpkm_Heart_ENCSR330WUK.bw" \
 				 "${bw_dir}/summary_scaled_rpkm_Heart_ENCSR374ZRW.bw" \
				 "${bw_dir}/summary_scaled_rpkm_Heart_ENCSR390SLL.bw" \
				 "${bw_dir}/summary_scaled_rpkm_Heart_ENCSR437OOJ.bw" \
				 "${bw_dir}/summary_scaled_rpkm_Heart_ENCSR439TZT.bw" \
				 "${bw_dir}/summary_scaled_rpkm_Heart_ENCSR451JSB.bw" \
				 "${bw_dir}/summary_scaled_rpkm_Heart_ENCSR452OSK.bw" \
				 "${bw_dir}/summary_scaled_rpkm_Heart_ENCSR522FGI.bw" \
				 "${bw_dir}/summary_scaled_rpkm_Heart_ENCSR710SMN.bw" \
				 "${bw_dir}/summary_scaled_rpkm_Heart_ENCSR731ODJ.bw" \
				 "${bw_dir}/summary_scaled_rpkm_Heart_ENCSR736CIN.bw" \
 				 "${bw_dir}/summary_scaled_rpkm_Heart_ENCSR769DGC.bw" \
				 "${bw_dir}/summary_scaled_rpkm_Heart_ENCSR822BAD.bw" \
				 "${bw_dir}/summary_scaled_rpkm_Heart_ENCSR848TMJ.bw" \
				 "${bw_dir}/summary_scaled_rpkm_Heart_ENCSR851EBF.bw" \
				 "${bw_dir}/summary_scaled_rpkm_Heart_ENCSR855BMI.bw" \
				 "${bw_dir}/summary_scaled_rpkm_Heart_ENCSR925LGW.bw" \ ## <-- This is our file example 
                  --binSize 20 --outFileFormat bigwig -o "${final_bw_dir}/summary_rpkm_Heart.bw"
                  
                  