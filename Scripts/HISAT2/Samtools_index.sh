#!/bin/bash
#SBATCH --job-name=FL_ctrl1_bai_redo
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=3:00:00
#SBATCH --export=NONE
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=user@uga.edu 
#SBATCH --mail-type=ALL

cd $SLURM_SUBMIT_DIR

ml purge
module load SAMtools/1.10-iccifort-2019.5.281

echo -e "\n** Script started on `date` **\n"

time samtools index hisat2_FL_ctrl1_redo.bam

echo -e "\n** Script ended on `date` **\n"
echo Done!
