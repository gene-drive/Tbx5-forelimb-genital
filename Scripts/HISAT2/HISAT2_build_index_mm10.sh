### The script below was ran on the Sapelo2 cluster of the Georgia Advanced Computing Resource Center (GACRC) ###

#!/bin/bash
#SBATCH --job-name=mm10_index_HISAT2
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --export=NONE
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=user@uga.edu
#SBATCH --mail-type=ALL

cd $SLURM_SUBMIT_DIR

GENOME=/home/user/path-to-genome-file/mm10.fa

ml purge

module load HISAT2/2.1.0-foss-2019b

echo -e "\n** Script started on `date` **\n"

time hisat2-build ${GENOME} mm10_index_HISAT2

echo -e "\n** Script ended on `date` **\n"

echo Done!
