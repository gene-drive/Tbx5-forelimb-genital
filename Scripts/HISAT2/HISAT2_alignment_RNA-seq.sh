### The script below was ran on the Sapelo2 cluster of the Georgia Advanced Computing Resource Center (GACRC) ###

## First, we need to build an index of the mm10 genome using the following script:
# https://github.com/gene-drive/Tbx5-forelimb-genital/blob/main/Scripts/HISAT2/HISAT2_build_index_mm10.sh

## After you run the below script, you may need to create an indexed bam file (.bam is the final output file of HISAT2 after you align your .fastq input). Use the following script to index the .bam:
# https://github.com/gene-drive/Tbx5-forelimb-genital/blob/main/Scripts/HISAT2/HISAT2_build_index_mm10.sh


### Below is the script to align .fastq reads to the mm10.fa genome using HISAT2:

#!/bin/bash
#SBATCH --job-name=hisat2_FL_ctrl1_v4_redo
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --time=24:00:00
#SBATCH --export=NONE
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=user@uga.edu 
#SBATCH --mail-type=ALL

cd $SLURM_SUBMIT_DIR

GENOME_INDEX=/lustre2/scratch/user/path-to-index/mm10_index_HISAT2
FASTQ=/lustre2/scratch/user/path-to-fastq-files/AA_mm_FL_E10_5_RNA_Tbx5_ctrl1_S1_L001_R1_001.fastq.gz

ml purge
module load HISAT2/2.1.0-foss-2019b
module load SAMtools/1.10-iccifort-2019.5.281

echo -e "\n** Script started on `date` **\n"
echo -e "\n** Mapping reads to reference genome **\n"

time hisat2 \
--summary-file hisat2_FL_ctrl1_redo_summary.txt \
-x ${GENOME_INDEX} \
-U ${FASTQ} \
-S hisat2_FL_ctrl1_redo.sam

echo -e "\n** convert output sam to sorted bam and extract mapped reads **\n"

time samtools view -Sb hisat2_FL_ctrl1_redo.sam | samtools sort - -T hisat2_FL_ctrl1_redo.sorted -o hisat2_FL_ctrl1_redo.bam

echo -e "\n** Script ended on `date` **\n"
echo Done!
