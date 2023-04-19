### The script below was ran on the Sapelo2 cluster of the Georgia Advanced Computing Resource Center (GACRC) ###

### Here is documentation for HOMER Motif analysis:
# http://homer.ucsd.edu/homer/motif/index.html
# http://homer.ucsd.edu/homer/ngs/peakMotifs.html

### Before running, you need to take your output from peak calling and create a .bed containing Column1-chromosome name, Column2-starting position, Column 3-ending position, and Column4-peak name. I centered peaks based on summits and extended each peak 50 bp in each direction (each final peak used for HOMER were 100 bp total).

### Below is the script using default lengths of motif to be found. If you need to specify a length, you can use the -len parameter.


#!/bin/bash
#SBATCH --job-name=v1_c_mm10_FL_Homer
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --time=4:00:00
#SBATCH --export=NONE
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=user@uga.edu
#SBATCH --mail-type=ALL

cd $SLURM_SUBMIT_DIR

ml Homer/4.11-foss-2019b

echo -e "\n** Script started on `date` **\n"

findMotifsGenome.pl \
/scratch/user/path-to-bed-files/mm10_FL_Tbx5_9153711_conservative_peak_for_HOMER_100bp_top2000.bed \
mm10 \
v1_c_mm10_FL_HomerOutput/ \
-size 100 \
-mask \
-S 15 \
-N 100000 \
-preparsedDir ./preparsed \

echo -e "\n** Script ended on `date` **\n"
echo Done!
