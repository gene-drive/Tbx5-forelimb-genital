######### This document describes how I used bedtools to find peak coordinate intersetions between ChIP-seq datasets.

# bedtools Documentation: https://bedtools.readthedocs.io/en/latest/index.html

# Install bedtools in Linux Mint
apt-get install bedtools

# The input files for this tool are .bed peak coordinates. Column 1 is chromosome (chr), Column 2 is starting position, Column 3 is ending position, and Column 4 is peak name. Other columns are optional.
# ***Make sure the 4th column in your .bed file has peak names!
# It's recommended to sort large files, so that the algorithm used to perform intersections occurs quicker.


### In Microsoft Excel, you can use the following function (paste it into a cell) to count duplicates in a column:
# =SUMPRODUCT((A2:A100<>"")/COUNTIF(A2:A100,A2:A100&"")-(COUNTIF(A2:A100,A2:A100&"")=1))


######### Using the intersect command

### below are command options to know
# -a is standard to choose file 1
# -b is standard for 2nd file to intersect
# -wa	Write the original entry in A for each overlap.
# -wb	Write the original entry in B for each overlap. Useful for knowing what A overlaps. Restricted by -f and -r
# -u	Write original A entry once if any overlaps found in B. In other words, just report the fact at least one overlap was found in B. 
# -wo	Write the original A and B entries plus the number of base pairs of overlap between the two features. Only A features with overlap are reported
# -sorted	For very large B files, invoke a “sweeping” algorithm that requires position-sorted (e.g., sort -k1,1 -k2,2n for BED files) input. When using -sorted, memory usage remains low even for very large files.
# -v	Only report those entries in A that have no overlap in B




### Sample command used to output which TBX5 mouse forelimb peaks (13,580 peaks) intersect with TBX5 mouse genital tubercle peaks (20,728 peaks):

bedtools intersect \
-a mm10_FL_Tbx5_9153711_idr_conservative_peak.bed \
-b mm10_GT_Tbx5_Sung_conservative_peak.bed \
-wa -u \
> mm10_FL_vs_GT_TBX5_intersected_v2.bed

# 4,388 TBX5 peaks overlap between the forelimb and genital tubercle

# The -wa command ensures only peaks from -a (forelimb peaks) are output. The -u command ensures the original entry from -a is output only once.
