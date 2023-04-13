######### This document describes how I used pyGenomeTracks to visualize genome browser tracks of ChIP-seq peaks of interest.

# pyGenomeTracks Documentation: https://pygenometracks.readthedocs.io/en/latest/

### First make sure you install Anaconda in Linux Mint: https://docs.anaconda.com/anaconda/install/linux/

### Then, install pyGenomeTracks:
pip install pyGenomeTracks

# Make sure to install bedtools, as it is required for pyGenomeTracks:
apt-get install bedtools


### Before running pyGenomeTracks, you first need to create a tracks.ini file.
# You can use this guide: https://pygenometracks.readthedocs.io/en/latest/content/usage.html
# Or see this file for an example: https://github.com/gene-drive/Tbx5-forelimb-genital/blob/main/Scripts/hs1448_Wnt5a_Plot_v1.ini
# The example contains tracks for: scalebar, TBX5 mouse forelimb ChIP-seq bigwig, TBX5 mouse forelimb ChIP-seq peak, TBX5 amniote-conserved forelimb ChIP-seq peak, VISTA-limb-positive-enhancer, Andrey_Forelimb-CaptureC, and phastCons Conservation.

### Now you can run this command to generate the plot:

pyGenomeTracks --tracks hs1448_Wnt5a_Plot_v1.ini \
--region chr14:28,970,650-28,975,348 \
--width 12 --dpi 300 --fontSize 6 \
--trackLabelFraction 0.1 \
-o hs1448_Wnt5a_Plot.svg

# See here for a description of all the arguments: https://pygenometracks.readthedocs.io/en/latest/content/usage.html




# After generation of plots, vector graphics were edited in Adobe Illustrator. Final figures were assembled for publication in Adobe InDesign.
