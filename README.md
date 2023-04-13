# Aaron Alcala's scripts for ChIP-seq and RNA-seq analyses

This repository contains scripts used in the following Ph.D. disseration:

**Limbs, Genitals, and Gene Regulatory Networks: Identifying Conserved Targets of TBX5 in Embryonic Forelimbs and Genitalia**

by Aaron Alcala

* [Website](https://aaronevodevo.wixsite.com/aaronevodevo)
* [SciComm Portfolio](https://aaronalcala.myportfolio.com/)
* [LinkedIn](https://www.linkedin.com/in/aaronalcala/)
* [Social Media Links - @AaronEvoDevo](https://linktr.ee/AaronAlcala)


Data Availability
----------------------------------------------------------------------
Genomic and transcriptomic data for this manuscript will be uploaded to the Gene Expression Omnibus database and updated accession numbers will be posted here.

Tools used
----------------------------------------------------------------------
For a list of all tools used in this manuscript, click [here](https://github.com/gene-drive/Tbx5-forelimb-genital/blob/main/Tools-Used.md).

* If you are a Windows user (like me) and need to run Linux within a virtual machine to perform analyses, check out [this](https://www.fosslinux.com/42789/how-to-install-linux-mint-on-virtualbox.htm) guide or [this](https://linuxhint.com/install_linux_mint_20_virtualbox/) one. I setup a shared folder this [this](https://helpdeskgeek.com/virtualization/virtualbox-share-folder-host-guest/) guide.

* After generation of figures, vector graphics were edited in Adobe Illustrator and assembled for publication in Adobe InDesign. Other photo edits were made in either Adobe Photoshop or Lightroom.

Scripts used to analyze datasets
----------------------------------------------------------------------
* [RNA-seq analysis pipeline comparing Tbx5 cKO mutant forelimbs to controls](https://github.com/gene-drive/Tbx5-forelimb-genital/blob/main/Scripts/RNA-seq%20pipeline_DESeq2%20analysis%20and%20visualization_Forelimb%20Tbx5%20mutants%20vs%20controls.R)
* [Using bedtools to find peak coordinate intersetions between ChIP-seq datasets](https://github.com/gene-drive/Tbx5-forelimb-genital/blob/main/Scripts/Peak-coordinate-intersections.sh)


Scripts used to generate figures
----------------------------------------------------------------------
* [RNA-seq heatmap of global expression patterns](https://github.com/gene-drive/Tbx5-forelimb-genital/blob/main/Scripts/Heatmaps_RNA-seq.R). Example:

<img width="700" src="https://user-images.githubusercontent.com/61433004/231803053-dba321dc-dd72-4ceb-a44d-68cfe83c39bb.jpg">

* [Genome browser tracks of ChIP-seq peaks of interest](https://github.com/gene-drive/Tbx5-forelimb-genital/blob/main/Scripts/peak-visualization_pyGenomeTracks) (and the [tracks file](https://github.com/gene-drive/Tbx5-forelimb-genital/blob/main/Scripts/hs1448_Wnt5a_Plot_v1.ini) used). Example:

<img width="700" src="https://user-images.githubusercontent.com/61433004/231871420-1f878d17-3a93-4389-a703-1f555134b266.jpg">
