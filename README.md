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

* If you are a Windows user (like me) and need to run Linux within a virtual machine to perform analyses, check out [this](https://www.fosslinux.com/42789/how-to-install-linux-mint-on-virtualbox.htm) guide or [this](https://linuxhint.com/install_linux_mint_20_virtualbox/) one. I setup a shared folder (to transfer files between my host Windows PC and Virtual Machine Linux) using [this](https://helpdeskgeek.com/virtualization/virtualbox-share-folder-host-guest/) guide.

* After generation of figures, vector graphics were edited in Adobe Illustrator and assembled for publication in Adobe InDesign. Other photo edits were made in either Adobe Photoshop or Lightroom.

Scripts used to analyze datasets
----------------------------------------------------------------------
* [Aligning reads from RNA-seq experiment using HISAT2](https://github.com/gene-drive/Tbx5-forelimb-genital/blob/main/Scripts/HISAT2/HISAT2_alignment_RNA-seq.sh)
* [RNA-seq analysis pipeline to detect differentially expressed genes between Tbx5 cKO mutant forelimbs and controls](https://github.com/gene-drive/Tbx5-forelimb-genital/blob/main/Scripts/RNA-seq%20pipeline_DESeq2%20analysis%20and%20visualization_Forelimb%20Tbx5%20mutants%20vs%20controls.R)
* [Using bedtools to find peak coordinate intersections between ChIP-seq datasets](https://github.com/gene-drive/Tbx5-forelimb-genital/blob/main/Scripts/Peak-coordinate-intersections.sh)
* [Performing motif analysis of ChIP-seq peaks using HOMER](https://github.com/gene-drive/Tbx5-forelimb-genital/blob/main/Scripts/HOMER_motif-analysis.sh)


Scripts used to generate figures
----------------------------------------------------------------------
* [RNA-seq heatmap of global expression patterns](https://github.com/gene-drive/Tbx5-forelimb-genital/blob/main/Scripts/Heatmaps_RNA-seq.R). Example:

<img width="700" src="https://user-images.githubusercontent.com/61433004/231803053-dba321dc-dd72-4ceb-a44d-68cfe83c39bb.jpg">

* [Genome browser tracks of ChIP-seq peaks of interest](https://github.com/gene-drive/Tbx5-forelimb-genital/blob/main/Scripts/peak-visualization_pyGenomeTracks.sh) (and the [tracks file](https://github.com/gene-drive/Tbx5-forelimb-genital/blob/main/Scripts/hs1448_Wnt5a_Plot_v1.ini) used). Mouse embryo photo is from the [VISTA Enhancer Browser](https://enhancer.lbl.gov/). Example:

<img width="700" src="https://user-images.githubusercontent.com/61433004/231871420-1f878d17-3a93-4389-a703-1f555134b266.jpg">


* [Enrichment heatmaps (comparing bigwigs) of ChIP-seq datasets](https://github.com/gene-drive/Tbx5-forelimb-genital/blob/main/Scripts/Heatmaps_ChIP-seq_deepTools.sh). Example:
<img width="300" src="https://user-images.githubusercontent.com/61433004/231881114-baa6c5f1-e411-4918-9377-bee69cd0ca36.jpg">


Figures generated using easy-to-use webtools (no coding background needed!)
----------------------------------------------------------------------
* Volcano plots of differentially expressed genes found via RNA-seq - [ggVolcanoR](https://ggvolcanor.erc.monash.edu/). Example:

<img width="400" src="https://user-images.githubusercontent.com/61433004/232326917-9c7b6e8a-9a43-48c0-a24c-998cef8b0561.jpg">


* Gene-set enrichment analysis of differentially expressed genes from RNA-seq and generation of lollipop plots using [ShinyGO](http://bioinformatics.sdstate.edu/go/). Example:

<img width="400" src="https://user-images.githubusercontent.com/61433004/232326643-6ea914b4-aac7-4328-911a-1c9596fc8864.jpg">


* Assign biological meaning to a set of ChIP-seq binding sites by analyzing the annotations of nearby genes using [GREAT (Genomic Regions Enrichment of Annotations Tool)](http://great.stanford.edu/public/html/). Example:

<img width="400" src="https://user-images.githubusercontent.com/61433004/232327104-d57d000d-c2f6-47f7-a4a3-09c803e5e0d8.jpg">

