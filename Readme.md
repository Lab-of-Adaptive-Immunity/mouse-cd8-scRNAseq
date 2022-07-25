---

---

## Single cell RNA sequencing data from Mouse PBMCs of various age

---

---

This directory contains files that allow you to reproduce the scRNA seq results from Paprckova et al. manuscript, titled:  

###Self-reactivity of CD8 T-cell clones determines their differentiation status rather than their responsiveness in infections  

Link to paper: COMING SOON!  

Link to data: [GEO208795](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208795)  

---

## Data Processing Scripts

### Requirements

The following software/packages are required for the reproduction of this analysis:  

- R 4.0.4 with package Seurat 4.0.3 and all its requirements (notably umap-learn)
- Cellranger 5.0.0 with:
- - GRCm38 v102 downloaded from Ensembl and prepared as reference 
- - VDJ reference built from IMGT data-base for mouse and including 5 extra sequences from file *IMGT_additional_VDJ.fasta* (the best way to build a reference is to generate .fasta following [10X page (--fetch-imgt)](https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/advanced/references)), then merge result with file IMGT_additional_VDJ.fasta and perform the rest of reference creation).
- MiXCR 3.0.12

Optional:

- Python 3.x.x for corrected barcode read demultiplexing. This is not necessary as the .fastq files were already demultiplexed, but you can have a look at the scripts.

### Data Treatment  

The data are treated in following way.  

1. Initial mapping to GRCm38 reference using Cell Rannger:  
**cellranger count --id=Target_dir --feature-ref=FeatureReference.csv --libraries=Source_libs.csv --transcriptome=GRCm38_v102 --localcores=20** (this is for processing **GEX and csp data**; *Target_dir* will store the final results, *Source_libs.csv* contains the paths, sample prefixes towards samples and data type (either 'Gene Expression' or 'Feature Barcoding'). PBMCs from young and old mice need to be mapped separately, and GEX has 2 technical replicates, so in total there shoudl be three paths (for each of 2 GEX technical replicates and for csp). *FeatureReference.csv* defines csp barcodes.)
2. Mapping to expanded mouse IMGT reference using Cell Ranger:  
**cellranger vdj --id=Target_VDJ_dir --reference=IMGT_Mouse --fastqs=Fastqs --sample=sample_name --localcores=20** (for processing **V(D)J files**; *Fastqs* contains fastq files and *sample_name* is prefix of files in Fastqs directory for given sample. Again, young and old mice samples have to be mapped separately.)
3. Using feature barcoding.





