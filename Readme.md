---

---

## Single cell RNA sequencing data from Mouse PBMCs of various age

---

---

This directory contains files that allow you to reproduce the scRNA seq results from Paprckova et al. manuscript, titled:  
### Self-reactivity of CD8 T-cell clones determines their differentiation status rather than their responsiveness in infections  

Link to paper: COMING SOON!  

Link to data: [GEO208795](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208795)  

License: MIT. See LICENCE.md file in directory which applies to all files. Just be aware that all scripts and data are provided without any warranty and the autors take no responsibility and/or liability for any claims.

---

## Data Processing Scripts

### Requirements

The following software/packages are required for the reproduction of this analysis:  

- R 4.0.4 with package Seurat 4.0.3 and all its requirements (notably umap-learn)
- Cellranger 5.0.0 with:
- - GRCm38 v102 downloaded from Ensembl and prepared as reference 
- - VDJ reference built from IMGT data-base for mouse and including 5 extra sequences from file *IMGT_additional_VDJ.fasta* (the best way to build a reference is to generate .fasta following [10X page (--fetch-imgt)](https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/advanced/references)), then merge result with file IMGT_additional_VDJ.fasta and perform the rest of reference creation).
- MiXCR 3.0.12.
- bamtofastq from 10x ([get it here](https://support.10xgenomics.com/docs/bamtofastq))

Optional:

- Python 3.x.x for corrected barcode read demultiplexing. This is not necessary as the .fastq files were already demultiplexed, but you can have a look at the scripts.

### Data Treatment  

The data were processed using following steps, which can be splitted itno two parts:

**Part 1**

1. Initial mapping to GRCm38 reference using Cell Ranger:  
**cellranger count --id=Target_dir --feature-ref=FeatureReference.csv --libraries=Source_libs.csv --transcriptome=GRCm38_v102 --localcores=20** (this is for processing **GEX and csp data**; *Target_dir* will store the final results, *Source_libs.csv* contains the paths, sample prefixes towards samples and data type (either 'Gene Expression' or 'Feature Barcoding'). PBMCs from young and old mice need to be mapped separately, and GEX has 2 technical replicates, so in total there shoudl be three paths (for each of 2 GEX technical replicates and for csp). *FeatureReference.csv* defines csp barcodes.)
2. Mapping to expanded mouse IMGT reference using Cell Ranger:  
**cellranger vdj --id=Target_VDJ_dir --reference=IMGT_Mouse --fastqs=Fastqs --sample=sample_name --localcores=20** (for processing **V(D)J files**; *Target_VDJ_dir* will store the final results, *Fastqs* contains fastq files and *sample_name* is prefix of files in Fastqs directory for given sample. Again, young and old mice samples have to be mapped separately.)
3. The feature barcoding counts are used to establish a limit for number of csp reads that was used to split the data used in the paper from an unrelated experiment. These scripts are uploaded in directory *part_1_1_Identify_cells* (each for each subset). 
4. Extract reads from analyzed experiment using scripts contained in directory *part_1_2_Extract_reads*. This allows to account for reads with corrected cellular barcodes. There are two steps, preparation of .fastq and .bam (obtained by mapping) and extraction itself, which is little different for csp. **NOTE:** The available data have been already processed **up to this point**, these scripts are here merely for those interested.
5. Re-map resulting fastqs using respective commands for GEX+csp and V(D)J. The only difference is the fastq paths, otherwise commands from (1.) and (2.) are identical.
6. Process GEX with MiXCR to extract supplementary V(D)J.  **Before running script that runs MiXCR you must configure it with paths to MiXCR and 10X bamtfastq software ** (look at the beginning of *main* function). Also it's better to copy .bam obtained by 10X and its index to different directory, because it will create new files. After that, run MiXCR with command *python Mixcr_pipeline_II.py source.bam result_dir*, where source.bam is copied bam along with its index. Result_dir then should contain directory mixcr where two .csv files should be present, one for TRA and TRB chains each. Create directory VDJ_GEX in directory where you'll run next step (7.) and copy them there. You need to run this script for both young and old B6 mice data, so in the end VDJ_GEX directory should contain 4 files that you should rename like this:  
* * B6_Young_blood_TRA.csv
* * B6_Young_blood_TRB.csv
* * B6_Old_blood_TRA.csv
* * B6_Old_blood_TRB.csv
7. Prepare data using scripts contained in directory *part_1_3_Adjacent_sources*. The main script is named *part_1_3_Data_set_preparation.Rmd*, which runs everything needed. **Precision:** since filtering by limit on gene present in number of cells was done on combined data-set also containing cells from unrelated experiment, this filtering is not possible on the data without these cells. For this reason we provide a list of genes which have to be present in the final data. You can also follow the fitlering as usual; the results will be slightly different, but the final conclusion should be the same.  

**Part 2**

In Part 2, we describe:

* Additional quality control of the data
* Integration of samples from young and aged mice
* Normalization, scaling, variable feature selection, dimensional reduction
* Differential expression testing
* Generation of figures for the manuscript

Please download and open the file `part_2_Final_analysis.html` to see all the final figures included in the manuscript with the code to produce them. After the file is downloaded and saved, it can be viewed easily in any web browser. Please note that it is crucial to download the file to your laptop first and open it from there. You will only see html code if you try to look at it directly from GitHub. 
---

Contact:  

Should you have any question, feel free to write to:  
- veronika.niederlova@img.cas.cz  
- juraj.michalik@img.cas.cz  

---


