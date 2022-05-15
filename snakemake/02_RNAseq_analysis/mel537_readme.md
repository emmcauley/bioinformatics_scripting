# Lui, Moore 2022: RNAseq Analysis, Continuation

## Scientific Summary
[Lui, Moore et al., 2022](https://onlinelibrary.wiley.com/doi/10.1111/pcmr.13013) differentiates YAP from its ortholog TAZ in melanoma progression. YAP plays a more crucial role than TAZ in driving pro-tumorigenic phenotypes in melanoma cells (Mel357), suggesting that YAP is a driver of melanoma progression, migration, and invasion. In this manuscript, the authors provide support that YAP drives melanoma migration through YAP-specific regulation of ARPC5, shifting ARP2/3 dynamics toward a pro-migratory phenotype.

## Technical Summary
[Snakemake](https://github.com/snakemake/snakemake) is a workflow manager to help researchers create reproducible and scalable data analyses. Here we present a Snakemake pipeline to analyze RNA-sequencing data in mel357 and skmel-5 cell lines. The goal of the pipeline is to examine differentially expressed genes as a result of YAP and/or TAZ inhibition. The pipeline includes: FastQC->Cut Adapt->FactQC->HiSat2->Bioconductor Feature Count->DESeq2. We first construct and then validate the pipeline by comparing Snakemake results to results published in Lui, Moore et al., 2022. In this documemntation we describe all steps for procuring data from NCBI GEO, constructing the Snakemake pipeline, and then comparing Snakemake-generated results to published results.

### Step 1. Fetch Relevant Files
Data in Lui, Moore et al., 2022 were made available in the NCBI GEO database under accession number [GSE146918](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146918). We exported the SRR Accession list from [SRA Run Selector](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA612430&o=acc_s%3Aa) and used `prefetch` from [sra-tools](https://hpc.nih.gov/apps/sratoolkit.html) in a Conda virtual environment to fetch the .sra files. We next used `fastq-dump` to convert to .fastq files.

###### Specific Bash Commands

1. Prefetch: `cat SRR_Acc_List.txt | while read line; do prefetch $line; done`

2. Fastq conversion:`cat SRR_Acc_List.txt | while read srr; do fastq-dump --outdir fastq --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip $srr; done`

### Step 2. Create Snakemake Pipeline Components
Herein we use the [Snakemake workflow manager](https://snakemake.readthedocs.io/) to execute the rest of our pipeline: FastQC, CutAdapt, HiSat2, FeatureCount, and DESeq2.  

#### Required Conda Environments
For each tool in the pipeline, a conda virtual environment was created to house that tool. A corresponding .yaml file was created from each environment and included in the Snakefile.

###### Specific Bash Commands

1. Environment creation: `conda create -n <name>`
2. YAML file creation: `conda env export > <name>.yml`

### Appendix

#### Conda environment
