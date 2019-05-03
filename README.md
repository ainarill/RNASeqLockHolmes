# RNASeqLockHolmes
## Pipeline for differential expression analysis of RNAseq data.

*Luisa Santus, Aina Rill and Altaïr C.Hernández *

## **Index**

<!-- TOC depthFrom:1 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->

- [Data Structure](#data-structure)

<!-- /TOC -->

## Data Structure

      RNASeqLockHolmes/
            README.md
            bibliography.bib
            index.Rmd
            QAanalysis.Rmd
            Makefile
            QAanalysis_files/
            DEanalysis_files/
            rawCounts/
                seCOAD.rds
            REFS/
            results/


* README.md
* bibliography.bib: BibTeX file format containing the references of this project.
* index.Rmd: R markdown file format for the index file.
* QAanalysis.Rmd
* Makefile: to compile the rmd into HTML format.
* QAanalysis_files/: Quality Assesment and Normalization results of the RNA-seq data.
* DEanalysis_files/: Differential Expression analysis of the RNA-seq data.
* rawCounts/: raw data of the whole experiment taken b
* REFS/
* results/

### RNASeqLockHolmes

RNASeqLockHolmes is a project of RNA-seq analysis of Colorectal Cancer(CRC) developed by **Luisa Santus**, **Aina Rill** and  **Altaïr C. Hernández** (*Master of Bioinformatics for Health Science*).

Colorectal Cancer (CRC) is the second most common cancer in women (614,000 cases per year) and the third most common in men(746,000 cases per year). Currently it is the most common malignant cancer in the gatrointestinal tract, representing 13\% of all malignant tumors, and it is considered the main cause of death in gastrointestinal cancer @Siegel2015. It can be arised from one or a combination of three differen mechanisms, namely chromosmomal inestability (CIN), CpG islands methylator phenotype (CIMP), and microsatellite inestability (MSI). 

Understanding the specific mechanisms of tumorigenesis and the underlaying genetic and epigenetic traits is crutial in the disease phenotype comprehension. [The Cancer Genome Atlas](https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga) (TCGA) has comprehensively profiled this type of cancer in a patient cohort. Here we analyze the expression profiles of those patients, accessible in the form of a raw RNA-seq counts produced by @rahman2015alternative using a pipeline based on the R/Bioconductor software package `r Biocpkg("Rsubread")`.





