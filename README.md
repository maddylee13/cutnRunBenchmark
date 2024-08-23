# CUT&Run Benchmark

This repository is split into 2 major parts: (i) the processing and analysis of CUT&Run (Cleavage Under Targets & Release using nucleotides) data and (ii) the benchmarking of CUT&Run data against other more established methods. Any analysis or benchmarking using R was done in R (v.4.3.1).

**CUT&Run Data Analysis** 

The processing and analysis of the CUT&Run data follows the tutorial provided by the Henikoff group which can be found at (https://dx.doi.org/10.17504/protocols.io.bjk2kkye). All code used for this analysis is enclosed and labelled in numerical order. Scripts used for this analysis were adapted from (https://github.com/clabanillas/cutnTag_processing). 

**Benchmarking**

The CUT&Run data was benchmarked against publicly available ChIPseq datasets. The code for functions used to generate the UpSet plot was generated from publicly available code from (https://github.com/neurogenomics/EpiCompare). Specifically, the "messager" and "overlap_upset_plot" functions were implemented from the EpiCompare package. 

**R Packages and Command Line Tools**

The following R packages were used:
* BiocManager (v.1.30.23)
* ChIPpeakAnno (v.3.36.0)
* ChIPseeker (v.1.38.0)
* chromVAR (v.1.24.0)
* ComplexUpset (v.1.3.3)
* DiffBind (v.3.12.0)
* dplyr (v.1.1.3)
* GenomicFeatures (v.1.54.1)
* GenomicRanges (v.1.54.1)
* ggplot2 (v.3.4.4)
* ggpubr (v.0.6.0)
* tidyr (v.1.3.0)
* viridis (v.0.6.5)

The following Anaconda evironments and command line tools were used:
| Env   |bowtie2            |cutnTag            |deeptools          | macs2         | seacr             |ucsc-tools        |
|-------|-------------------|-------------------|-------------------|---------------|-------------------|------------------|
|tools  |bowtie2 (v.2.5.3)  |bedtools (v.2.31.1)|bedtools (v.2.31.1)|macs2 (v.1.9.4)|bedtools (v.2.31.1)|ucsc-tools (v.357)|
|       |samtools (v.1.2)   |fastp (v.0.23.4)   |bowtie2 (v.2.5.3)  |               |seacr (v.1.3)      |                  |
|       |                   |fastqc (v.0.12.1)  |deeptools (v.3.5.5)|               |                   |                  |
|       |                   |multiqc (v.1.2.1)  |samtools (v.1.2)   |               |                   |                  |
|       |                   |picard (v.2.18.7)  |fastqc (v.0.12.1)  |               |                   |                  |
|       |                   |samtools (v.1.2)   |                   |               |                   |                  |
|       |                   |bioconda


**References**
