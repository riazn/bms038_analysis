# bms038 analysis

Data and associated code to perform analysis in manuscript
["Tumor and Microenvironment Evolution during Immune Checkpoint Blockade Therapy with 
Nivolumab"](https://doi.org/10.1016/j.cell.2017.09.028). Each sub-directory contains code to reproduce relevant portion of analysis. 
See READMEs in each sub-directory for additional details. 
Description of how the mutation list and RNAseq count matrix are
generated from raw BAM files are provided in the methods section of the manuscript. 

## Contents
1. [Exome](exome/)
<br>R scripts to generate results for mutation load and clonal mutation load

2. [Clonality_pre_post](clonality_pre_post/)
<br>Scripts to analyze change in clonality of tumor pre and on therapy

3. [RNASeq](rnaseq/)
<br>Scripts to analyze RNA-seq data

4. [TCR](tcr/)
<br>Code to produce results regarding TCR analysis

5. [Data](data/)
<br>Data from genomic pipelines (mutation calls, RNAseq matrix, etc)

6. [Output](output/)
<br>All scripts save results to this directory

7. [Packages](packages/)
<br>Custom code/packages necessary to reproduce different portions of analysis


