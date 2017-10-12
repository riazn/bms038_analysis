# RNAseq Differentially Expressed Gene (DEG) Analysis 

Code to reproduce key results of Table S6 tables from examination of differential
expression of genes between different sample groupings. Note, these scripts
require installation of package nrstats before running. Lastly, some scripts
include a variable PARENTDIR; which needs to be set appropriately for the script
to run. Results are outputted into "output" subdirectory. Please note version 1.12.4
of DESeq2 was used for these scripts, which is provided in the packages subdirectory

### Contents

1.) [TableS6.A.R](TableS6.A.R)
<br>Reproduces DEG results between responders and non-responders in Pre-therapy samples (Table S6A) 

2.) [TableS6.B.R](TableS6.B.R)
<br>Reproduces DEG results between genomic persistence and genomic contraction in Pre-therapy samples (Table S6B)
considering the previous ipilimumab treatment status (progressed or naive).

3.) [TableS6.C.R](TableS6.C.R)
<br>Reproduces DEG results between On- and Pre-therapy samples (Table S6C) considering the patient difference.

4.) [TableS6.D.R](TableS6.D.R)
<br>Reproduces DEG results between Pre- and On-therapy samples, considering genes that change differentially in long term benefit
versus no long term benefit patients (Table S6D)
