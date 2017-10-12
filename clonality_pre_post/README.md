# Pre and On therapy clonality analysis

Code to reproduce key plots shown in Figure 2 of manuscript. 
Please note R version 3.3.0 or greater with the following packages are required:
RColorBrewer, survival, and ks 
<br>Some scripts include a variable PARENTDIR; which needs to be set
appropriately for the script to run. Results are outputted into "output"
subdirectory. Please note custom package nrstats is also required and available
in the packages subdirectory


### Contents

1.) [clonality.analysis.R](clonality.analysis.R)
<br>Code to generate Fig 2a, 2b, 2d

2.) [clonality.surv.R](clonality.surv.R)
<br>generates Fig 2c

3.) [clonality.utils.R](clonality.utils.R)
<br>utility functions used by both clonality.analyis.R and clonality.surv.R