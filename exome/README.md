# Baseline Exome Analysis 

Code to reproduce key results from examination of baseline exome data,
specifically evaluations looking at baseline mutation load and baseline
clonal mutation load. Note, these scripts require installation of
package nrstats before running. Lastly, some scripts include a variable
PARENTDIR; which needs to be set appropriately for the script to run.
Results are outputted into "output" subdirectory.


### Contents

1.) [baseline_analysis.r](baseline_analysis.r)
<br>Reproduces portion of Fig 1B, and Supp Fig 1A, along with several
p-values discussed in main text. 

2.) [baseline_clonality.r](baseline_clonality.r)
<br>Reproduces portion of Fig 1B, and Supp Fig 1D, 1B, and 1C. 
must run baseline_analysis.r beforehand. 
