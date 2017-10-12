###### Analyze TCR seq data ######
#
# Code to reproduce TCR results
#   figure 5B and 5C
rm(list=ls())
# CHANGE This to directory downloaded to
PARENTDIR <- "/mnt/skimcs/Melanoma/BMS_038/forGithub/bms038_analysis/"
setwd(PARENTDIR)

library(nrstats)
library(reshape2)


# Figure 5B
# Note data here is placed in one csv file for convience
# Data were parsed from the original adaptive files using JS_splitter_Imsqv4.5.py 
# Richness & evenness were computed using JS_entropy_v7.7.py 
# data were subsequently summarized using linepick_7.7.2.py and the file used here
# includes data from linepick_entropy.out. 
#The column headers of linepick_entropy.out are as follows: 
#  FileName
#Hcdr3	
#Hvj	
#Htot	
#CLcdr3	
#CLvj	
#CLtot	
#Hcdr3_max	
#Hvj_max	
#Htot_max	
#Num_CDR3	
#Num_VJ	
#Num_totCDR3
#
#Values from the columns titled ???Num_CDR3??? and ???CLcdr3??? and were used to generate Figure 5B. Note that clonality values (CLcdr3) were converted to evenness using the relationship:
#  
#  evenness = 1- CLcdr3

tcr_sum_data <- read.csv("data/tcr/TCRmetrics.csv", stringsAsFactors = FALSE)

# Create NB vs LB (PRCR & SD)
tcr_sum_data$Response <- rep("Benefit", dim(tcr_sum_data)[1])
i <- which(tcr_sum_data$myBOR == "PD")
tcr_sum_data$Response[i] <- "NoBenefit"

tcrIpiProg <- subset(tcr_sum_data, Cohort == "NIV3-PROG")
tcrIpiNaive <- subset(tcr_sum_data, Cohort == "NIV3-NAIVE")

p<-tTestandGraph(tcrIpiProg, "Response", c("Benefit","NoBenefit"),"Num_CDR3_FoldChange", c("Benefit", "No Benefit"), tit="Ipi Prog: Richness")
ggsave(file="output/Fig5B_part1.png", p)

p<-tTestandGraph(tcrIpiNaive, "Response", c("Benefit","NoBenefit"),"Num_CDR3_FoldChange", c("Benefit", "No Benefit"), tit="Ipi Naive: Richness")
ggsave(file="output/Fig5B_part2.png", p)

# Evaluations of evenness removed two outliers per Grubs tests (see manuscript for details)
i <- which(tcrIpiProg$Patient == "pt103")
tcrIpiProg <- tcrIpiProg[-i,]
p<-tTestandGraph(tcrIpiProg, "Response", c("Benefit","NoBenefit"),"evenness.cdr3_On.Pre", c("Benefit", "No Benefit"), tit="Ipi Prog: Evenness")
ggsave(file="output/Fig5B_part3.png", p)

i <- which(tcrIpiNaive$Patient == "pt28")
tcrIpiNaive <- tcrIpiNaive[-i,]
p<-tTestandGraph(tcrIpiNaive, "Response", c("Benefit","NoBenefit"),"evenness.cdr3_On.Pre", c("Benefit", "No Benefit"), tit="Ipi Naive: Evenness")
ggsave(file="output/Fig5B_part4.png", p)





### Figure 5C
#Run JS_aaCDR3perVJ_2.5.py for each only.productive.tsv file. This will generate a new file containing the count, Shannon entropy, and normalized Shannon entropy (i.e. evenness) of CDR3 sequences for each VJ cassette combination in a given patient sample. Note that the patient ID must be added to each output file name before running the script for another patient file. The summary statistics from these output files were used to generate Figures 5C

tcr_sum_data$diffEveness_perVJ <-  tcr_sum_data$Hcdr3PerVJ_Median_On - tcr_sum_data$Hcdr3PerVJ_Median_Pre 
tcrIpiProg <- subset(tcr_sum_data, Cohort == "NIV3-PROG")
tcrIpiNaive <- subset(tcr_sum_data, Cohort == "NIV3-NAIVE")

# RICHNESS per VJ
tc <- c("Patient","Cohort","NumCDR3perVJ_Median_Pre",	"NumCDR3perVJ_Median_On")
tmpDat <- melt(tcrIpiProg[,tc], id=c("Patient","Cohort"))
p <- ggplot(tmpDat, aes(x=variable, y=value)) + geom_point() + geom_line(aes(group=Patient))+ggtitle("Ipi-Prog: Richness per VJ")
ggsave(file="output/Fig5C_part1.png", p)

tmpDat <- melt(tcrIpiNaive[,tc], id=c("Patient","Cohort"))
p <- ggplot(tmpDat, aes(x=variable, y=value)) + geom_point() + geom_line(aes(group=Patient))+ggtitle("Ipi-Naive: Richness per VJ")
ggsave(file="output/Fig5C_part2.png", p)


# Evenness per VJ
tc <- c("Patient","Cohort", "Hcdr3PerVJ_Median_On", "Hcdr3PerVJ_Median_Pre")
tmpDat <- melt(tcrIpiProg[,tc], id=c("Patient","Cohort"))
p <- ggplot(tmpDat, aes(x=variable, y=value)) + geom_point() + geom_line(aes(group=Patient))+ggtitle("Ipi-Prog: Evenness per VJ")
ggsave(file="output/Fig5C_part3.png", p)


tmpDat <- melt(tcrIpiNaive[,tc], id=c("Patient","Cohort"))
p <- ggplot(tmpDat, aes(x=variable, y=value)) + geom_point() + geom_line(aes(group=Patient))+ggtitle("Ipi-Naive: Evenness per VJ")
ggsave(file="output/Fig5C_part4.png", p)

