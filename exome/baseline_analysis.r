###### Analyze pre-treatment Exome data ######
#
# Code for
#   figure 1B, Supp Figure 1A
rm(list=ls())
# CHANGE This to directory downloaded to
PARENTDIR <- "/mnt/skimcs/Melanoma/BMS_038/forGithub/bms038_analysis/"
setwd(PARENTDIR)

# Need to install package nrstats in packages directory 
# R CMD install packages/nrstats-0.1.0.tgz
######################################################
library(nrstats)
library(ggplot2)
library(grid)
library(sqldf)

##### CONSTANTS
ExomeN <- 68



####### READ IN Patient Data #######
# read in patient data 
ptdat <- read.csv("data/bms038_clinical_data.csv",row.names=1, stringsAsFactors = FALSE)
mutdat <- read.csv("data/pre_therapy_nonsynonmous_mutations.csv", stringsAsFactors = FALSE)
genomic_data_samp <- read.csv("data/genomic_data_per_case.csv", stringsAsFactors = FALSE)
rownames(genomic_data_samp) <- genomic_data_samp$Patient


# Determine non-synonmous mutation load per case from mutation
mutcntdat <-sqldf("select Patient, COUNT(*) as NonSynMut from mutdat GROUP by Patient")
rownames(mutcntdat) <- mutcntdat$Patient
mutcntdat$Sample <- paste(mutcntdat$Patient, "_pre", sep="")
write.csv(mutcntdat, file="output/MutCnt_pretheray.csv")

ptdat <- merge(ptdat, genomic_data_samp, by.x=0, by.y=0, all.x=TRUE)
rownames(ptdat) <- ptdat$Patient
ptdat <- merge(ptdat, mutcntdat, by.x=0, by.y=0, all.x=TRUE)
ptdat$Patient <- ptdat$Patient.x
rownames(ptdat) <- ptdat$Patient

# For analysis here, only keep patient that have pre-therapy exome data
ptdat <- subset(ptdat, Pre.treatment.Exome == 1)


####### Simple statistics in results section  #######
# Display simple stats on IQL; min, max (in results)
summary(ptdat$NonSynMut)


# Frequency of Subtypes
table(ptdat$SubtypeEZ) / ExomeN

# Compute p-value of TripleWt in Ipi-Naive vs. Ipi-Prog
tbl1<-table(ptdat$SubtypeEZ, ptdat$Cohort)
m <- matrix(NA,2,2)
m[1,] <- c(11+2+9, 8+2+5)
m[2,] <- c(11,20)
cat(sprintf("P-value between TripleWt in ipi-naive vs. ipi-prog %f\n", fisher.test(m)$p.value))



# Comparison of mutload between ipi-Niaive and ipi-Prog
naive_i <- which(ptdat$Cohort=="NIV3-NAIVE")
prog_i <- which(ptdat$Cohort=="NIV3-PROG")
mn_mut_load_naive <- 
cat(sprintf("Median mutation load: Ipi-Naive: %.2f, Ipi-Prog %.2f", median(ptdat$NonSynMut[naive_i]), median(ptdat$NonSynMut[prog_i])))
pv_wilcox <- wilcox.test(ptdat$NonSynMut[naive_i], ptdat$NonSynMut[prog_i])$p.value
pv_ttest <- t.test(log10(ptdat$NonSynMut[naive_i]+1), log10(ptdat$NonSynMut[prog_i]+1))$p.value
cat(sprintf("Comparison of mutaiton loads between ipi-naive & prog: wilcox: p=%.2f, t-test p=%.2f", pv_wilcox,pv_ttest))


####### Survival analysis of Figure 1B, Sup Fig1A,D  #######

# Check mutation load and survival (using cutpoint from NEJM)
ptdat$mut_cut <- ptdat$NonSynMut > 100
with(ptdat, qkm2(OSWK/52, !OS_SOR,3, mut_cut, c("low","high"),"All Pts"))

# Overall Survival by Mutation load in Ipi-Naive Patients (Figure 1B)
png("output/Fig1B_part1.png")
with(ptdat[naive_i,], qkm2(OSWK/52, !OS_SOR,3, mut_cut, c("low","high"),"Mutation Load OS: Naive Pts"))
dev.off()
# Overall Survival by Mutation load in Ipi-Prog Patients (Sup Figure 1A)
png("output/SupFig1A.png")
with(ptdat[prog_i,], qkm2(OSWK/52, !OS_SOR,3, mut_cut, c("low","high"),"Mutation Load OS: Prog Pts"))
dev.off()

