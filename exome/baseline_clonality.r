###### Analyze pre-treatment Exome data ######
#
# Code for
#   figure 1B, Supp Figure 1A
rm(list=ls())
# CHANGE This to directory downloaded to
PARENTDIR <- "/mnt/skimcs/Melanoma/BMS_038/forGithub/bms038_analysis/"
setwd(PARENTDIR)

# need to run baseline_analysis.r first to generate precursor files
# Need to install package nrstats in packages directory 
# R CMD install packages/nrstats-0.1.0.tgz
######################################################
library(nrstats)
library(ggplot2)
library(grid)
library(sqldf)



#################### FUNCTIONS ######################

readPyCloneDat <- function(dirLoc) {
  fls <- dir(dirLoc, "*ci95.tsv");
  #print(sprintf("here: reading %s", fls))
  ids <- c();
  median_ccfs <- c();
  muts <- c();
  clusters <- c();
  clusters_2muts <- c();
  thresh95 <- c()
  thresh90 <- c()
  thresh85 <- c()

  for (f in fls) {
    # Read in Cellular fraciton data for sample and determine median CCF for patient 
    print(f)
    f2<-paste(dirLoc,f,sep="")
    dat<-read.csv(f2, stringsAsFactors=FALSE, sep="\t")
    
    mc <- median(dat$cellular_prevalence);
    median_ccfs <- c(mc, median_ccfs);
    
    # Parse out patient id from file name
    id <- strsplit(f, "\\.")[[1]][1]
    ids <- c(id, ids)
    
    # Number of clusters (raw)
    num_clus <- length(unique(dat$cluster_id))
    clusters <- c(num_clus, clusters);
    
    
    # Number of clusters with at least 2 mutations
    tmp<-sqldf("select COUNT(*) as res from dat group by cluster_id")
    num_clus2 <- sum(tmp > 1);
    clusters_2muts <- c(num_clus2, clusters_2muts);
    
    # Fraction Clonal/Subclonal
    n <- dim(dat)[1]
    muts <- c(n, muts)
    thresh95 <- c(sum(dat$CI95 >= 0.95)/n, thresh95)
    thresh90 <- c(sum(dat$CI95 >= 0.90)/n, thresh90)
    thresh85 <- c(sum(dat$CI95 >= 0.85)/n, thresh85)
    
  }
  
  ccfdat <- data.frame(ids=ids, ccf=median_ccfs, clusters=clusters, clusters_2muts = clusters_2muts,
                       thresh95=thresh95, thresh90=thresh90, 
                       thresh85=thresh85, pyclone_muts=muts, stringsAsFactors=FALSE);
  rownames(ccfdat)<-ccfdat$ids
  return(ccfdat);
}

######################################################

# Create clonal mutation load data from pyclone output
indir <- paste(PARENTDIR, "data/pyclone/pre/", sep="");
dat <- readPyCloneDat(indir)

# read in patient clinicl data
ptdat <- read.csv("data/bms038_clinical_data.csv",row.names=1, stringsAsFactors = FALSE)
mutdat <- read.csv("output/MutCnt_pretheray.csv")
combdat <- merge(ptdat, dat, by.x="Sample", by.y="ids")
combdat <- merge(combdat, mutdat, by.x="Sample", by.y="Sample")
naive_i <- which(combdat$Cohort=="NIV3-NAIVE")
prog_i <- which(combdat$Cohort=="NIV3-PROG")


combdat$thresh95.muts <- log10(combdat$NonSynMut) * combdat$thresh95
combdat$thresh90.muts <- log10(combdat$NonSynMut) * combdat$thresh90
combdat$thresh85.muts <- log10(combdat$NonSynMut) * combdat$thresh85




# Check clonal mutation load and survival (cut point at median)
# Overall Survival by Clonal Mutation load in Ipi-Naive Patients (Figure 1B)
png("output/Fig1B_part2.png")
med_val <- median(combdat$thresh95.muts[naive_i], na.rm=TRUE)
combdat$clon_cut[naive_i] <- combdat$thresh95.muts[naive_i] >= med_val
with(combdat[naive_i,], qkm2(OSWK/52, !OS_SOR,3, clon_cut, c("low","high"),"Clonal Mut Load: OS Naive Pts"))
dev.off()

png("output/SupFig1D.png")
med_val <- median(combdat$thresh95.muts[prog_i], na.rm=TRUE)
combdat$clon_cut[prog_i] <- combdat$thresh95.muts[prog_i] >= med_val
with(combdat[prog_i,], qkm2(OSWK/52, !OS_SOR,3, clon_cut, c("low","high"),"Clonal Mut Load: OS Prog Pts"))
dev.off()
######################################################



p<-tTestandGraph(combdat, "Cohort", c("NIV3-NAIVE","NIV3-PROG"), "thresh95",c("Ipi-Naive","Ipi-Prog"))
ggsave(file="output/SupFig1b.png", p)

p<-tTestandGraph(combdat[naive_i,], "myBOR", c("PRCR","PD"), "thresh95.muts", c("PR/CR","PD"), tit="Ipi-Naive")
ggsave(file="output/SupFig1c_part1.png", p)

p<-tTestandGraph(combdat[prog_i,], "myBOR", c("PRCR", "PD"), "thresh95.muts", c("PR/CR","PD"), tit="Ipi-Prog")
ggsave(file="output/SupFig1c_part2.png", p)

