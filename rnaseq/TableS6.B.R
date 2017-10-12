###### Analyze pre-treatment Exome data ######
#
# Code for
#   figure 1B, Supp Figure 1A
rm(list=ls())
# CHANGE This to directory downloaded to
PARENTDIR <- "/mnt/skimcs/Melanoma/BMS_038/forGithub/bms038_analysis/"
setwd(PARENTDIR)

# Need to install these packages if not available
library("DESeq2")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("ggplot2")
library(annotate)
library(org.Hs.eg.db)
library("genefilter")
library("Biobase")
library("matrixStats")
library("IHW")

# Gene Filter function
filterthegenes <- function(DeseqObjevt){
	chch=DeseqObjevt
	chch <- chch[, !(is.na(chch$Response)) ]
	colData(chch)=droplevels(colData(chch))
	chch <- estimateSizeFactors(chch)

	#Filter expression measure above A in at least k samples
	#Remove low read counts (gene must have 4 or more samples with more or equal to 2 read count)
	mat=assay(chch)
	f1 <- kOverA(k = 4, A = 2,na.rm = T)
	ffun <- filterfun(f1)
	wh1 <- genefilter(mat, ffun)
	byKoA=which(wh1 == "FALSE")
	byKoA=names(byKoA)

	if(length(byKoA) == 0){ mock = chch}
	if(length(byKoA) > 0){ mock = chch[-(which(row.names(chch) %in% byKoA)), ]}
	# Filter based on delta between last 1/3 and first 1/2 samples
	# the gene with (variances of delta) < 1 got filtered
	mock=estimateSizeFactors(mock)
	normcounts=assay(mock,normalized=T)
	dlta=normcounts[,(ncol(mock)/2+1):(ncol(mock))]-normcounts[,1:(ncol(mock)/2)]
	vars=log10(rowVars(dlta))
	names(vars)=row.names(dlta)
	#hist(vars,100)
	byVar=unique(names(which(vars < 0)))

	# Filter based on read count variances
	# the gene with (variances of expression) < 1 or standard deviation < 1 got filtered
	mock <- chch[-(which(row.names(chch) %in% c(byKoA,byVar))), ]
	mock=estimateSizeFactors(mock)
	mat=assay(mock,normalized=T)
	vars=rowVars(mat)
	names(vars)=row.names(mat)
	madds=rowMads(mat)
	names(madds)=row.names(mat)
	#hist(log10(vars),100)
	#hist(log10(madds),100,main="Histogram of mads")
	byVar2=c(names(which(log10(vars) < 0)),names(which(log10(madds) <0)))

	myfilter=unique(c(byKoA,byVar,byVar2))
	return(myfilter)
}

# Gene symbol search function
ENTREZtoSYMBOL <- function(x){
  vec=mapIds(org.Hs.eg.db,keys=x,column="SYMBOL",keytype="ENTREZID",multiVals="first")
  repl=which(is.na(vec))
  for(y in repl){vec[y]=x[y]}
  return(vec)
}

# Result filter function: no NA & FDR < 0.2
ResFiltered <- function(x){
  temp=na.omit(x)
  temp=temp[which(temp$padj < 0.2),]
  temp$symbol=ENTREZtoSYMBOL(row.names(temp))
  return(temp)
}

##########################################################
# Build object and run DEG
##########################################################
# Read in matrix and Sample annotation
mat <- read.delim("data/CountData.BMS038.txt")
SampleTableCorrected <- read.csv("data/SampleTableCorrected.9.19.16.csv", row.names=1)

# Find the overlapping samples
inter <- intersect(colnames(mat),rownames(SampleTableCorrected))
SampleTableCorrected <- SampleTableCorrected[inter,]
mat <- mat[,match(rownames(SampleTableCorrected),colnames(mat))]

# Create the DESeq2 object
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
ebg <- exonsBy(txdb, by="gene")
intersection <- intersect(rownames(mat),ebg@partitioning@NAMES)
ebg2 <- ebg[ebg@partitioning@NAMES %in% intersection]

# sort by ID
ebg2 <- ebg2[order(names(ebg2))]

# Sort by gene model order
mat <- mat[match(names(ebg2), rownames(mat)),]

# Create object
ddsMat <- DESeqDataSetFromMatrix(countData = mat,
                                  colData = SampleTableCorrected,
                                  design = ~ 1)
# Assign genomic range
rowRanges(ddsMat) <- ebg2

# Get teh gene symbol
mcols(ddsMat)$symbols <-  mapIds(org.Hs.eg.db,
                                       keys=rownames(ddsMat),
                                       column="SYMBOL",
                                       keytype="ENTREZID",
                                       multiVals="first")

# Contraction and Expansion def
genomic <- read.table("data/survival.data.clonal.pyclone.0.95.TRUE.txt")
genomic$molicpheno <- ifelse(genomic$score_freq > 0, 'Expansion',ifelse(genomic$score_freq < 0,'Contraction','Neutral'))

# Find the overlapping samples
inter2 <- intersect(rownames(genomic),colnames(ddsMat))
genomic <- genomic[inter2,]
dds.genom <- ddsMat[,match(inter2,colnames(ddsMat))]

# Set the genomic phenotye levels
dds.genom$molicpheno <- factor(genomic$molicpheno, levels=c('Expansion','Contraction'))

# Drop samples "Pt48_Pre" "Pt52_Pre" who are not included
dds.genom <- dds.genom[,-which(colnames(dds.genom) %in% c("Pt48_Pre","Pt52_Pre"))]

# Filter out genes
gene.excluded <- read.csv('data/TableS6B_GeneExcluded.csv',row.names=1)
myfilter <- filterthegenes(dds.genom)
dds.genom <- dds.genom[ -which(rownames(dds.genom) %in% myfilter), ]
dds.genom <- dds.genom[ -which(rownames(dds.genom) %in% rownames(gene.excluded)), ]

# Size factor
dds.genom <- estimateSizeFactors(dds.genom)

# Set the design with controling cohort
design(dds.genom) <- ~ Cohort + molicpheno

# Run the DEG analysis
dds.genom <- DESeq(dds.genom)

# Get the Expansion (Persistence) vs. Contraction result
Response.result.genom <- ResFiltered(results(dds.genom, contrast=c("molicpheno","Expansion","Contraction")))
Response.result.genom <- Response.result.genom[order(Response.result.genom$padj),]

# Save result
write.table(Response.result.genom , 'output/TableS6.B.csv', sep=",", quote=F, col.names=NA)

