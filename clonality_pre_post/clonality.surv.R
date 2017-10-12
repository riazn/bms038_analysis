###### Analyze Pre and On-therapy clonality changes ######
#
# Code for
#   Code to generate Fig 2c

rm(list=ls())
# CHANGE This to directory downloaded to
PARENTDIR <- "/mnt/skimcs/Melanoma/BMS_038/forGithub/bms038_analysis/"
setwd(PARENTDIR)

# Need to install package nrstats in packages directory 
# R CMD install packages/nrstats-0.1.0.tgz
######################################################
library(nrstats)
source(paste(PARENTDIR,'clonality_pre_post/clonality.utils.R', sep="/"))

dd = read.table(paste(PARENTDIR,'/data/clonality.data.table.pyclone.95.ci.txt',sep=""), sep = '\t')

qkm_var_quant <- function(event_time, event, max, var, tit) {
  # Create variable that splits value at median and plot survival curve
  quants = cut(var, breaks = quantile(var, probs = seq(0,1, by=.25)), include.lowest = T)
  qkm2(event_time, event, max, quants, c("q1", "q2", "q3","q4"), tit);
}

make.survival.plot = function(dd, plot.type, ccf.type, threshold, ci=F) {
	dd = update.ccf(dd, ccf.type, threshold, ci)
	dd = dd[!is.na(dd$pre.ccf) & !is.na(dd$on.ccf),]

	samples = as.character(unique(dd$sample))
	sample.count = table(dd$sample)[samples]
	samples = paste(class[samples], samples)
	idx = order(substr(samples, 1, 2), sample.count)

	data = sapply(split(as.factor(dd$clonal2), dd$sample), table)
	rel.data = t(t(data)/colSums(data))
	
	var2 = rel.data['no selection',] - rel.data['negative selection',]
	var3 = data['no selection',] - data['negative selection',]
	var5 = var2 > 0
	
	#rel and absolute data give same order -> same surviavl curves
	
	names(var2) = names(var3) = names(var5) = paste(colnames(data), '_pre', sep = '')
	oswk = clin$OSWK
	os_sor = clin$OS_SOR
	pfswk = clin$PFSWK
	pfs_sor = clin$PFS_SOR
	names(pfswk) = names(pfs_sor) = names(oswk) = names(os_sor) = as.character(clin$Sample)
	
	tmp = dd[dd$pre.ccf != 0,]
	pre.count = table(tmp$sample)
	tmp = dd[dd$on.ccf != 0,]
	on.count = table(tmp$sample)

	idx = gsub('_.*','', names(var5))
	dt = cbind(oswk[names(var5)]/52, !os_sor[names(var5)], var3, var2, data['positive selection', idx], data['no selection',idx], data['negative selection',idx], pre.count[idx], on.count[idx])
	colnames(dt) =c('oswk', 'os_sor', 'score_count', 'score_freq', 'new/increased', 'stable', 'lost/decreased', 'pre.count', 'on.count')
	dt = dt[order(dt[,4]),]
	write.table(dt, file = paste('survival.data',plot.type, ccf.type, threshold, ci, 'txt',sep = '.'), sep = '\t')
	pdf(paste(PARENTDIR,paste('/output/Fig2C_survival.data',plot.type, ccf.type, threshold, ci, 'pdf',sep = '.'), sep=""))
	qkm2(oswk[names(var5)]/52, !os_sor[names(var5)], 3, var5, tit = 'neg - pos, absolute freq, OS')
	qkm2(pfswk[names(var5)]/52, !pfs_sor[names(var5)], 3, var5, tit = 'neg - pos, absolute freq, PFS')
	vec = var2; names(vec) = names(var2)
	vec = sort(vec) 
	barplot(vec, las = 2, col = c('red', 'darkgreen', 'orange')[as.factor(class[gsub('_.*', '', names(vec))])])
	dev.off()
}

print('Survival plots (Fig 2C)')
make.survival.plot(dd, plot.type = 'clonal', ccf.type = 'pyclone', threshold = .95, ci = T)





