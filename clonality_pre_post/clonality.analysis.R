###### Analyze Pre and On-therapy clonality changes ######
#
# Code for
#   Code to generate Fig 2a, 2b, 2d

rm(list=ls())
# CHANGE This to directory downloaded to
PARENTDIR <- "/mnt/skimcs/Melanoma/BMS_038/forGithub/bms038_analysis/"
setwd(PARENTDIR)

# Need to install package nrstats in packages directory 
# R CMD install packages/nrstats-0.1.0.tgz
######################################################
require(RColorBrewer)
require(survival)
require(ks)

source(paste(PARENTDIR,'clonality_pre_post/clonality.utils.R', sep="/"))



dd = read.table(paste(PARENTDIR,'/data/clonality.data.table.pyclone.95.ci.txt', sep=""), sep = '\t')

sample.pairs = unique(as.character(dd$sample))


reds = brewer.pal(6, 'PuRd')
greens = brewer.pal(6, 'Blues')
yellows = brewer.pal(7, 'Greys')

cc2 = c(reds[4], yellows[3], greens[4])
cc1 = c(reds[4:2], yellows[3:2], greens[3:1])

type.order = c('novel clonal in on', 'novel subclonal in on', 'increased pre2on', 'clonal both', 'subclonal both', 'lost clonal pre2on', 'lost subclonal pre2on', 'decreased pre2on')

update.plot.data = function(dd, plot.type) {
	if(plot.type == 'clonal') {
		data = sapply(split(as.factor(dd$clonal), dd$sample), table)
		data = data[type.order, idx]
	}
	else if(plot.type == 'selection') {
		data = sapply(split(as.factor(dd$clonal2), dd$sample), table)
	}
	data
}


# barplots to vizualize the different changes in CCF between pre and on treatment samples
make.barplot = function(dd, plot.type, ccf.type, threshold, ci = F, relative = F) {
	dd = update.ccf(dd, ccf.type, threshold, ci)
	type.order = c('novel clonal in on', 'novel subclonal in on', 'increased pre2on', 'clonal both', 'subclonal both', 'lost clonal pre2on', 'lost subclonal pre2on', 'decreased pre2on')
	samples = as.character(unique(dd$sample))
	sample.count = table(dd$sample)[samples]
	samples = paste(class[samples], samples)
	idx = order(substr(samples, 1, 2), sample.count)

	if(plot.type == 'clonal') {
		data = sapply(split(as.factor(dd$clonal), dd$sample), table)
		data = data[type.order, idx]
		if(relative) data = t(t(data)/colSums(data))
		barplot(data, beside = F, col = cc1, las = 2, border = NA)
		legend('topleft', legend = rev(rownames(data[type.order,])), fill = rev(cc1))
	}
	if(plot.type == 'selection') {
		data = sapply(split(as.factor(dd$clonal2), dd$sample), table)
		if(relative) data = t(t(data)/colSums(data))
		barplot(data[c('positive selection', 'no selection', 'negative selection'), idx], beside = F, col = cc2, las = 2, border = NA)
		legend('topright', legend = rev(rownames(data[c('positive selection', 'no selection', 'negative selection'),])), fill = rev(cc2))
	}
    return(data[,idx])
}

print('Barplots (Fig 2A) ..')
pdf(paste(PARENTDIR,'/output/Fig2A_part1_pyclone.95.clonal.ci.pdf', sep="")); make.barplot(dd, plot.type = 'clonal', ccf.type = 'pyclone', threshold = .95, ci = T); dev.off();
pdf(paste(PARENTDIR,'/output/Fig2A_part2_pyclone.95.clonal.relative.ci.pdf', sep="")); make.barplot(dd, plot.type = 'clonal', ccf.type = 'pyclone', threshold = .95, ci = T, relative = T); dev.off();
pdf(paste(PARENTDIR,'/output/Fig2A_part3_pyclone.95.selection.relative.ci.pdf', sep="")); tmp = make.barplot(dd, plot.type = 'selection', ccf.type = 'pyclone', threshold = .95, ci = T, relative = T); dev.off();


# boxplots to compare sets of variants (contraction/persistence/expansion)
make.boxplot = function(dd, plot.type, ccf.type, threshold, ci = F) {
	dd = update.ccf(dd, ccf.type, threshold, ci)
	type.order = c('novel clonal in on', 'novel subclonal in on', 'increased pre2on', 'clonal both', 'subclonal both', 'lost clonal pre2on', 'lost subclonal pre2on', 'decreased pre2on')
	samples = as.character(unique(dd$sample))
	sample.count = table(dd$sample)[samples]
	samples = paste(class[samples], samples)
	idx = order(substr(samples, 1, 2), sample.count)

	data = sapply(split(as.factor(dd$clonal2), dd$sample), table)
	data = t(t(data)/colSums(data))
	colnames(data) = paste(class[colnames(data)], colnames(data))
	cl = gsub(' .*', '', colnames(data))
	sapply(rownames(data), function(x) {
		       y = split(data[x,], cl); 
		       pv1 = round(t.test(y$PD, y$SD, alternative = 'greater')$p.value, digits = 5); 
		       pv2 = round(t.test(y$PD, y$SD, alternative = 'less')$p.value, digits = 5); 
		       boxplot(list(PRCR = y$PRCR, SD = y$SD, PD = y$PD), border = c('darkgreen', 'orange', 'red'), outline = F, ylim = c(0,1), main = paste(x, ', p.vals = ', pv1, pv2))
	})
}


print('Boxplots (Fig 2B) ..')
pdf(paste(PARENTDIR,'/output/Fig2B_pyclone.95.ci.boxplot.pdf', sep="")); 
make.boxplot(dd, plot.type = 'clonal', ccf.type = 'pyclone', threshold = .95, ci = T); 
dev.off();


# density plots for each pair of pre and on treatment samples
plot.density = function(dd, sample, ccf.type, threshold, ci = F) {
	dd = update.ccf(dd, ccf.type, threshold, ci)
	idx = as.character(dd$sample) == sample 
	data = cbind(dd$pre.ccf[idx], dd$on.ccf[idx]) 
	data = data[rowSums(is.na(data)) == 0,]
	data[,1] = jitter(data[,1])
	data[,2] = jitter(data[,2])
	fhat = kde(x=as.data.frame(data));
	plot(fhat, main = paste(class[sample], sample), xlim=c(-.3, 1.3), ylim = c(-.3, 1.3), display='filled.contour2', cont = seq(0,100,1),xlab = 'pre ccf', ylab = 'on ccf')
	points(data, cex=.2, pch=16)
}


print('Density plots (Fig 2D)..')
pdf(paste(PARENTDIR,'/output/Fig2D_pyclone.95.densities.pdf', sep="")); 
sapply(sample.pairs, function (sample) {
	plot.density(dd, sample, 'pyclone', threshold = .95, ci = F) 
})
dev.off()


