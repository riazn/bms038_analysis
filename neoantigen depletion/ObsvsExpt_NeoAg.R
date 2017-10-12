### Assessment of on-therapy neoantigen depletion ###

### See README file for details

rm(list = ls())

dat = read.table("Mut_NeoAg_Expected_vs_Observed_BMS038.csv",header=T,sep=",")

matriz <- matrix(c(log10(dat$NeoAgCnt_On),log10(dat$E_NeoAgCnt_On)), nrow=nrow(dat), ncol=2)

rownames(matriz) <- dat$Patient

matriz=t(matriz)
barplot(matriz, beside=TRUE, legend=TRUE, col=c("gray","blue"),las=2,names.arg=dat$Patient, ylab="# on-therapy neoantigens (log10)",ylim=c(0,3.5),cex.names=1, cex.lab=1.5, cex.axis = 1.2)
legend(2,3,legend=c("Observed", "Expected"),bty = "n", fill= c("gray", "blue"), cex = 1.2)
