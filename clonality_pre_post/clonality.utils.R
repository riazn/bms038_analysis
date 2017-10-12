#clinical annotation
clin = read.csv(paste(PARENTDIR,"data/bms038_clinical_data.csv", sep=""), sep = ',',  head = T)

#response class per sample
response = unique(data.frame(gsub('_.*', '', clin$Sample), clin$myBOR))
class = as.character(response[,2])
names(class) = response[,1]

#classify variants based on change between pre and on treatment samples
update.clonality = function(dd, threshold = .95) {
	dd$clonal = rep(0, nrow(dd))
	dd$clonal[is.na(dd$pre.ccf) | is.na(dd$on.ccf)] = NA
	dd$clonal[dd$on.ccf > threshold & dd$pre.ccf > threshold] = 'clonal both'
	dd$clonal[dd$pre.ccf == 0 & dd$on.ccf > threshold] = 'novel clonal in on'
	dd$clonal[dd$pre.ccf == 0 & dd$on.ccf <= threshold] = 'novel subclonal in on'
	dd$clonal[dd$on.ccf == 0 & dd$pre.ccf <= threshold] = 'lost subclonal pre2on'
	dd$clonal[dd$on.ccf == 0 & dd$pre.ccf > threshold] = 'lost clonal pre2on'
	dd$clonal[dd$clonal == 0] = 'subclonal both'
	dd$clonal[dd$clonal == 'subclonal both'& dd$on.ccf > threshold & dd$pre.ccf <= threshold] = 'increased pre2on'
	dd$clonal[dd$clonal == 'subclonal both'& dd$on.ccf <= threshold & dd$pre.ccf > threshold] = 'decreased pre2on'
	dd$clonal = as.factor(dd$clonal)

	dd$clonal2 = rep(NA, nrow(dd))
	dd$clonal2[dd$clonal == 'clonal both'] = 'no selection'
	dd$clonal2[dd$clonal == 'novel clonal in on'] = 'positive selection'
	dd$clonal2[dd$clonal == 'novel subclonal in on'] = 'positive selection'
	dd$clonal2[dd$clonal == 'lost subclonal pre2on'] = 'negative selection'
	dd$clonal2[dd$clonal == 'lost clonal pre2on'] = 'negative selection'
	dd$clonal2[dd$clonal == 'subclonal both'] = 'no selection'
	dd$clonal2[dd$clonal == 'increased pre2on'] = 'positive selection'
	dd$clonal2[dd$clonal == 'decreased pre2on'] = 'negative selection'
	dd
}


# auxiliary function to switch between classifications
update.ccf = function(dd, ccf.type, threshold = .95, ci = T) {
	if(ccf.type == 'absolute' ) {
		dd$on.ccf = dd$on.ccf.abs
		dd$pre.ccf = dd$pre.ccf.abs
	}
	else if (ccf.type == 'pyclone') {
		if(ci == T) {
			dd$on.ccf = dd$on.ccf.pyc.ci
			dd$pre.ccf = dd$pre.ccf.pyc.ci
		}
		else{
			dd$on.ccf = dd$on.ccf.pyc
			dd$pre.ccf = dd$pre.ccf.pyc
		}
	}
	dd = update.clonality(dd, threshold)
}


