library(GenomicRanges)
library(Hmisc)
library(sequenza)
library(bootstrap)
library(boot)

###Load data
load('Data/OneDrive_12Jul/mutations.RData')
load('Data/OneDrive_12Jul/cnaSegments.RData')
load('Data/OneDrive_12Jul/purityRes.RData')
cnData <- res

output_directory <- 'Results/OneDrive_12Jul'

source('Scripts/ccfFunctions.R')

###Add MutId
muts$mutId <- paste(muts$id, seqnames(muts), start(muts), muts$ref, muts$alt, sep='_')

###Add copy number to mutations
muts$minorCN <- NA
muts$majorCN <- NA
muts$totalCN <- NA
muts$purity <- NA
muts$ploidy <- NA

###Add result columns
muts$vaf_lower <- NA
muts$vaf_upper <- NA
muts$mult_point <- NA
muts$mult_lower <- NA
muts$mult_upper <- NA
muts$ccf_point <- NA
muts$ccf_lower <- NA
muts$ccf_upper <- NA
muts$prob_clonal <- NA
muts$prob_subclonal <- NA
muts$timing <- NA

muts <- muts[!seqnames(muts)%in%c('chrY')]
###Add copy number data
for (id in unique(muts$id)){
	if (id %in% cnas$id){
	print(id)
	theseMuts <- muts[muts$id==id]
	theseCnas <- cnas[cnas$id==id]
	
	ovlp <- findOverlaps(theseMuts, theseCnas)
	muts$totalCN[match(theseMuts$mutId[queryHits(ovlp)], muts$mutId)] <- theseCnas$tumourTotal[subjectHits(ovlp)]
	muts$majorCN[match(theseMuts$mutId[queryHits(ovlp)], muts$mutId)] <- theseCnas$tumourMajor[subjectHits(ovlp)]
	muts$minorCN[match(theseMuts$mutId[queryHits(ovlp)], muts$mutId)] <- theseCnas$tumourMinor[subjectHits(ovlp)]
	
	###problematic
	outOfRange <- which(theseMuts$totalCN==0)

	muts$purity[muts$id==id] <- cnData$purity[match(id, cnData$id)]
	muts$ploidy[muts$id==id] <- cnData$ploidy[match(id, cnData$id)]
	
	#res <- earlyLate(theseMuts)
 	#muts$timing[match(theseMuts$mutId, muts$mutId)] <- res
	}
}



###Compute CCFs and things
for (i in 1:length(muts)){
	print(round(i/length(muts),2))
	res <- computeCCF(vaf=muts$vaf[i], depth=muts$depth[i], tumCN=muts$totalCN[i], purity=muts$purity[i])
	muts$vaf_lower[i] <- res['vaf_CI_lower']
	muts$vaf_upper[i] <- res['vaf_CI_upper']
	muts$mult_point[i] <- res['mult_point']
	muts$mult_lower[i] <- res['mult_lower']
	muts$mult_upper[i] <- res['mult_upper']
	muts$ccf_point[i] <- res['ccf_point']
	muts$ccf_lower[i] <- res['ccf_lower']
	muts$ccf_upper[i] <- res['ccf_upper']
	muts$prob_clonal[i] <- res['prob_clonal']
	muts$prob_subclonal[i] <- res['prob_subclonal']
}

muts$clonalStatus <- 'SUBCLONAL'
muts$clonalStatus[muts$ccf_upper==1] <- 'CLONAL'


###WRITE TO FILE
save(muts,res, file=paste0(output_directory, '/ccfMuts.RData'))

res <- data.frame(chr=seqnames(muts), start=start(muts), end=end(muts), ref=muts$ref, alt=muts$alt, sample=muts$id, depth=muts$depth, vaf=muts$vaf, gene=muts$gene, class=muts$class, minorCN=muts$minorCN, majorCN=muts$majorCN, totalCN=muts$totalCN, vaf_lower=muts$vaf_lower, vaf_upper=muts$vaf_upper, mult_point=muts$mult_point, mult_lower=muts$mult_lower, mult_upper=muts$mult_upper, ccf_point=muts$ccf_point, ccf_lower=muts$ccf_lower, ccf_upper=muts$ccf_upper, prob_clonal=muts$prob_clonal, prob_subclonal=muts$prob_subclonal, samplePurity=muts$purity, samplePloidy=muts$ploidy, stringsAsFactors=F)
write.table(res, file=paste0(output_directory,'/ccfMuts.txt'), sep='\t', quote=F, row.names=F)

mathScore <- sapply(unique(muts$id), function(x) mad(muts$vaf[muts$id==x])/median(muts$vaf[muts$id==x]))
devScore <- sapply(unique(muts$id), function(x) sd(muts$ccf_point[muts$id==x], na.rm=T))
ithScores <- cbind(mathScore, devScore)
save(ithScores, file=paste0(output_directory,'/ithScores.RData'))
ith2 <- data.frame(sample=rownames(ithScores), mathScore=ithScores[,'mathScore'], devScore=ithScores[,'devScore'])
write.table(ith2, file=paste0(output_directory,'/ithScores.txt'), sep='\t', quote=F, row.names=F)
