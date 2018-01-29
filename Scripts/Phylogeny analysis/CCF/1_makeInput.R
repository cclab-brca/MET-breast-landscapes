library(GenomicRanges)

segments_input_dir <- 'Data/OneDrive_12Jul/ascat/segments'
segments_output_file <- 'Data/OneDrive_12Jul/cnaSegments.RData'
purity_input_dir <- 'Data/OneDrive_12Jul/ascat/cellularity'
purity_output_file <- 'Data/OneDrive_12Jul/purityRes.RData'
mut_input_file <- 'Data/OneDrive_12Jul/all-filtered.maf'
mut_output_file <- 'Data/OneDrive_12Jul/mutations.RData'

####Make copy number object
CN_files <- dir(segments_input_dir, pattern='txt', recursive=T)
patientNames <- sapply(strsplit(CN_files, '_segments.txt'), function(x) x[1])

allDat <- GRangesList()
for (i in 1:length(CN_files)){
	tmp <- read.table(paste0(segments_input_dir, '/',CN_files[i]), header=T, as.is=T, sep='\t')
	tmpRanges <- GRanges(seqnames=paste0('chr', tmp$chr), ranges=IRanges(start=tmp$startpos, end=tmp$endpos))
	tmpRanges$id <- patientNames[i]
	tmpRanges$tumourMajor <- tmp$nMajor
	tmpRanges$tumourMinor <- tmp$nMinor
	tmpRanges$tumourTotal <- tmpRanges$tumourMajor+tmpRanges$tumourMinor
	allDat[[i]] <- tmpRanges
}

cnas <- unlist(allDat)
save(cnas, file=segments_output_file)


####Make purity object
purity_files <- dir(purity_input_dir, pattern='txt', recursive=T)
patientNames <- sapply(strsplit(purity_files, '_cellularity.txt'), function(x) x[1])

res <- data.frame(id=patientNames, purity=1, ploidy=2, goodnessOfFit=NA)
for (i in 1:length(purity_files)){
	tmp <- read.table(paste0(purity_input_dir, '/', purity_files[i]), header=F, as.is=T)
	res$purity[i] <- tmp$V2
	res$ploidy[i] <- tmp$V3
}

save(res, file=purity_output_file)

####Make mutation object 
tmp <- read.table(mut_input_file, sep='\t', as.is=T, header=T, quote='')

muts <- GRanges(seqnames=paste0('chr', tmp$Chr), ranges=IRanges(start=tmp$Start, end=tmp$End))
muts$id <- tmp$sample
muts$ref <- tmp$Ref_Allele
muts$alt <- tmp$Tumor_Alt_Allele
muts$depth <- tmp$tumour.depth
muts$vaf <- tmp$tumour.vaf
muts$exon <- tmp$Exon
muts$gene <- tmp$Hugo_Symbol
muts$class <- tmp$Ensembl_Variant

save(muts, file=mut_output_file)













