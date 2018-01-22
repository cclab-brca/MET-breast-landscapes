rawdatadir <- '~/Documents/Caldas/Data/dataRaw/SNVs_WES/'
outdatadir <- '~/Documents/Caldas/Data/dataProcessed/treeomics/input_WES_all/'

dir.create(outdatadir,showWarnings = FALSE)

cases <- c(288,'290exOvary',291,'298Brain','298main','298All',308,315,323,328,330,'DET52')

for (case in cases) {
  
  if (case=='298Brain'|case=='298main'|case=='298All') {
    baseCase=298
  } else if (case=="290exOvary") {
    baseCase=290
  } else {
    baseCase=case
  }
  
  if (case=='298Brain') {
    SNPfile <- list.files(path=rawdatadir,
                          pattern='298-023.filtered.maf',
                          full.names = TRUE)
  } else if (case=='298main') {
    SNPfile <- list.files(path=rawdatadir,
                          pattern='298-023en.filtered.maf',
                          full.names = TRUE)
  } else if (case=='298All') {
    SNPfile <- list.files(path=rawdatadir,
                          pattern='298.filtered.maf',
                          full.names = TRUE)
  } else {
    SNPfile <- list.files(path=rawdatadir,
                          pattern=paste0(baseCase,"-[0-9|normal]*.filtered.maf"),
                          full.names = TRUE) 
  }
  dat <- read.table(SNPfile,sep = '\t',header=TRUE,stringsAsFactors = FALSE,quote="")
  
  idxAlt <- grep('alt_count',colnames(dat))
  idxDepth <- grep('depth',colnames(dat))
  
  metadata <- data.frame(Chromosome=dat$Chr,
                         Position=dat$Start,
                         Change=paste0(dat$Ref_Allele,'>',dat$Tumor_Alt_Allele),
                         Gene=dat$Hugo_Symbol)
  
  altCount <- data.frame(metadata,dat[,idxAlt])
  covCount <- data.frame(metadata,dat[,idxDepth])
  
  
  if (baseCase==290) {  
    colnames(altCount)[-(1:4)] <- colnames(covCount)[-(1:4)] <- c("X004","X005","X007","X008","X014","X015","X016A","X016BIDC","X016BWT","X016Bmuc","X017","X018","X019","X021","X024")
  } else {
    colnames(altCount)[-(1:4)] <- colnames(covCount)[-(1:4)] <- paste0('X',gsub(paste0("X?",baseCase,"\\.([mt0-9|normal]*)\\.alt_count"),"\\1",colnames(altCount)[-(1:4)]))#paste0('X',gsub("X[0-9]*\\.([0-9]*)\\.alt_count","\\1",colnames(altCount)[-(1:4)]))
  }
  
  if (case=="DET52") {
    idx <- grep("Xt",colnames(altCount))
    altCountA <- altCount[,-idx] 
    covCountA <- covCount[,-idx]
    write.table(covCountA,paste0(outdatadir,'/X',case,'main_covCount.txt'),sep='\t',row.names=FALSE,quote=FALSE,na = "n/a")
    write.table(altCountA,paste0(outdatadir,'/X',case,'main_altCount.txt'),sep='\t',row.names=FALSE,quote=FALSE,na = "n/a")
  }

  ## write to file
  write.table(covCount,paste0(outdatadir,'/X',case,'_covCount.txt'),sep='\t',row.names=FALSE,quote=FALSE,na = "n/a")
  write.table(altCount,paste0(outdatadir,'/X',case,'_altCount.txt'),sep='\t',row.names=FALSE,quote=FALSE,na = "n/a")
}