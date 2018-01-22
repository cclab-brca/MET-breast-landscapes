rawdatadir <- '~/Documents/Caldas/Data/dataRaw/SNVs_WES/'
treeomicsDir_WES <- '~/Documents/Caldas/Data/dataProcessed/treeomics/output_WES_all/'
outdatadir <- '~/Documents/Caldas/Data/dataProcessed/treeomics/input_WES_selected/'

dir.create(outdatadir,showWarnings = FALSE)

library("data.table")

cases <- c(288,"290exOvary",291,"298Brain","298main",308,315,323,328,330,'DET52')

for (case in cases) {
  if (case=='298Brain'|case=='298main') {
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
  
  ## read posterior probabilities from previous treeomics run (used for filtering)
  postFile <- list.files(path=paste0(treeomicsDir_WES,'/X',case,'/'),pattern='posterior.txt',full.names = TRUE)
  posterior <- read.csv(postFile,sep='\t',stringsAsFactors = FALSE)
  ## remove rows (purities, VAFs and priors)
  posterior <- posterior[-(1:4),]
  colnames(posterior)[1:4] <- c('Chromosome','Position','Position2','Gene.Name')
  
  ## match rows between posterior, covCount and altCount
  posterior <- as.data.table(posterior)
  setkey(posterior,Chromosome,Position)
  posterior$ID <- paste0(posterior$Chromosome,'_', posterior$Position)
  posterior$ID2 <- paste0(posterior$Chromosome,'_', posterior$Position2)
  
  altCount <- as.data.table(altCount)
  setkey(altCount,Chromosome,Position)
  altCount$ID <- paste0(altCount$Chromosome,'_', altCount$Position)
  altCount <- altCount[altCount$ID%in%posterior$ID | altCount$ID%in%posterior$ID2,]
  altCount[,ID := NULL]
  
  covCount <- as.data.table(covCount)
  setkey(covCount,Chromosome,Position)
  covCount$ID <- paste0(covCount$Chromosome,'_', covCount$Position)
  covCount <- covCount[covCount$ID%in%posterior$ID | covCount$ID%in%posterior$ID2,]
  covCount[,ID := NULL]
  
  posterior[,c("ID","ID2") := NULL]
  
  if (baseCase==290) {
    
    colnames(altCount)[-(1:4)] <- colnames(covCount)[-(1:4)] <- c("X004","X005","X007","X008","X014","X015","X016A","X016BIDC","X016BWT","X016Bmuc","X017","X018","X019","X021","X024")
    ## eclude samples to make data set smaller
    exclude <- colnames(altCount)%in%c('X016BIDC', 'X016BWT', 'X016Bmuc')
    altCount <- altCount[,!exclude,with=FALSE]
    covCount <- covCount[,!exclude,with=FALSE]
    posterior <- posterior[,which(!colnames(posterior)%in%c('X016BIDC', 'X016BWT', 'X016Bmuc')),with=FALSE]
  } else {
    colnames(altCount)[-(1:4)] <- colnames(covCount)[-(1:4)] <- paste0('X',gsub(paste0("X?",baseCase,"\\.([mt0-9|normal]*)\\.alt_count"),"\\1",colnames(altCount)[-(1:4)]))#paste0('X',gsub("X[0-9]*\\.([0-9]*)\\.alt_count","\\1",colnames(altCount)[-(1:4)]))
  }
  

  if (case==308) {
    exclude <- colnames(altCount)%in%paste0('X',c('007','009','015','020','022'))    
    altCount <- altCount[,!exclude,with=FALSE] 
    covCount <- covCount[,!exclude,with=FALSE]
    posterior <- posterior[,!colnames(posterior)%in%paste0('X',c('007','009','015','020','022'))  ,with=FALSE]
  }
  
  
  if (case=="DET52") {
    colnames(altCount)[-(1:4)] <- colnames(covCount)[-(1:4)] <- colnames(posterior)[-(1:4)]
    
    idx <- grep("Xt",colnames(altCount))
    altCountA <- altCount[,-idx,with=FALSE] 
    covCountA <- covCount[,-idx,with=FALSE]
    posteriorA <- posterior[,-grep("Xt",colnames(posterior)),with=FALSE]
    
    ## FILTER MUTATIONS
    ## number of mutated samples
    Nmut <- apply(posteriorA,1,function(x) sum(x[-(1:4)]>0.5))
    ## remove sites where no sample has a probability of being mutated above 0.5
    altCountA <- altCountA[Nmut>=1,]
    covCountA <- covCountA[Nmut>=1,]
    
    write.table(covCountA,paste0(outdatadir,'/X',case,'main_covCount.txt'),sep='\t',row.names=FALSE,quote=FALSE,na = "n/a")
    write.table(altCountA,paste0(outdatadir,'/X',case,'main_altCount.txt'),sep='\t',row.names=FALSE,quote=FALSE,na = "n/a")
  }
  
  ## FILTER MUTATIONS
  ## number of mutated samples
  Nmut <- apply(posterior,1,function(x) sum(x[-(1:4)]>0.5))
  ## remove sites where no sample has a probability of being mutated above 0.5
  altCount <- altCount[Nmut>=1,]
  covCount <- covCount[Nmut>=1,]
  posterior <- posterior[Nmut>=1,]

  ## write to file
  write.table(covCount,paste0(outdatadir,'/X',case,'_covCount.txt'),sep='\t',row.names=FALSE,quote=FALSE,na = "n/a")
  write.table(altCount,paste0(outdatadir,'/X',case,'_altCount.txt'),sep='\t',row.names=FALSE,quote=FALSE,na = "n/a")
}