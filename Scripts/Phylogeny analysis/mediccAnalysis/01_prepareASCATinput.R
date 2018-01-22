calcBAF=TRUE
calcLogR=TRUE
merge=TRUE

##First read in the arguments listed at the command line
args=(commandArgs(TRUE))

##args is now a list of character vectors
## First check to see if arguments are passed.
if(length(args)==0){
  stop("No arguments supplied.")
} else{
  baseCase = args[1]
  case = args[2]
  samplesFile = args[3]
  normal_sample = args[4]
  dataDir = args[5]
  QDNAseqDir = args[6]
  SNP_pos_with_names = args[7]
  count_files = args[8]
  normal_count_file = args[9]
  chrCount = as.numeric(args[10])
}

nonGenderChrs = chrCount - 2

################################################################################
#### load packages ####

library(data.table)
library(GenomicRanges)

################################################################################
#### create dirs and load sample information ####

ascatInDir = file.path(dataDir,paste0('ascatInput/',case))
ascatOutDir = file.path(dataDir,paste0('ascatOutput/',case))
dir.create(ascatInDir,recursive = TRUE,showWarnings = FALSE)
dir.create(ascatOutDir,recursive = TRUE,showWarnings = FALSE)

## read in sample information
samples <- read.table(samplesFile,stringsAsFactors = FALSE)$V1

################################################################################
#### Step 1: calculate BAFs ####

if (calcBAF) {
  
  ### modified from runASCAT.R (part of ascatNGS pipeline)
  normalcounts = read.table(normal_count_file,sep="\t")
  
  SNPpos = matrix(nrow = dim(normalcounts)[1],ncol = 2)
  rownames(SNPpos) = paste("snp",1:dim(SNPpos)[1],sep="")
  colnames(SNPpos) = c("Chr","Position")
  SNPpos[,1] = as.vector(normalcounts[,1])
  SNPpos[,2] = normalcounts[,2]
  
  ## This gives every SNP a name (to be able to do GC correction)
  SNPposWithNames = read.table(SNP_pos_with_names,sep="\t",header=T,row.names=1)
  
  ctrans = 1:chrCount
  names(ctrans)=c(1:nonGenderChrs,"X","Y")
  newnames = ctrans[as.vector(SNPposWithNames[,1])]*1000000000+SNPposWithNames[,2]
  newnamesSEQ = ctrans[as.vector(SNPpos[,1])]*1000000000+as.numeric(SNPpos[,2])
  
  namestrans = rownames(SNPposWithNames)
  names(namestrans)=newnames
  
  rownames(SNPpos) = newnamesSEQ
  SNPpos = SNPpos[rownames(SNPpos)%in%newnames,]
  rownames(SNPpos)=namestrans[rownames(SNPpos)]
  
  rownames(normalcounts) = newnamesSEQ
  normalcounts = normalcounts[rownames(normalcounts)%in%newnames,]
  rownames(normalcounts)=namestrans[rownames(normalcounts)]
  
  for (sampleID in samples) {
    count_file <- paste0(eval(parse(text=count_files)),collapse="")
    counts = read.table(count_file,sep="\t")
    
    rownames(counts) = newnamesSEQ
    counts = counts[rownames(counts)%in%newnames,]
    rownames(counts)=namestrans[rownames(counts)]
    
    BAF = matrix(nrow = dim(normalcounts)[1],ncol = 1)
    rownames(BAF) = rownames(SNPpos)
    colnames(BAF) = sampleID
    acgt = counts[,c(3:6)]
    acgts = t(apply(acgt,1,sort))
    BAF[,1] = acgts[,4]/(acgts[,3]+acgts[,4])
    BAF[,1] = ifelse(runif(length(BAF[,1]))<0.5,BAF[,1],1-BAF[,1])
    BAF[is.nan(BAF)]=NA
    
    # limit the number of digits:
    BAF = round(BAF,4)
    
    # write output to files
    write.table(cbind(SNPpos,BAF),file.path(ascatInDir,paste0(sampleID,".BAF.txt")),sep="\t",row.names=T,col.names=NA,quote=F)
  }
}

################################################################################
#### merge BAFs for all samples ####

if (merge) {
  
  #### BAF #### 
  dat <- fread(file.path(ascatInDir,paste0(samples[1],'.BAF.txt')),colClasses = c('character','character','numeric','numeric'))
  setnames(dat, c(colnames(dat)[-ncol(dat)],samples[1]))
  
  for (sampleID in samples[-1]) {
    
    dats <- fread(file.path(ascatInDir,paste0(sampleID,'.BAF.txt')),colClasses = c('character','character','numeric','numeric'))
    setnames(dats, paste0('V',1:ncol(dats)))
    dat[,(paste0(sampleID)) := dats[,V4]]
    
  }
  
  write(paste0(c('','Chr','Position',samples),collapse = '\t'),file=file.path(ascatInDir,'merged.tumour.BAF.txt'),sep = '\t')
  write.table(dat,file=file.path(ascatInDir,'merged.tumour.BAF.txt'),append=TRUE,quote=FALSE,sep='\t',row.names = FALSE,col.names = FALSE)
  
}

################################################################################
#### Step 2: prepare LogR files from QDNAseq output ####

if (calcLogR) {
  
  ## logr fomr QDNAseq
  ## read in SNP positions
  pos <- fread(SNP_pos_with_names,colClasses = c('character','character','numeric'),skip=1)
  pos[,V4:=V3]
  setnames(pos,c('name','chromosome','start','end'))
  posRanges <- makeGRangesFromDataFrame(pos)
  
  ## select ratio values from QDNAseq data
  logr <- fread(paste0(QDNAseqDir,'/',baseCase,'_QDNAseq.txt'),nrows = 4)
  logr <- fread(paste0(QDNAseqDir,'/',baseCase,'_QDNAseq.txt'),colClasses = c("character","character",rep("numeric",ncol(logr)-2)))
  idx <- c(1:4,grep('ratio',colnames(logr)))
  logr <- logr[,idx,with=FALSE]
  setnames(logr,gsub('.ratio','',colnames(logr)))
  if (case=='DET52') {
    logr[,'DET52_normal':=rep(0.003,nrow(logr))]
    setnames(logr,gsub("_","-",colnames(logr)))
  }
  ## select subset of samples if current case excludes samples
  idx <- which(colnames(logr)%in%c('feature','chromosome','start','end',samples))
  logr <- logr[,idx,with=FALSE]
  
  logRanges <- makeGRangesFromDataFrame(logr,
                                        keep.extra.columns=TRUE)
  
  ## find overlap with SNP positions
  logRX <- as.data.table(matrix(as.numeric(NA),nrow=nrow(pos),ncol=ncol(logr)-4),colClasses='numeric')
  hits <- findOverlaps(posRanges, logRanges)
  logRX[queryHits(hits)] <- logr[subjectHits(hits),-(1:4)]
  logR <- cbind(pos,logRX)
  setnames(logR,colnames(logr))
  logR[,end:=NULL]
  
  write(paste0(c('','Chr','Position',samples),collapse = '\t'),file=file.path(ascatInDir,'merged.tumour.LogR.txt'),sep = '\t')
  write.table(logR,file=file.path(ascatInDir,'merged.tumour.LogR.txt'),append=TRUE,quote=FALSE,sep='\t',row.names = FALSE,col.names = FALSE)
  
  logRNormal <- logR[,c(1:3,grep(normal_sample,colnames(logR))),with=FALSE]
  write(paste0('\tChr\tPosition\t',sampleID),file=paste0(ascatInDir,'/',normal_sample,'.LogR.txt'),sep = '\t')
  write.table(logRNormal,file=paste0(ascatInDir,'/',normal_sample,'.LogR.txt'),append=TRUE,quote=FALSE,sep='\t',row.names = FALSE,col.names = FALSE)
  
}



################################################################################
#### Continue with multi sample segmentation ####