rawdatadir <- '~/Documents/Caldas/Data/dataRaw/targetedSeq/'
treeomicsDir_targeted <- '~/Documents/Caldas/Data/dataProcessed/treeomics/output_targeted_all/'
outdatadir <- '~/Documents/Caldas/Data/dataProcessed/treeomics/input_targeted_selected/'
sampledir <- '~/Documents/Caldas/Data/dataRaw'

dir.create(outdatadir,showWarnings = FALSE)

library("data.table")

CASES <- c(288,290,'290exOvary',291,298,'298main',308,315,323,328,330,"DET52") ## TODO add DET52
SAMPLESETS <- c('All','Biopsies','Fluids','WESsamples')

normalSamples <- read.table(paste0(sampledir,'/samples_germline.txt'),stringsAsFactors = FALSE)$V1
normalSamples <- data.frame(case=gsub('(.*)-.*','\\1',normalSamples),
                            normal=gsub('.*-(.*)','\\1',normalSamples))

for (case in CASES) {
    
  if (case=='290exOvary') {
    baseCase=290
  } else if (case=='298main') {
    baseCase=298
  } else {
    baseCase=case
  }
  
  ## read total depth data
  if (case=='DET52') {
    depthFile <- paste0(rawdatadir,'merge.depth/ncomms2015.0.panel.snv.sum.dep.txt')
  } else {
    depthFile <- paste0(rawdatadir,'merge.depth/',baseCase,'.sum.depth.txt')
  }
  depthDat <- read.table(depthFile,sep = '\t',header=TRUE,stringsAsFactors = FALSE)#,quote="")
  ## read alt depth data
  if (case=='DET52') {
    altFile <- paste0(rawdatadir,'merge.reads/ncomms2015.0.panel.snv.sum.reads.txt')
  } else {
    altFile <- paste0(rawdatadir,'merge.reads/',baseCase,'.sum.reads.txt')  
  }
  altDat <- read.table(altFile,sep = '\t',header=TRUE,stringsAsFactors = FALSE)#,quote="")
  
  if (case=='DET52') {
    colnames(depthDat) <- colnames(altDat) <- c('Mutation_Classification','Consequence','amplicon','ID',
                                                'DET52.Control_G1','DET52.Control_G2','DET52.Control_G3',
                                                'DET52.Control_G4','DET52.Control_G5','DET52.Control_P1',
                                                'DET52.Control_P2','DET52.Control_P3','DET52.Control_P4',
                                                'DET52.Control_P5','DET52.Control_P6','DET52.Control_P7',
                                                'DET52.Control_P8','DET52.DCIS','DET52.MT2','DET52.MT3',
                                                'DET52.MT4','DET52.MT6','DET52.MT5','DET52.MT7','DET52.MT8',
                                                'DET52.normal','DET52.T1','DET52.T2','DET52.T3','DET52.M6',
                                                'DET52.T4','DET52.T5','DET52.T6','DET52.T7','DET52.T8','DET52.T9')
    depthDat <- depthDat[,!colnames(depthDat)%in%c('Mutation_Classification','Consequence')]
    altDat <- altDat[,!colnames(altDat)%in%c('Mutation_Classification','Consequence')]
  }  else {
    colnames(depthDat)[1:2] <- colnames(altDat)[1:2] <- c('amplicon','ID')
  }
  
  if (case=='DET52') {
    metadata <- data.frame(Chromosome=gsub(pattern = "(chr[0-9X]*).*",replacement = "\\1",x = depthDat$ID),
                           Position=gsub(pattern = "chr[0-9X]*:([0-9]*):.*",replacement = "\\1",x = depthDat$ID),
                           Change=gsub(pattern = "chr[0-9X]*:[0-9]*:(.*):(.*)", replacement = "\\1>\\2", x = depthDat$ID),
                           Gene=gsub(pattern = "_D0057_00[12]", replacement = "", x = depthDat$amplicon))
  } else {
    metadata <- data.frame(Chromosome=gsub(pattern = "(chr[0-9X]*).*",replacement = "\\1",x = depthDat$ID),
                           Position=gsub(pattern = "chr[0-9X]*:([0-9]*):.*",replacement = "\\1",x = depthDat$ID),
                           Change=gsub(pattern = "chr[0-9X]*:[0-9]*:(.*):(.*)", replacement = "\\1>\\2", x = depthDat$ID),
                           Gene=gsub(pattern = "[0-9]*_([^_]*).*", replacement = "\\1", x = depthDat$amplicon))
  }
  
  ## only use sites that are in read depth filtered file.
  
  ## get index of these sites
  # raw data dir (used for filtering)
  filteredFile <- list.files(path=paste0(sampledir,'/targetedSeq/merge.calls.sum.500-exomeSelect-normalAsControl'),pattern=paste0(baseCase,'.500..*sum.variants.txt'),full.names = TRUE) 
  filtered <- fread(filteredFile)
  ID <- paste0(gsub("chr","",metadata$Chromosome),":",metadata$Position)
  filterIdx <- ID%in%filtered$ID
  
  ## remove all other sites from data
  metadata <- metadata[filterIdx,]
  depthDat <- depthDat[filterIdx,]
  altDat <- altDat[filterIdx,]

  ## remove anything but count data
  excludeIdx <- c(1,2,grep("Control",colnames(depthDat)))
  covCount <- depthDat[,-excludeIdx]
  altCount <- altDat[,-excludeIdx]
  colnames(covCount) <- colnames(altCount) <- gsub(paste0("X?",baseCase,"[_.]","(.*)"),"X\\1",colnames(covCount)) ## \\.
  colnames(covCount) <- colnames(altCount) <- gsub("\\.","_",colnames(covCount))
  
  altCount <- data.frame(metadata,altCount)
  covCount <- data.frame(metadata,covCount)
 
  for (sampleset in SAMPLESETS) {
  
    ## read posterior probabilities from previous treeomics run (used for filtering)
    postFile <- list.files(path=paste0(treeomicsDir_targeted,'/X',case,'_',sampleset),pattern='posterior.txt',full.names = TRUE)
    posterior <- read.csv(postFile,sep='\t',stringsAsFactors = FALSE)
    columnNames <- colnames(posterior)
    posterior <- read.csv(postFile,sep='\t',stringsAsFactors = FALSE,header = TRUE,skip = 4)
    ## remove rows (purities, VAFs and priors)
    colnames(posterior) <- c('Chromosome','Position','Position2','Gene.Name',columnNames[-(1:4)])
    ## filter posterior
    posteriorFilterIdx <- paste0(posterior$Chromosome,":",posterior$Position)%in%filtered$ID
    posterior <- posterior[posteriorFilterIdx,]

    # select samples
    altCountCurr <- as.data.table(altCount[,c('Chromosome','Position','Change','Gene',colnames(posterior)[-(1:4)],paste0('X',normalSamples$normal[normalSamples$case==baseCase]))])
    covCountCurr <- as.data.table(covCount[,c('Chromosome','Position','Change','Gene',colnames(posterior)[-(1:4)],paste0('X',normalSamples$normal[normalSamples$case==baseCase]))])
    
    ## reduce data set size for some of the cases
    exclude <- NA
    if (case==288 & sampleset=='All') {
      exclude <- c('X021','Xascitic','X13A_0267A28','X13A_0267A45','XPL_aut','XPL_OnTTM','XCSF')
    }
    if (baseCase==290 & sampleset=='All') {
      exclude <- c('XPL_A','XPL_B','XPL_C','XCSF','XPR1_Primary','Xascitic')
    }
    if (case==291 & sampleset=='All') {
      exclude <- c('XPL_aut','Xperic','Xpleural','X019','X13A_0287A45_Lung_4','X13A_0287A8_Ovary_Left')
    }
    if (case==308 & sampleset%in%c('WESsamples')) {
      exclude <- c('X006','X015','X017','X018','X019','X020','XPR2_Primary','Xpleural','XCSF')
    }
    if (case==308 & sampleset%in%c('Biopsies')) {
      exclude <- c('X006','X015','X016','X017','X018','X019','X020','XPR2_Primary','Xpleural','XCSF')
    }
    if (case==308 & sampleset%in%c('All')) {
      exclude <- c('X006','X015','X016','X017','X018','X019','X020','X022','XPR2_Primary','Xpleural','XCSF','XPL_aut')
    }
    if (case==330 & sampleset=='All') {
      exclude <- c('XPL_aut','X003')
    }
    if (case=='DET52' & sampleset=='All') {
      exclude <- c('XM6','XT3','XT5','XT8')
    }
    if (!all(is.na(exclude))) {
      altCountCurr <- altCountCurr[,!colnames(altCountCurr)%in%exclude,with=FALSE]
      covCountCurr <- covCountCurr[,!colnames(covCountCurr)%in%exclude,with=FALSE]
      posterior <- posterior[,!colnames(posterior)%in%exclude]
    }
    
    ## match rows between posterior, covCount and altCount
    posterior <- as.data.table(posterior)
    
    ## FILTER MUTATIONS
    ## number of mutated samples
    Nmut <- apply(posterior[,-(1:4),with=FALSE],1,function(x) sum(as.numeric(x)>0.5))
    # Nmut <- apply(posterior[,-(1:4),with=FALSE],1,function(x) sum(as.numeric(x[-(1:4)])>0.5))
    ## remove sites where no sample has a probability of being mutated above 0.5
    posterior <- posterior[Nmut>=1,]
    # altCountCurr <- altCountCurr[Nmut>=1,]
    # covCountCurr <- covCountCurr[Nmut>=1,]
    
    setkey(posterior,Chromosome,Position)
    posterior$ID <- paste0('chr',posterior$Chromosome,'_', posterior$Position)
    posterior$ID2 <- paste0('chr',posterior$Chromosome,'_', posterior$Position2)
    
    setkey(altCountCurr,Chromosome,Position)
    altCountCurr$ID <- paste0(altCountCurr$Chromosome,'_', altCountCurr$Position)
    altCountCurr <- altCountCurr[altCountCurr$ID%in%posterior$ID | altCountCurr$ID%in%posterior$ID2,]
    altCountCurr[,ID := NULL]
    
    setkey(covCountCurr,Chromosome,Position)
    covCountCurr$ID <- paste0(covCountCurr$Chromosome,'_', covCountCurr$Position)
    covCountCurr <- covCountCurr[covCountCurr$ID%in%posterior$ID | covCountCurr$ID%in%posterior$ID2,]
    covCountCurr[,ID := NULL]
    
    posterior[,c("ID","ID2") := NULL]
    

    
    ## write to file
    write.table(covCountCurr,paste0(outdatadir,'/X',case,'_',sampleset,'_covCount.txt'),sep='\t',row.names=FALSE,quote=FALSE,na = "n/a")
    write.table(altCountCurr,paste0(outdatadir,'/X',case,'_',sampleset,'_altCount.txt'),sep='\t',row.names=FALSE,quote=FALSE,na = "n/a")
  
  }
}