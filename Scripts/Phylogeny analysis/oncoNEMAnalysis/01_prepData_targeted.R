rawdatadir <- '~/Documents/Caldas/Data/dataRaw/targetedSeq/'
outdatadir <- '~/Documents/Caldas/Data/dataProcessed/treeomics/input_targeted_all/'
annotDatadir <- '~/Documents/Caldas/Data/dataAnnotation/'
options(stringsAsFactors = FALSE)

library(XLConnect)
cases <- c(288,290,'290exOvary',291,298,'298main',308,315,323,328,330,'DET52') ## 298Brain only has 1 sample for targeted deep seq

annotation <- readWorksheetFromFile(paste0(annotDatadir,"/Metadata_for_Edith_V4_annotated_mod2.xls"), 
                                    sheet=1, 
                                    startRow = 1,
                                    endCol = 9)
annotation$Sample.Name <- paste0("X",gsub("-",".",annotation$Sample.Name))


for (case in cases) {
  
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
  metadata <- data.frame(Chromosome=gsub(pattern = "(chr[0-9X]*).*",replacement = "\\1",x = depthDat$ID),
                         Position=gsub(pattern = "chr[0-9X]*:([0-9]*):.*",replacement = "\\1",x = depthDat$ID),
                         Change=gsub(pattern = "chr[0-9X]*:[0-9]*:(.*):(.*)", replacement = "\\1>\\2", x = depthDat$ID),
                         Gene=gsub(pattern = "[0-9]*_([^_]*).*", replacement = "\\1", x = depthDat$amplicon))
  
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
  
  ## use anything but control when selecting samples
  ##
  excludeIdx <- c(1,2,grep("Control",colnames(depthDat)))
  if (case=='290exOvary') {
    excludeIdx <- c(excludeIdx,grep("022",colnames(depthDat)))
  } else if (case==291) {
    excludeIdx <- c(excludeIdx,
                   grep("291_13A.0287A36_Liver",colnames(depthDat)))
  } else if (case=='298main') {
    excludeIdx <- c(excludeIdx,
                    grep("009",colnames(depthDat)),
                    grep("010",colnames(depthDat)))
  } else if (case==323) {
    excludeIdx <- c(excludeIdx,
                    grep("323_14A076.A8_Systemic",colnames(depthDat)))
  } else if (case==330) {
    excludeIdx <- c(excludeIdx,
                    grep("330_14A105.A1_Lung",colnames(depthDat)))
  }
  
  depth <- depthDat[,-excludeIdx]
  alt <- altDat[,-excludeIdx]
  if (case=='DET52') {
    colnames(depth) <- colnames(alt) <- gsub("DET52\\.","XDET52.",colnames(depth))
  } else {
    colnames(depth) <- colnames(alt) <- gsub(paste0("X.",baseCase,"[_]","(.*)\\."),paste0("X",baseCase,"_\\1"),colnames(depth)) ## column names of data set changed
    colnames(depth) <- colnames(alt) <- gsub(paste0("X.",baseCase,"[.]","(.*)\\."),paste0("X",baseCase,".\\1"),colnames(depth))
  }
  
  annotationIdx <- sapply(colnames(depth),function(x) grep(strtrim(x,19),x=annotation$Sample.Name)) 
  if (any(duplicated(annotationIdx))) {
    stop("Duplication!")
  }
  colnames(depth) <- colnames(alt) <- paste0("X",gsub(paste0("X",baseCase,"."),"",colnames(depth)))
  
  ## filter out sites that have missing values in normal sample (otherwise treeomics fails)
  normalIdx <- which(annotation$Tissue.Origin[annotationIdx] == 'Blood')
  keepIdx <- which(!is.na(alt[,normalIdx])&!is.na(depth[,normalIdx]))
  depth <- depth[keepIdx,]
  alt <- alt[keepIdx,]
  metadata <- metadata[keepIdx,]
  
  ## join with metadata
  depthAll <- data.frame(metadata,depth)
  altAll <- data.frame(metadata,alt)
  colnames(depthAll) <- colnames(altAll) <- gsub('\\.','_',colnames(depthAll))
  
  write.table(depthAll,paste0(outdatadir,'/X',case,'_covCount_All.txt'),sep='\t',row.names=FALSE,quote=FALSE,na = "n/a")
  write.table(altAll,paste0(outdatadir,'/X',case,'_altCount_All.txt'),sep='\t',row.names=FALSE,quote=FALSE,na = "n/a")
  
  ## select fluids
  fluidIdx <- which(annotation$Tissue.Origin[annotationIdx] %in% c('Fluid','Blood'))
  depthFluid <- data.frame(metadata,depth[,fluidIdx])
  altFluid <- data.frame(metadata,alt[,fluidIdx])
  colnames(depthFluid) <- colnames(altFluid) <- gsub('\\.','_',colnames(depthFluid))
  
  write.table(depthFluid,paste0(outdatadir,'/X',case,'_covCount_Fluids.txt'),sep='\t',row.names=FALSE,quote=FALSE,na = "n/a")
  write.table(altFluid,paste0(outdatadir,'/X',case,'_altCount_Fluids.txt'),sep='\t',row.names=FALSE,quote=FALSE,na = "n/a")
  
  ## select biopsies
  biopsyIdx <- which(!annotation$Tissue.Origin[annotationIdx] %in% c('Fluid'))
  depthBiopsy <- data.frame(metadata,depth[,biopsyIdx])
  altBiopsy <- data.frame(metadata,alt[,biopsyIdx])
  colnames(depthBiopsy) <- colnames(altBiopsy) <- gsub('\\.','_',colnames(depthBiopsy))
  
  write.table(depthBiopsy,paste0(outdatadir,'/X',case,'_covCount_Biopsies.txt'),sep='\t',row.names=FALSE,quote=FALSE,na = "n/a")
  write.table(altBiopsy,paste0(outdatadir,'/X',case,'_altCount_Biopsies.txt'),sep='\t',row.names=FALSE,quote=FALSE,na = "n/a")
  
  ## select WES samples
  wesIdx <- which(annotation$WES.available[annotationIdx]=='Yes'|annotation$Tissue.Origin[annotationIdx]=='Blood')
  depthWESsamples <- data.frame(metadata,depth[,wesIdx])
  altWESsamples <- data.frame(metadata,alt[,wesIdx])
  colnames(depthWESsamples) <- colnames(altWESsamples) <- gsub('\\.','_',colnames(depthWESsamples))
  
  write.table(depthWESsamples,paste0(outdatadir,'/X',case,'_covCount_WESsamples.txt'),sep='\t',row.names=FALSE,quote=FALSE,na = "n/a")
  write.table(altWESsamples,paste0(outdatadir,'/X',case,'_altCount_WESsamples.txt'),sep='\t',row.names=FALSE,quote=FALSE,na = "n/a")
  
}