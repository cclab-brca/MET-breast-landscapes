# Further filtering on neoantigen calls - collapse highly similar epitopes into one call

pvacSeqDir<-"/analysis/pvacseq/new-mutations/output/summary/03-filt-exp-vaf/"
pvacSeqFiles <- list.files(pvacSeqDir, pattern = "*.txt")

outputDir<-"/analysis/pvacseq/new-mutations/output/summary/04-GLD-OneHLA/"

for (f in pvacSeqFiles){
  cat("Processing: ", f, "\n")
  sample<-sapply(strsplit(f,"\\."),"[",1)
  
  pvacseqTable<-read.table(paste0(pvacSeqDir,f), header = T, sep = "\t", stringsAsFactors = F)
  pvacseqTable.orig<-pvacseqTable
  nrow(pvacseqTable)
  if (nrow(pvacseqTable)>0){
    
    #remove duplicates
    pvacseqTable$id=paste(pvacseqTable$Chromosome, pvacseqTable$Start,pvacseqTable$Stop, pvacseqTable$Transcript, pvacseqTable$HLA.Allele, pvacseqTable$MT.Epitope.Seq, sep="_")
    e=table(pvacseqTable$id)
    e[e>1]
    for (u in unique(pvacseqTable$id)){
      subsetDf <- pvacseqTable[pvacseqTable$id==u,]
      if (nchar(subsetDf$Reference)==1 && nchar(subsetDf$Variant)==1){
        
        pvacseqTable <- pvacseqTable[ !(pvacseqTable$id==u & pvacseqTable$Best.MT.Score!= min(subsetDf$Best.MT.Score)),]
        subsetDf <- pvacseqTable[pvacseqTable$id==u,]
        
        if (nrow(subsetDf)>1){
          pvacseqTable <- pvacseqTable[ !rownames(pvacseqTable) %in% rownames(subsetDf)[2:nrow(subsetDf)] , ]
        }
      }
    }
    nrow(pvacseqTable)
    
    pvacseqTable$id=paste(pvacseqTable$Chromosome, pvacseqTable$Start,pvacseqTable$Stop, pvacseqTable$MT.Epitope.Seq, pvacseqTable$HLA.Allele, sep="_")
    
    for (u in unique(pvacseqTable$id)){
      subsetDf <- pvacseqTable[pvacseqTable$id==u,]
      if (nrow(subsetDf)>1){
        pvacseqTable <- pvacseqTable[ !(pvacseqTable$id==u & pvacseqTable$Best.MT.Score!= min(subsetDf$Best.MT.Score)),]
        subsetDf <- pvacseqTable[pvacseqTable$id==u,]
        if (nrow(subsetDf)>1){
          pvacseqTable <- pvacseqTable[ !rownames(pvacseqTable) %in% rownames(subsetDf)[2:nrow(subsetDf)] , ]
        }
      }
    }
    nrow(pvacseqTable)
    
    pvacseqTable$id=paste(pvacseqTable$Chromosome, pvacseqTable$Start,pvacseqTable$Stop, pvacseqTable$MT.Epitope.Seq, sep="_")
    
    for (u in unique(pvacseqTable$id)){
      subsetDf <- pvacseqTable[pvacseqTable$id==u,]
      if (nrow(subsetDf)>1){
        pvacseqTable <- pvacseqTable[ !(pvacseqTable$id==u & pvacseqTable$Best.MT.Score!= min(subsetDf$Best.MT.Score)),]
        subsetDf <- pvacseqTable[pvacseqTable$id==u,]
        if (nrow(subsetDf)>1){
          pvacseqTable <- pvacseqTable[ !rownames(pvacseqTable) %in% rownames(subsetDf)[2:nrow(subsetDf)] , ]
        }
      }
    }
    nrow(pvacseqTable)
    
    pvacseqTable$status="KEEP"
    #remove highly similar sequences
    for (r in 1:nrow(pvacseqTable)){
      
      similar=pvacseqTable[agrep(pvacseqTable[r,"MT.Epitope.Seq"],pvacseqTable$MT.Epitope.Seq, max.distance = 0.3),]
      if (nrow(similar)>1) { 
        
        sim=pvacseqTable[agrep(pvacseqTable[r,"MT.Epitope.Seq"],pvacseqTable$MT.Epitope.Seq),"MT.Epitope.Seq"]
        bestBinderScore <- min(similar$Best.MT.Score)
        toRemove <- rownames(similar[ similar$Best.MT.Score!= bestBinderScore,])
        pvacseqTable[rownames(pvacseqTable) %in% toRemove,"status"]="REMOVE"
        similar<- similar[!rownames(similar) %in% toRemove & similar$status=="KEEP",]
        
        #there might be cases where score is identical. In this case pick highest Predicted.Stability
        if (nrow(similar)>1){
          mostStableTranscript <- max(similar$Predicted.Stability)
          toRemove <- rownames(similar[ similar$Predicted.Stability!= mostStableTranscript,])
          pvacseqTable[rownames(pvacseqTable) %in% toRemove,"status"]="REMOVE"
        }
      }
    }
    pvacseqTable=pvacseqTable[pvacseqTable$status=="KEEP",]
    pvacseqTable=pvacseqTable[,-ncol(pvacseqTable)]
    pvacseqTable=pvacseqTable[,-ncol(pvacseqTable)]
    
    nrow(pvacseqTable)
  }
  write.table(pvacseqTable, file=paste0(outputDir,sample,".txt"), row.names=F, quote = F, sep = "\t")
  
}
