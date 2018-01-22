# Filters neoantigen calls using criteria suggested by pvacseq

library(IRanges)
library(VariantAnnotation)

#Prior to looping, load gene and transcript FPKM files
geneFpkmFile <- "/analysis/datasets/rna/counts/all-rpkm.txt"
geneFpkmMasterData <- read.table(geneFpkmFile,header = T, row.names = 1, sep = "\t", stringsAsFactors = F)

transcriptFpkmFile <- "/analysis/datasets/rna/cufflinks-transcript/cufflinks_isoforms.fpkm_tracking_counts_all.RData"
load(transcriptFpkmFile)
transcriptFpkmMasterData<-cufflinks_isoforms.fpkm_tracking_counts_all
rownames(transcriptFpkmMasterData)<-transcriptFpkmMasterData$TranscriptID
transcriptFpkmMasterData<-transcriptFpkmMasterData[,-1]

rnaVcfDir <- "/analysis/datasets/rna/rna/mutations/dna-muts-in-rna/"
dnaVcfDir <- "analysis/datasets/mutations/exome/vcf-filtered/"

#Identify neoantigen files
pvacSeqDir <- "/analysis/datasets/rna/neoantigens/pvacseq-netmhcpan-only/output/unfiltered/original-output/"
pvacSeqFiles <- list.files(pvacSeqDir, pattern = "*.final.tsv")

output.annotated<-"/analysis/datasets/rna/neoantigens/pvacseq-netmhcpan-only/output/unfiltered/annotated-output/"
output.filtered<-"/analysis/datasets/rna/neoantigens/pvacseq-netmhcpan-only/output/filtered/01-expression-vaf/"

.matrixOfListsToArray <- function(x) {
  n <- elementNROWS(x)
  maxn <- max(n)
  idx <- n < maxn
  x[idx] <- lapply(x[idx], function(a){c(a, rep(NA, maxn-length(a)))})
  x <- array(unlist(x), dim=c(maxn, nrow(x), ncol(x)),
             dimnames=list(NULL, rownames(x), colnames(x)))
  x <- aperm(x, c(2,3,1))t
  x
}


for (f in pvacSeqFiles){
  cat("Processing: ", f, "\n")
  sample<-sapply(strsplit(f,"\\."),"[",1)
  pvacseqTable<-read.table(paste0(pvacSeqDir,f), header = T, sep = "\t", stringsAsFactors = F)
  pvacseqTable[,sapply(pvacseqTable,class) == "logical"] <-
    sapply(pvacseqTable[,sapply(pvacseqTable,class) == "logical"],
           function(i) substr(as.character(i),1,1))

  #have any neoantigens been predicted?
  if ( nrow(pvacseqTable) > 0 ) {
    
    # has this sample been rna sequenced?
    if ( sample %in% colnames(geneFpkmMasterData) ) {
      # Merge Gene FPKM values
      geneFpkmData <- geneFpkmMasterData[, grepl(sample,colnames(geneFpkmMasterData)),drop=F]
      geneFpkmData= geneFpkmData[  rownames(geneFpkmData) %in% pvacseqTable$Ensembl.Gene.ID,,drop=F]
      geneFpkmData$Ensembl.Gene.ID<-rownames(geneFpkmData)
      pvacseqTable <- merge(x=pvacseqTable, y=geneFpkmData, by="Ensembl.Gene.ID")
      pvacseqTable$Gene.Expression<- pvacseqTable[,ncol(pvacseqTable)]
      pvacseqTable<- pvacseqTable[,-ncol(pvacseqTable)]
      
      t1<-pvacseqTable[,c(1,28)]
      t1<-t1[order(t1$Ensembl.Gene.ID),]
      t1<-t1[!duplicated(t1),]
      t1$master<-geneFpkmData[,1]
      stopifnot(sum(t1$Gene.Expression-t1$master)==0)
      
      # Merge Transcript FPKM values
      transcriptFpkmData <- transcriptFpkmMasterData[, grepl(sample,colnames(transcriptFpkmMasterData)),drop=F]
      transcriptFpkmData= transcriptFpkmData[  rownames(transcriptFpkmData) %in% pvacseqTable$Transcript,,drop=F]
      transcriptFpkmData$Transcript<-rownames(transcriptFpkmData)
      # this will remove rows without a known transcript
      pvacseqTable <- merge(x=pvacseqTable, y=transcriptFpkmData, by="Transcript")
      
      pvacseqTable$Transcript.Expression<- pvacseqTable[,ncol(pvacseqTable)]
      pvacseqTable<- pvacseqTable[,-ncol(pvacseqTable)]
      pvacseqTable$Transcript.Expression<-as.numeric(pvacseqTable$Transcript.Expression)
      
      t1<-pvacseqTable[,c(1,29)]
      t1<-t1[order(t1$Transcript),]
      t1<-t1[!duplicated(t1),]
      t1$master<-as.numeric(transcriptFpkmData[,1])
      stopifnot(sum(t1$Transcript.Expression-t1$master)==0)
      
      
      # load combined DNA and RNA VCF file and filter
      patient<-sapply(strsplit(sample,"_"),"[",1)
      param <- ScanVcfParam()
      vcf <- readVcf(paste0(rnaVcfDir,patient,".dna-rna.vcf"), "hg19", param)
      
      pvacseqTable$len<-nchar(pvacseqTable$Reference)-nchar(pvacseqTable$Variant)
      pvacseqTable[pvacseqTable$len==0,"Start"]<-pvacseqTable[pvacseqTable$len==0,"Stop"]
      pvacseqTable<-pvacseqTable[,-c(ncol(pvacseqTable))]
      pvacseqTable$variantID<-paste0(pvacseqTable$Chromosome,":",pvacseqTable$Start,"_",pvacseqTable$Reference,"/",pvacseqTable$Variant)
      vcf <- vcf[rownames(vcf)%in% unique(pvacseqTable$variantID)]
      
      #check that all the same
      vcfDifference=setdiff(unique(pvacseqTable$variantID),rownames(vcf))
      stopifnot(vcfDifference==0)
      
      vcf<- vcf[,colnames(vcf) %in% grep(paste0("BC|",sample),colnames(vcf),value=T)]
      depth<-geno(vcf)$DP
      altCounts <- .matrixOfListsToArray(geno(vcf)$AD)[,,2]
      refCounts <- .matrixOfListsToArray(geno(vcf)$AD)[,,1]
      vaf<-altCounts/depth
      colnames(vaf)<-paste0("vaf.",colnames(vaf))
      colnames(depth)<-paste0("depth.",colnames(depth))
      
      mutMatrix <- data.frame(cbind(depth,vaf))
      
      colnames(mutMatrix)[grepl("(?=.*BC)(?=.*DNA)(?=.*vaf)", colnames(mutMatrix), perl = TRUE)] <- "Normal.VAF.1"
      colnames(mutMatrix)[grepl("(?=.*BC)(?=.*DNA)(?=.*depth)", colnames(mutMatrix), perl = TRUE)] <- "Normal.Depth.1"
      colnames(mutMatrix)[grepl("(?=.*PT)(?=.*DNA)(?=.*depth)", colnames(mutMatrix), perl = TRUE)] <- "Tumor.DNA.Depth.1"
      colnames(mutMatrix)[grepl("(?=.*PT)(?=.*RNA)(?=.*depth)", colnames(mutMatrix), perl = TRUE)] <- "Tumor.RNA.Depth.1"
      colnames(mutMatrix)[grepl("(?=.*PT)(?=.*DNA)(?=.*vaf)", colnames(mutMatrix), perl = TRUE)] <- "Tumor.DNA.VAF.1"
      colnames(mutMatrix)[grepl("(?=.*PT)(?=.*RNA)(?=.*vaf)", colnames(mutMatrix), perl = TRUE)] <- "Tumor.RNA.VAF.1"
      
      mutMatrix <- mutMatrix[,c("Tumor.DNA.Depth.1","Tumor.DNA.VAF.1", "Tumor.RNA.Depth.1","Tumor.RNA.VAF.1", "Normal.Depth.1","Normal.VAF.1")]
      
      mutMatrix$variantID<-rownames(mutMatrix)
      
      pvacseqTable <- merge(x=pvacseqTable, y=mutMatrix, by="variantID")
      pvacseqTable$Tumor.DNA.Depth<-pvacseqTable$Tumor.DNA.Depth.1
      pvacseqTable$Tumor.DNA.VAF<-pvacseqTable$Tumor.DNA.VAF.1
      pvacseqTable$Tumor.RNA.Depth<-pvacseqTable$Tumor.RNA.Depth.1
      pvacseqTable$Tumor.RNA.VAF<-pvacseqTable$Tumor.RNA.VAF.1
      
      pvacseqTable$Normal.Depth<-pvacseqTable$Normal.Depth.1
      pvacseqTable$Normal.VAF<-pvacseqTable$Normal.VAF.1
      
      pvacseqTable<-pvacseqTable[,-c( (ncol(pvacseqTable)-6) : ncol(pvacseqTable))]
      pvacseqTable <- pvacseqTable[,c(4:8,2,3,9:ncol(pvacseqTable))]
      
      pvacseqTable<- pvacseqTable[order(pvacseqTable$Chromosome),]
      
      # Save this file as an unfiltered file, prior to more processing
      write.table(pvacseqTable,paste0(output.annotated,sample,".annotated.tsv"), sep = "\t", quote = F,col.names = T,row.names = F)
      
      
      # Commence filtering
      
      # A) keep expressed genes
      pvacseqTable=pvacseqTable[pvacseqTable$Gene.Expression>=1,]
      # B) keep expressed transcripts
      pvacseqTable=pvacseqTable[pvacseqTable$Transcript.Expression>=1,]
      # C) Normal coverage >= 5 as per pvacseq paper
      pvacseqTable <- pvacseqTable[pvacseqTable$Normal.Depth>5,]
      # D) Tumor DNA VAF Cutoff > 0 
      pvacseqTable <- pvacseqTable[pvacseqTable$Tumor.DNA.VAF>0,]
      # E) Tumor RNA VAF Cutoff > 0 
      pvacseqTable<-pvacseqTable[!is.na(pvacseqTable$Tumor.RNA.VAF),]
      pvacseqTable <- pvacseqTable[pvacseqTable$Tumor.RNA.VAF>0,]
      
      write.table(pvacseqTable,paste0(output.filtered,sample,".filtered-exp-vaf.txt"), sep = "\t", quote = F,col.names = T,row.names = F)
    }
    else {
      # there is a neoantigen file, but no RNA, so can't filter on RNA.
      # in this case, we can do rudimental RNA filtering based on the expression of other cases
      patient<-sapply(strsplit(sample,"_"),"[",1)
      
      #include DNA data
      param <- ScanVcfParam()
      vcf <- readVcf(paste0(dnaVcfDir,patient,".combined.vcf"), "hg19", param)
      
      pvacseqTable$len<-nchar(pvacseqTable$Reference)-nchar(pvacseqTable$Variant)
      pvacseqTable[pvacseqTable$len==0,"Start"]<-pvacseqTable[pvacseqTable$len==0,"Stop"]
      pvacseqTable<-pvacseqTable[,-c(ncol(pvacseqTable))]
      pvacseqTable$variantID<-paste0(pvacseqTable$Chromosome,":",pvacseqTable$Start,"_",pvacseqTable$Reference,"/",pvacseqTable$Variant)
      vcf <- vcf[rownames(vcf)%in% unique(pvacseqTable$variantID)]
      
      #check that all the same
      vcfDifference=setdiff(unique(pvacseqTable$variantID),rownames(vcf))
      stopifnot(vcfDifference==0)
      
      vcf<- vcf[,colnames(vcf) %in% grep(paste0("BC|",sample),colnames(vcf),value=T)]
      depth<-geno(vcf)$DP
      altCounts <- .matrixOfListsToArray(geno(vcf)$AD)[,,2]
      refCounts <- .matrixOfListsToArray(geno(vcf)$AD)[,,1]
      vaf<-altCounts/depth
      
      colnames(vaf)<-paste0("vaf.",colnames(vaf))
      colnames(depth)<-paste0("depth.",colnames(depth))
      
      mutMatrix <- data.frame(cbind(depth,vaf))
      
      colnames(mutMatrix)[grepl("(?=.*BC)(?=.*vaf)", colnames(mutMatrix), perl = TRUE)] <- "Normal.VAF.1"
      colnames(mutMatrix)[grepl("(?=.*BC)(?=.*depth)", colnames(mutMatrix), perl = TRUE)] <- "Normal.Depth.1"
      colnames(mutMatrix)[grepl("(?=.*PT)(?=.*depth)", colnames(mutMatrix), perl = TRUE)] <- "Tumor.DNA.Depth.1"
      colnames(mutMatrix)[grepl("(?=.*PT)(?=.*vaf)", colnames(mutMatrix), perl = TRUE)] <- "Tumor.DNA.VAF.1"
      mutMatrix$Tumor.RNA.Depth.1 <- ""
      mutMatrix$Tumor.RNA.VAF.1 <-""
      mutMatrix$variantID<-rownames(mutMatrix)
      
      pvacseqTable <- merge(x=pvacseqTable, y=mutMatrix, by="variantID")
      pvacseqTable$Tumor.DNA.Depth<-pvacseqTable$Tumor.DNA.Depth.1
      pvacseqTable$Tumor.DNA.VAF<-pvacseqTable$Tumor.DNA.VAF.1
      pvacseqTable$Tumor.RNA.Depth<-pvacseqTable$Tumor.RNA.Depth.1
      pvacseqTable$Tumor.RNA.VAF<-pvacseqTable$Tumor.RNA.VAF.1
      
      pvacseqTable$Normal.Depth<-pvacseqTable$Normal.Depth.1
      pvacseqTable$Normal.VAF<-pvacseqTable$Normal.VAF.1
      
      pvacseqTable<-pvacseqTable[,-c( 1,(ncol(pvacseqTable)-6) : ncol(pvacseqTable))]
      pvacseqTable<- pvacseqTable[order(pvacseqTable$Chromosome),]
      
      write.table(pvacseqTable,paste0(output.annotated,sample,".no-rna.annotated.tsv"), sep = "\t", quote = F,col.names = T,row.names = F)
      
      
      #are there any related samples? 
      otherSamples=grep(patient,colnames(geneFpkmMasterData),value=T)
      if (length(otherSamples)>0){
        otherSamples=otherSamples[1]
        geneFpkmData <- geneFpkmMasterData[, grepl(otherSamples,colnames(geneFpkmMasterData)),drop=F]
        geneFpkmData= geneFpkmData[  rownames(geneFpkmData) %in% pvacseqTable$Ensembl.Gene.ID,,drop=F]
        geneFpkmData$Ensembl.Gene.ID<-rownames(geneFpkmData)
        pvacseqTable <- merge(x=pvacseqTable, y=geneFpkmData, by="Ensembl.Gene.ID")
        pvacseqTable$Gene.Expression<- pvacseqTable[,ncol(pvacseqTable)]
        pvacseqTable<- pvacseqTable[,-ncol(pvacseqTable)]
        
        # Merge Transcript FPKM values
        transcriptFpkmData <- transcriptFpkmMasterData[, grepl(otherSamples,colnames(transcriptFpkmMasterData)),drop=F]
        transcriptFpkmData= transcriptFpkmData[  rownames(transcriptFpkmData) %in% pvacseqTable$Transcript,,drop=F]
        transcriptFpkmData$Transcript<-rownames(transcriptFpkmData)
        # this will remove rows without a known transcript
        pvacseqTable <- merge(x=pvacseqTable, y=transcriptFpkmData, by="Transcript")
        pvacseqTable$Transcript.Expression<- pvacseqTable[,ncol(pvacseqTable)]
        pvacseqTable<- pvacseqTable[,-ncol(pvacseqTable)]
        pvacseqTable$Transcript.Expression<-as.numeric(pvacseqTable$Transcript.Expression)
        nrow(pvacseqTable)
        # Commence filtering
        
        pvacseqTable=pvacseqTable[pvacseqTable$Gene.Expression>=1,]
        pvacseqTable=pvacseqTable[pvacseqTable$Transcript.Expression>=1,]
        pvacseqTable <- pvacseqTable[pvacseqTable$Normal.Depth>5,]
        pvacseqTable <- pvacseqTable[pvacseqTable$Tumor.DNA.VAF>0,]
        pvacseqTable$Gene.Expression=""
        pvacseqTable$Transcript.Expression=""
        nrow(pvacseqTable)
        pvacseqTable <- pvacseqTable[,c(3:7,1,2,8,9:ncol(pvacseqTable))]
        write.table(pvacseqTable,paste0(output.filtered,sample,".filtered-exp-vaf-rnasurrogate.txt"), sep = "\t", quote = F,col.names = T,row.names = F)
      } else {
        write.table(pvacseqTable,paste0(output.filtered,sample,".filtered-exp-vaf-NO-rnasurrogate.txt"), sep = "\t", quote = F,col.names = T,row.names = F)
      }
    }
  }
  else {
    pvacseqTable<-pvacseqTable[,c(1:ncol(pvacseqTable)-1)]
    write.table(pvacseqTable,paste0(output.annotated,sample,".no-NA.annotated.tsv"), sep = "\t", quote = F,col.names = T,row.names = F)
    write.table(pvacseqTable,paste0(output.filtered,sample,".filtered-exp-vaf.txt"), sep = "\t", quote = F,col.names = T,row.names = F)
  }
