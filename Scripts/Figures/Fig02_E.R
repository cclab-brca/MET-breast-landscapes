# Compare expression of stem vs clade vs private mutations

suppressMessages(library(VariantAnnotation))

inputVCFDir="/mutations/rna/1rna-muts-from-dna/"

muts.file<-"/mutations/exome/somatic/all-filtered.maf"
mutations <- read.table(muts.file, sep="\t",header=TRUE, fill=TRUE,comment.char="#",quote="",
                        blank.lines.skip=T, row.names = NULL, stringsAsFactors = FALSE)
mutations$patient<-sapply(strsplit(mutations$sample,"-"),"[",1)
mutations$ID=paste0(mutations$Chr,":",mutations$Start,"_",mutations$Ref_Allele,"/",mutations$Tumor_Alt_Allele)
mutations<-mutations[,c(30,4,12)]
mutations<-mutations[!duplicated(mutations),]

mutClassifierCols<-c("#FC4A1A","#F7B733","#4ABDAC")

geneExpressionFile="/htseq/exon-gene_id-union/tpm-hugo.txt"
gene_expression<-read.table(geneExpressionFile,row.names=1,header=TRUE, sep="\t", dec = ".",check.names=FALSE)
colnames(gene_expression)<-gsub("X","",colnames(gene_expression))
colnames(gene_expression)<-gsub("\\.","-",colnames(gene_expression))


vcfFiles <- list.files(inputVCFDir, pattern = "*.dna-rna.vcf")

allMutantTpmScores<-numeric()
sampleSet<-character()

.matrixOfListsToArray <- function(x) {
  n <- elementNROWS(x)
  maxn <- max(n)
  idx <- n < maxn
  x[idx] <- lapply(x[idx], function(a){c(a, rep(NA, maxn-length(a)))})
  x <- array(unlist(x), dim=c(maxn, nrow(x), ncol(x)),
             dimnames=list(NULL, rownames(x), colnames(x)))
  x <- aperm(x, c(2,3,1))
  x
}


library(reshape)
for (vcfFile in vcfFiles){
  cat (vcfFile,"\n")
  vcfPath <- paste0(inputVCFDir,vcfFile)
  patient <-  sapply(strsplit( vcfFile,"\\."),"[[",1)
  param <- ScanVcfParam()
  vcf <- readVcf(vcfPath, "hg19", param)
  vcf <- unique(vcf)
  normalSamples<-c("288-024","290-004","291-005","298-023","308-024","315-022","323-007","328-010","330-007","DET52-normal")
  vcf=vcf[, !grepl(paste(normalSamples,collapse = "|"),colnames(vcf)) ]
  samples=colnames(vcf)
  
  #extract depth and vaf for all samples
  depth<-geno(vcf)$DP
  altCounts <- .matrixOfListsToArray(geno(vcf)$AD)[,,2]
  vaf<-altCounts/depth
  
  #classify mutations into Stem/clade/private
  dnaVaf<-vaf[,grep("DNA",colnames(vaf)),drop=F]
  colnames(dnaVaf)<-gsub("_DNA","",colnames(dnaVaf))
  noSamples=ncol(dnaVaf)
  vafCeil=ceiling(dnaVaf)
  vafCeil[is.na(vafCeil)] <- 0
  vafCeil <- rowSums(vafCeil)
  
  if (noSamples>2){
    mutation.stem <- names(vafCeil[vafCeil==noSamples])
    mutation.clade <- names(vafCeil[vafCeil<noSamples & vafCeil>1])
    mutation.private <- names(vafCeil[vafCeil==1])
  } else if (noSamples ==2){
    mutation.stem <- names(vafCeil[vafCeil==noSamples])
    mutation.clade <- c()
    mutation.private <- names(vafCeil[vafCeil==1])
  } else if (noSamples ==1){
    mutation.stem <- c()
    mutation.clade <- c()
    mutation.private <- names(vafCeil[vafCeil==1])
  }
  
  mutationStatus<-data.frame(mutationID=c(mutation.stem,mutation.clade,mutation.private), 
                             status=c(rep("Stem",length(mutation.stem)), 
                                      rep("Clade",length(mutation.clade)),
                                      rep("Private",length(mutation.private))))
  
  
  
  #filter RNA vaf, keep only mutations that have at least one read in the RNA
  rnaDepth=depth[,grep("RNA",colnames(depth)),drop=F]
  colnames(rnaDepth)<-gsub("_RNA","",colnames(rnaDepth))
  head(rnaDepth)
  library(reshape2)
  rnaDepthMelt <- melt(rnaDepth)
  rnaDepthMelt<-rnaDepthMelt[complete.cases(rnaDepthMelt),]
  rnaDepthMelt<-rnaDepthMelt[rnaDepthMelt$value>1,]
  rnaDepthMelt$masterID<-paste0(rnaDepthMelt$X1,"_",rnaDepthMelt$X2)
  head(rnaDepthMelt)
  nrow(rnaDepthMelt)
  
  rnaVaf<-vaf[,grep("RNA",colnames(vaf)),drop=F]
  colnames(rnaVaf)<-gsub("_RNA","",colnames(rnaVaf))
  head(rnaVaf)
  
  exV1<-melt(rnaVaf)
  colnames(exV1)<-c("ID","sample","rnaVaf")
  exV1=exV1[complete.cases(exV1),]
  
  #retain only mutations that have RNA coverage
  exV1=exV1[paste0(exV1$ID,"_",exV1$sample) %in% as.character(rnaDepthMelt$masterID),]
  dnaVafMelt<- melt(dnaVaf)
  dnaVafMelt<-dnaVafMelt[dnaVafMelt$value>0,]
  exV1=exV1[paste0(exV1$ID,"_",exV1$sample) %in% paste0(dnaVafMelt$X1,"_",dnaVafMelt$X2),]
  
  rnaVaf=merge(exV1,mutationStatus, by.x="ID",by.y="mutationID")
  
  #merge metadata file
  rnaVaf=merge(rnaVaf,mutations,by.x="ID",by.y="ID")
  rnaVaf$sample<-as.character(rnaVaf$sample)
  sampleSet=c(sampleSet,unique(rnaVaf$sample))
  
  for (s in unique(rnaVaf$sample)){
    genes<- rnaVaf[rnaVaf$sample==s,]
    if (length(gene_expression[[s]])!=0){
      GE <- gene_expression[rownames(gene_expression) %in% genes$Hugo_Symbol ,s,drop=F]
      GE <- merge(x=GE,y=genes,by.x=0,by.y="Hugo_Symbol", sort=F)
      GE <- GE[,c(1,2,5,6,7)]
      colnames(GE)<-c("Hugo_Symbol","tpm","vaf","status","variant")
      GE$correctedTpm<-GE$tpm*GE$vaf
      mGE<-mean(GE$correctedTpm)
      sGE<-sd(GE$correctedTpm)
      GE$Z1<-(GE$correctedTpm-mGE)/sGE
      k<-1
      MIG<-NA
      for (gen in c("Stem","Clade","Private")) {
        MIG[k]<-mean(GE[GE$status==gen,"Z1"],na.rm=T)
        k<-k+1
      }
      cat(MIG,"\n")
      allMutantTpmScores=rbind(allMutantTpmScores,MIG)
    }
  }
}

rownames(allMutantTpmScores)<-sampleSet
colnames(allMutantTpmScores)<-c("Stem","Clade","Private")
head(allMutantTpmScores)

allMutantTpmScoresDf<-data.frame(allMutantTpmScores)
allMutantTpmScoresDf$sample<-sapply(strsplit(rownames(allMutantTpmScoresDf),"-"),"[",2)
allMutantTpmScoresDf$patient<-sapply(strsplit(rownames(allMutantTpmScoresDf),"-"),"[",1)
allMutantTpmScoresDf<-melt(allMutantTpmScoresDf, id.vars = c("patient","sample"))

allMutantTpmScoresDf=allMutantTpmScoresDf[allMutantTpmScoresDf$patient!="291",]

pdf ("~/fig2e.pdf", height=4, width=4)

mutClassifierCols<-c("darkblue","tomato","darkgreen")


print(
  ggplot(allMutantTpmScoresDf)+
    geom_boxplot(aes(x=variable,y=value, fill=variable))+
    labs(x="Mutation classification",y="Zscore normalised  mutant allele TPM")+
    scale_fill_manual(name="Mutation\nclassification",values=mutClassifierCols)+
    theme_bw(base_size = 11)+
    theme(legend.key = element_blank(),
          plot.title =element_text(face="bold"),
          axis.title.x = element_text(face="bold"),
          axis.title.y = element_text(face="bold", margin=margin(0,10,0,0)),
          legend.title = element_text(face="bold"),
          legend.position="bottom",
          panel.border = element_rect(fill=NA, colour = "black", 
                                      size=1))+
    guides(fill=F)
)

dev.off()