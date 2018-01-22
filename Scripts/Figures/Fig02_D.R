#  Using the TS data plot percentage of stem/clade mutations is detected in the primary

library(ggplot2)
suppressMessages(library(VariantAnnotation))

mutClassifierCols<-c("darkblue","tomato","darkgreen")

inputVCFDir="/mutations/exome/somatic/03b-lookback-vcf-normfilter/"
outDir <- "/plots/"
vcfFiles <- list.files(inputVCFDir, pattern = "*.vcf")


vcfFiles=grep("290|291|308|323|328|DET52",vcfFiles,value = T)
vcfFiles=vcfFiles[!grepl("290-004ov.combined.vcf|298-023en.combined.vcf",vcfFiles)]

primarySamples<- c("290-PR1_Primary","290-PR2_Primary",
                   "291_13-106703-4",
                   "308-PR1_Primary","308-PR2_Primary",
                   "323-PR_Primary",
                   "328-PR_Prima",
                   "P-1","P-3")
statMx<-matrix(ncol=2, nrow=length(primarySamples),dimnames = list(primarySamples,c("Stem","Clade")))

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



allStatus<-c()
for (vcfFile in vcfFiles){
  
  cat (vcfFile,"\n")
  vcfPath <- paste0(inputVCFDir,vcfFile)
  patient <-  sapply(strsplit( vcfFile,"\\."),"[[",1)
  
  param <- ScanVcfParam()
  vcf <- readVcf(vcfPath, "hg19", param)
  
  vcf <- unique(vcf)
  
  rd=c("290-022","288-021","288-022","288-023",
              "291-019",
              "298-004","298-005","298-006","298-007","298-012","298-020","298-021","298-022",
              "330-007",
              "DET52-mt3","DET52-mt4","DET52-mt8","DET52-t1","DET52-t2","DET52-t9","288-024","290-004","291-005","298-023","308-024","315-022","323-007","328-010","330-007",
                    "DET52-normal")
  vcf=vcf[, !grepl(paste(rd,collapse = "|"),colnames(vcf)) ]

  samples=colnames(vcf)
  
  depth<-geno(vcf)$DP
  altCounts <- .matrixOfListsToArray(geno(vcf)$AD)[,,2]
  depth<-geno(vcf)$DP
  vaf<-altCounts/depth
  
  #classify mutations into Stem/clade/private
  noSamples=ncol(depth)
  vafCeil=ceiling(vaf)
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
  
  mutationStatus$case<-sapply(strsplit(patient,"-"),"[",1)
  mutationStatus$mutationID=sapply(strsplit(as.character(mutationStatus$mutationID),"_"),"[",1)
  
  raindanceDir<- "/targetted-seq/Multi-regional_20170627/seq_results/merge.freq/merge.calls.sum.500-exomeSelect-normalAsControl/"
  raindanceFile<-  Sys.glob(file.path(paste0(raindanceDir,sapply(strsplit(patient,"-"),"[",1), "*zscore.txt")))
  raindanceMuts<-read.table(raindanceFile,header=T, sep="\t",stringsAsFactors = F, quote ="\"")
  colnames(raindanceMuts)[1]<-"ID"
  colnames(raindanceMuts)<-gsub("X","",colnames(raindanceMuts))
  colnames(raindanceMuts)<-gsub("\\.","-",colnames(raindanceMuts))
  
  primaryMuts<-raindanceMuts[,colnames(raindanceMuts) %in% c("ID",primarySamples),drop=F]

  primaryMuts=merge(x=primaryMuts,y=mutationStatus,by.x="ID",by.y="mutationID")
  
  primaryMuts[,colnames(primaryMuts)%in% primarySamples]<-primaryMuts[,colnames(primaryMuts)%in% primarySamples]>=3
  
  
  for (s in colnames(primaryMuts)[colnames(primaryMuts)%in% primarySamples]){
    e=table(primaryMuts[,s],primaryMuts$status)
    e=e["TRUE",]/(e["TRUE",]+e["FALSE",])
    statMx[s,"Stem"]=e["Stem"]
    statMx[s,"Clade"]=e["Clade"]
  }
}

library(reshape)

rownames(statMx)<-c("290 (P1)","290 (P2)","219",
                    "308 (P1)",
                    "308 (P2)",
                    "323",
                    "328",
                    "DET52 (P1)",
                    "DET52 (P3)")

f=melt(statMx)
f=f[complete.cases(f),]
f$X2<-factor(f$X2,levels=c("Stem","Clade"))


pdf ("~/fig2d.pdf", height=4, width=4)
ggplot(f, aes(x=X1,y=value,fill=X2))+
  geom_bar(stat="identity")+
  facet_grid(~X2, scales="free_x", space = "free_x")+
  labs(x="Primary tumor",y="% deeply sequenced metastatic\nmutations detected in primary")+
  theme_bw()+
  scale_fill_manual(values=mutClassifierCols)+
  theme_classic(base_size = 11)+
  theme(legend.key = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background =element_rect(colour="white",fill="white"),
        strip.text = element_text(face = "bold", colour="black"),
        axis.title.x = element_text(face="bold"),
        axis.title.y = element_text(face="bold", margin=margin(0,10,0,0)),
        axis.text.x = element_text(angle = 45, hjust = 1, margin=margin(5,5,10,0)),
        legend.title = element_text(face="bold"),
        legend.position="bottom",
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1))+
  guides(fill=FALSE)

dev.off()
