# Load up DNA and RNA mutations from VCF and derive a VAF file

library(ggplot2)
library(ggrepel)
library(GGally)

mutClassifierCols<-c("darkblue", "tomato",  "darkgreen")


genesOfInterest <- scan("/gene lists/autopsy", what=character(),skip = 1)

muts.file<-"/mutations/exome/somatic/all-filtered.maf"


mutations <- read.table(muts.file, sep="\t",header=TRUE, fill=TRUE,comment.char="#",quote="",
                        blank.lines.skip=T, row.names = NULL, stringsAsFactors = FALSE)
mutations$patient<-sapply(strsplit(mutations$sample,"-"),"[",1)
mutations$ID=paste0(mutations$Chr,":",mutations$Start,"_",mutations$Ref_Allele,"/",mutations$Tumor_Alt_Allele)

suppressMessages(library(VariantAnnotation))

inputVCFDir="/mutations/exome/somatic/03a-lookback-vcf/"
vcfFiles <- list.files(inputVCFDir, pattern = "*.vcf")


.matrixOfListsToArray <- function(x) {
  # find number of elements of each cell of x
  n <- elementNROWS(x)
  maxn <- max(n)
  # for cells with less than the max number of elements, add NAs
  idx <- n < maxn
  x[idx] <- lapply(x[idx], function(a){c(a, rep(NA, maxn-length(a)))})
  # unlist and convert to array
  x <- array(unlist(x), dim=c(maxn, nrow(x), ncol(x)),
             dimnames=list(NULL, rownames(x), colnames(x)))
  x <- aperm(x, c(2,3,1))
  
  x
}

for (vcfFile in vcfFiles){
  
  cat (vcfFile,"\n")
  vcfPath <- paste0(inputVCFDir,vcfFile)
  patient <-  sapply(strsplit( vcfFile,"\\."),"[[",1)
  pt=sapply(strsplit( vcfFile,"-"),"[[",1)
  
  param <- ScanVcfParam()
  vcf <- readVcf(vcfPath, "hg19", param)
  vcf <- unique(vcf)
  normalSamples<-c("288-024","290-004","291-005","298-023","308-024","315-022","323-007","328-010","330-007","DET52-normal")
  nonBreast=c("290-022","298-004","298-005","298-006","298-007","298-012","298-020","298-021","298-022")
  
  vcf=vcf[, !grepl(paste(nonBreast,collapse = "|"),colnames(vcf)) ]
  vcf=vcf[, !grepl(paste(normalSamples,collapse = "|"),colnames(vcf)) ]
  
  samples=colnames(vcf)
  
  depth<-geno(vcf)$DP
  altCounts <- .matrixOfListsToArray(geno(vcf)$AD)[,,2]
  vaf<-altCounts/depth
  
  #classify mutations into Stem/clade/private
  noSamples=length(samples)
  vafCeil=ceiling(vaf)
  vafCeil[is.na(vafCeil)] <- 0
  vafCeil <- rowSums(vafCeil)
  noSamples=length(samples)
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
  
  
  
  vaf=merge(x=vaf,y=mutationStatus, by.x=0, by.y="mutationID")
  
  #is the mutation in a driver?
  mutations$patient=sapply(strsplit(mutations$patient,"\\."),"[",1)
  subMuts<-mutations[mutations$patient==pt,]
  vaf=merge(x=vaf,y=subMuts, by.x="Row.names", by.y="ID")
  vaf$Driver="Non-driver"
  vaf[vaf$Hugo_Symbol %in% genesOfInterest,"Driver"]<-"Driver"
  
  vaf[vaf$Driver!="Driver","Hugo_Symbol"]<-""
  vaf[vaf$status!="Stem","Hugo_Symbol"]<-""
  vaf$status<-factor(vaf$status, levels=c("Stem","Clade","Private"))
  
  colnames(vaf)
  vaf=vaf[,c(1,2:7,16:ncol(vaf))]
  e=data.frame(vaf)
  
  e=e[,c(1,grep("^X|^DE|sample|patient|status|Driver|Hugo_Symbol",colnames(e)))]
  
  fig3a=ggpairs(e, columns=c(2:length(grep("^X|^DE",colnames(e)))),diag = "blank", upper = "blank", 
                mapping=aes(color=status),
                lower = list(continuous = wrap("points", size=0.3), 
                             combo = wrap("dot", size=0.3) )
  )+
    ggtitle(pt)+
    theme_bw(base_size = 15)+
    theme(strip.background =element_rect(colour="white",fill="white"),
          strip.text = element_text(face = "bold", colour="black",margin = margin(0, 25, 0, 25, "cm")),
          plot.title =element_text(face="bold"),
          axis.title.x = element_text(face="bold"),
          axis.title.y = element_text(face="bold", margin=margin(0,10,0,0)),
          
          legend.title = element_text(face="bold"))

  
  for(i in 1:fig3a$nrow) {
    for(j in 1:fig3a$ncol){
      g=e
      if(unique(g$patient!="DET52")){
        g$sample=paste0("X",gsub("-","\\.",g$sample))
      }
      g=g[g$sample %in% c(fig3a$xAxisLabels[i], fig3a$yAxisLabels[j]),]
      g=g[!duplicated(g$Row.names),]
      
      fig3a[i,j] <- fig3a[i,j] + 
        scale_color_manual(name="Mutation Type:",values=mutClassifierCols) +
        geom_abline(slope=1, intercept = 0, col="gray10",size=0.3,linetype="dotted")+
        scale_x_continuous(limits = c(0, 1), breaks=seq(0,1,0.5))+
        scale_y_continuous(limits = c(0, 1), breaks=seq(0,1,0.5))+ 
        geom_text_repel(data=g, aes(label=Hugo_Symbol),color="black",
                        size = 3,
                        min.segment.length = unit(0, "lines"),
                        point.padding = unit(0.5, "lines"),
                        box.padding = unit(0.5, "lines"))
    }
  }
}

pdf("~/fig3a.pdf",height=9,width=9, useDingbats = F)
print(fig3a)  
dev.off()




