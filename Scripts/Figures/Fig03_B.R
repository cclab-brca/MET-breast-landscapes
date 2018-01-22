suppressMessages(library(VariantAnnotation))

inputVCFDir="/mutations/exome/somatic/03b-lookback-vcf-normfilter/"
outDir <- "/plots/"
vcfFiles <- list.files(inputVCFDir, pattern = "*.vcf")
vcfFiles=vcfFiles[!grepl("290-004ov.combined.vcf|298-023en.combined.vcf",vcfFiles)]

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

load("ccfMuts.RData")
res$patient<-sapply(strsplit(res$sample,"-"),"[",1)

vcfFile=vcfFiles[2]

cat (vcfFile,"\n")
vcfPath <- paste0(inputVCFDir,vcfFile)
patient <-  sapply(strsplit( vcfFile,"\\."),"[[",1)
pt=sapply(strsplit( vcfFile,"-"),"[[",1)

param <- ScanVcfParam()
vcf <- readVcf(vcfPath, "hg19", param)

vcf <- unique(vcf)

rd=c("290-022","298-004","298-005","298-006","298-007","298-012","298-020","298-021","298-022",
     "288-024","290-004","291-005","298-023","308-024","315-022","323-007","328-010","330-007","DET52-normal")
vcf=vcf[, !grepl(paste(rd,collapse = "|"),colnames(vcf)) ]

#get sample names
samples=colnames(vcf)


#classify mutations
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
res=res[res$patient==sapply(strsplit(patient,"-"),"[",1),]
res$chr=gsub("chr","",res$chr)
res$mutationID<-paste0(res$chr,":",res$start,"_",res$ref,"/",res$alt)
f=merge(x=res,y=mutationStatus, by.x="mutationID", by.y="mutationID")
f=f[,c(1,7,20:24,27,28)]
f=f[complete.cases(f),]
f$color="1"
f[f$ccf_upper==1,"color"]<-"2"

f=f[order(f$color,f$status,f$ccf_point, f$ccf_lower, f$ccf_upper,decreasing = F),]

f=f[f$status!="Private",]

for (s in unique(f$sample)){
  
  for (d in c("Stem","Clade")){
    noSamp<-(nrow(f[f$patient==sapply(strsplit(patient,"-"),"[",1) & f$status==d & f$sample==s,]))
    f[f$status==d & f$sample==s,"loc"]<- c( 1 : noSamp)
    
  }
  
}
f$status<-factor(f$status,levels=c("Stem","Clade","Private"))

head(f)
pt=sapply(strsplit(patient,"-"),"[",1)

if (pt!="DET52") {
  f$sample<-substring(as.character(f$sample),5)
} else {
  f$sample<-substring(as.character(f$sample),7)
}




cols <- c("#e31a1c","#1f78b4")
pdf("~/fig3b.pdf",width=8,height=8)

f=f[f$sample!="022" & f$sample!="016-B-muc" & f$sample!="016-B-WT" & f$sample!="016-B-IDC",]

ggplot(f,aes(x=loc, color=color))+
  geom_pointrange(aes(y=ccf_point,ymin = ccf_lower, ymax = ccf_upper), size=0.1,position = "jitter")+
  facet_grid(sample~status, space = "free", scales = "free")+
  labs(x="Mutations",y="Cancer Cell Fraction")+
  theme_classic(base_size = 14)+
  theme(legend.key = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background =element_rect(colour="white",fill="white"),
        strip.text = element_text(face = "bold", colour="black"),
        plot.title =element_text(face="bold"),
        axis.title.y = element_text(face="bold", margin=margin(0,10,0,0)),
        axis.title.x=element_text(face="bold"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_text(face="bold"),
        legend.position="bottom",
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1))+
  scale_y_continuous(breaks=c(0,0.5,1))+
  scale_color_manual(name="",values=cols,breaks=c(1,2), labels=c("Subclonal","Clonal"))+
  guides(fill=F)
dev.off()



