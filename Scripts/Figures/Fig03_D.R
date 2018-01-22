library(grid)
library(ggplot2)
library(reshape2)

wd<-"/pyclone-exomes/pyclone-use/workspace/"
outDir<-"/pyclone-exomes/pyclone-use/plots/"


raindanceDir<- "/targetted-seq/Multi-regional_20170627/seq_results/merge.freq/merge.calls.sum.500/"
patients <- basename(list.dirs(wd, recursive = F))

cancerGenes <- scan("/gene lists/List_Drivers_NZ_Per.csv", what=character(),skip = 1)
mafDir<-"/mutations/exome/somatic/05-maf/"


sampleOrder<-list()
sampleOrder[["288"]]<-c("005","006","022","008","016","017","020","021","023")
sampleOrder[["290"]]<-c("014","015","016-A","016-B-IDC","016-B-muc","016-B-WT","021","005","007","008","017","018","019","024","022")
sampleOrder[["291"]]<-c("015","019")
sampleOrder[["308"]]<-c("002","003","004","005","001","015","016","017","018","019","020","021","022","006","007","008","009","010","014")
sampleOrder[["323"]]<-c("003","004")
sampleOrder[["298"]]<-c("004","005","006","007","012","020","021","022","009","010")
sampleOrder[["315"]]<-c( "001","002","003","004","005","007","009","012","014","015","017","018")
sampleOrder[["DET52"]]<-c("mt3","mt4","mt5","mt6","mt7","mt8","t1","t2","t9")
sampleOrder[["328"]]<-c("003","004","005")
sampleOrder[["330"]]<-c("001","002","003","004","005")

normalSamples<-c("288-024","290-004","291-005","298-023","308-024","315-022","323-007","328-010","330-007","DET52-normal")


patient=290
pycloneOutput <- paste0(wd,"/",patient,"/tables/")

setwd(pycloneOutput)
clusters <- read.table(paste0(pycloneOutput,"cluster.tsv"), header = T, sep="\t",stringsAsFactors = F, colClasses = c("character","character",rep("numeric",3)))
clusters<-clusters[clusters$size>2,]
head(clusters)

samples<-unique(clusters$sample_id)

clustersToKeep<-unique(clusters$cluster_id)

#keep clusters with at least 3 mutations in them
py <- read.table(paste0(pycloneOutput,"loci.tsv"), header = T, sep="\t",stringsAsFactors = F, colClasses = c("character","character",rep("numeric",3)))
py<-py[py$cluster_id %in% clustersToKeep,]
py<-py[!py$cluster_id%in% c(11,12),]
py<-py[,c(1,3)]
py<-py[!duplicated(py),]

#have these mutations been detected in TS data?
raindanceFile<-  Sys.glob(file.path(paste0(raindanceDir,patient, "*")))

raindanceMuts<-read.table(raindanceFile,header=T, sep="\t",stringsAsFactors = F, quote ="\"")
head(raindanceMuts)
raindanceMuts$ID<-gsub("chr","",raindanceMuts$ID)

colnames(raindanceMuts)=gsub("X","",colnames(raindanceMuts))
colnames(raindanceMuts)=gsub("\\.","-",colnames(raindanceMuts))

raindanceMuts$ID<-paste0(sapply(strsplit(raindanceMuts$ID,":"),"[",1),"_",sapply(strsplit(raindanceMuts$ID,":"),"[",2))
chr<-sapply(strsplit(py$mutation_id,":"),"[",1)
loc<-sapply(strsplit(py$mutation_id,":"),"[",2)
loc<-sapply(strsplit(loc,"_"),"[",1)
py$ID<-paste0(chr,"_",loc)
py<-py[,-1]
head(py)

raindanceMuts=raindanceMuts[raindanceMuts$ID %in% py$ID ,]
raindanceMuts=raindanceMuts[,!colnames(raindanceMuts) %in% c("amplicon","mutation_id")]
nrow(raindanceMuts)

head(raindanceMuts)

if (patient=="DET52"){
  raindanceMuts=raindanceMuts[,-c(1,2,4:16)]
}

raindanceMuts=merge(raindanceMuts,py,by.x="ID",by.y="ID")
head(raindanceMuts)
raindanceMuts<-raindanceMuts[,!grepl("Control",colnames(raindanceMuts))]
raindanceMuts=melt(raindanceMuts,id.vars = c("ID","cluster_id"))
head(raindanceMuts)
raindanceMuts$variable<-as.character(raindanceMuts$variable)
raindanceMuts$type<-"Tumor"
raindanceMuts[raindanceMuts$variable %in% normalSamples,"type"]<-"Normal"

raindanceMuts[raindanceMuts$variable %in% normalSamples,"variable"]<-"Normal"

head(raindanceMuts)
raindanceMuts=raindanceMuts[complete.cases(raindanceMuts),]

raindanceMuts=raindanceMuts[raindanceMuts$variable %in% c("Normal",samples),]
numSamples=length(c("Normal",samples))
clusterWithMoreThanOneMut=names(table(raindanceMuts$cluster_id)[table(raindanceMuts$cluster_id)>numSamples])
raindanceMuts=raindanceMuts[raindanceMuts$cluster_id %in% clusterWithMoreThanOneMut,]

so=sampleOrder[[as.character(patient)]]
if (patient!="DET52") {
  raindanceMuts[raindanceMuts$variable=="Normal","variable"]<-"xxx-Normal"
  raindanceMuts$variable<-substring(as.character(raindanceMuts$variable),5)
} else {
  raindanceMuts[raindanceMuts$variable=="Normal","variable"]<-"xxxxx-Normal"
  raindanceMuts$variable<-substring(as.character(raindanceMuts$variable),7)
}

raindanceMuts$variable<-factor(raindanceMuts$variable,levels=c("Normal",so))

e=sort(unique(raindanceMuts$cluster_id))
raindanceMuts$cluster_id<-paste0("Cluster ",raindanceMuts$cluster_id)
raindanceMuts$cluster_id<-factor(raindanceMuts$cluster_id,levels=paste0("Cluster ",e))

cols<-c("#a6cee3","#33a02c","#fb9a99","#e31a1c","#cab2d6")

fig3d<-ggplot(raindanceMuts,aes(y=(value),x=variable,fill=as.factor(cluster_id)))+
  geom_boxplot(outlier.colour = "white")+
  labs(y="Deep sequencing mutation AF",x="Metastasis",title=patient)+
  facet_wrap(~as.factor(cluster_id), scales = "free_x")+
  #geom_hline(aes(yintercept = (0.01)),col="red", linetype="dotted")+
  theme_classic(base_size = 16)+
  theme(legend.key = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background =element_rect(colour="white",fill="white"),
        strip.text = element_text(face = "bold", colour="black"),
        plot.title =element_text(face="bold"),
        axis.title.x = element_text(face="bold"),
        axis.title.y = element_text(face="bold", margin=margin(0,10,0,0)),
        axis.text.x = element_text(angle = 45, hjust = 1, margin=margin(5,5,10,0), size = 12),
        legend.title = element_text(face="bold"),
        legend.position="bottom",
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1))+
  guides(fill=F)+
  scale_fill_manual(values = cols)

dev.off()
