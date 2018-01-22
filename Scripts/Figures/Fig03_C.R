#get stem mutation clusters and plot

library(grid)
library(ggplot2)
library(reshape)

wd<-"/pyclone-exomes/pyclone-use/workspace/"
outDir<-"/pyclone-exomes/pyclone-use/plots/"


patients <- basename(list.dirs(wd, recursive = F))

cancerGenes <- scan("/gene lists/List_Drivers_NZ_Per.csv", what=character(),skip = 1)


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


#load mutation data
maf <- read.table("/mutations/exome/somatic/all-filtered.maf", sep="\t",header=TRUE, fill=TRUE,comment.char="#",quote="",
                  blank.lines.skip=T, row.names = NULL, stringsAsFactors = FALSE)
maf$patient=sapply(strsplit(maf$sample,"-"),"[",1)
maf$mutation_id=paste0(maf$Chr,":",maf$Start,"_",maf$Ref_Allele,"/",maf$Tumor_Alt_Allele)

patient=290

#load pyclone clusters, remove anything with a size less than 2 mutations
pycloneOutput <- paste0(wd,"/",patient,"/tables/")
py <- read.table(paste0(pycloneOutput,"cluster.tsv"), header = T, sep="\t",stringsAsFactors = F, colClasses = c("character","character",rep("numeric",3)))
py$cluster_lab <- paste0(py$cluster_id," (",py$size,")")
py <- py[py$size>2,]


#remove private mutation clusters
privateClusters <- cast(py, sample_id~cluster_id,value="mean")
rownames(privateClusters)<-privateClusters[,1]
privateClusters<-privateClusters[,-1]
clustersToKeep<-colnames(privateClusters)[apply(privateClusters, 2, function (x) sum(x>0.1))>2]
py<-py[py$cluster_id %in% clustersToKeep,]

if (nrow(py)>0){
  #load loci
  loci <- read.table(paste0(pycloneOutput,"loci.tsv"), header = T, sep="\t",stringsAsFactors = F, colClasses = c("character","character",rep("numeric",3)))
  loci$id=paste0(gsub(paste0(patient,"-"),"",loci$sample_id)," - ",loci$Site)
  loci <- loci[ loci$cluster_id %in% py$cluster_id,c(1:3)]
  
  #subset mutation data and retain mutations that have been clustered by PyClone
  mafSubset<- maf[maf$patient==patient,]
  mafSubset <-mafSubset[mafSubset$mutation_id %in% loci$mutation_id,]
  
  #select genes in cosmic
  inCosmic=merge(loci,maf,by="mutation_id")
  inCosmic=inCosmic[inCosmic$Hugo_Symbol %in% cancerGenes,]
  inCosmic<- inCosmic[,c(1,3,7)]
  inCosmic<-inCosmic[!duplicated(inCosmic),]
  inCosmic=merge(inCosmic,py,by.x="cluster_id",by.y="cluster_id")
  so=sampleOrder[[as.character(patient)]]
  
  
  g=inCosmic[,c(1,3)]  
  g=g[!duplicated(g),]
  g=aggregate(Hugo_Symbol ~., g, toString)
  g$Hugo_Symbol=paste0("atop(atop(textstyle('",g$Hugo_Symbol,"')))")
  g$Hugo_Symbol=gsub(", ","'),textstyle('", g$Hugo_Symbol)
  
  
  py$cluster_id<-as.numeric( py$cluster_id)
  py$label1=paste0("bold(",py$cluster_id,")")
  pyl=merge(py,g,by="cluster_id", all.x=T)  
  
  #remove patient ID
  if (patient!="DET52") {
    pyl$sample_id<-substring(as.character(pyl$sample_id),5)
    inCosmic$sample_id<-substring(as.character(inCosmic$sample_id),5)
    
  } else {
    pyl$sample_id<-substring(as.character(pyl$sample_id),7)
    inCosmic$sample_id<-substring(as.character(inCosmic$sample_id),7)
    
  }
  
  pyl$sample_id<-factor(pyl$sample_id,levels=so)
  inCosmic$sample_id<-factor(inCosmic$sample_id,levels=so)
  pyl$cluster_lab<-factor(pyl$cluster_lab,levels=unique(pyl$cluster_lab))
  
  
  
  pyl=pyl[pyl$mean>0.02,]
  
  cols<-c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6")
  fig3c=
    ggplot() +
    geom_bar(data=pyl, stat="identity", aes(x=sample_id, y=mean, fill=as.factor(cluster_lab)), 
             alpha=1,linetype="solid")+
    scale_x_discrete(name="Metastasis", expand = c(0.03, 0)) +
    scale_y_continuous(name="Mean Cellular Prevalence", expand=c(0.1,0),
                       breaks=seq(0,1,0.5))+
    ggtitle(patient)+
    theme_classic(base_size = 16)+
    theme(legend.key = element_blank(),
          strip.background =element_rect(colour="white",fill="white"),
          strip.text = element_text(face = "bold", colour="black",margin = margin(0, 25, 0, 25, "cm")),
          plot.title =element_text(face="bold"),
          axis.title.x = element_text(face="bold"),
          axis.title.y = element_text(face="bold", margin=margin(0,10,0,0)),
          axis.text.x = element_text(angle = 45, hjust = 1, margin=margin(5,5,10,0)),
          legend.title = element_text(face="bold"),
          legend.position="bottom",
          axis.line.x = element_line(color="black", size = 1),
          axis.line.y = element_line(color="black", size = 1))+
    guides(fill=F)+
    scale_fill_manual(values = cols)+
    facet_grid(cluster_id+Hugo_Symbol~., 
               labeller=labeller(.rows=label_parsed,.multi_line=TRUE))+
    theme(panel.border=element_blank(), axis.line=element_line()) + 
    theme(panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          panel.background=element_blank()) 
  
  pdf("~/fig3c.pdf",height=8,width=8, useDingbats = F, )
  print(fig3c)
  dev.off()
  
  
}







