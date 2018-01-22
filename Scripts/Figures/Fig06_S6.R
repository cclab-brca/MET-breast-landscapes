### Heatmpas based on Z scores of Immunophenograms ## 

let1 <- read.delim("..../IPS_autopsy_scores_matrix.txt",header=TRUE, sep="\t", check.names=FALSE) #  in TPM
head(let1)
rownames(let1) <- let1$GENE
let2 <- as.matrix(let1[,-1])
let2 <-let2[ ,1:5] # 288
let2 <-let2[ ,6:12] # 290 
let2 <-let2[ ,15:24] # 298
let2 <-let2[ ,25:42] # 308
let2 <-let2[ ,43:52] # 315
let2 <-let2[ ,55:57] # 328
let2 <-let2[ ,58:61]# 330
let2[1:7,1:7]
dim(let2)
head(let2)

# Heatmap3
install.packages("heatmap3")
library(heatmap3)

## Euclidean distance 
hc <- hclust(dist(t(let2), method = "euclidean"), method="ward.D")
summary(hr)
hr[3]
hr[4]
hr <- hclust(dist(let2, method = "euclidean"), method="complete")

# RColorBrewer
library(RColorBrewer)
my_palette1 <- colorRampPalette(c("purple", "black", "yellow"))(n = 299)

col_breaks = c(seq(-10,-0.5,length=100), 
               seq(-0.5,0.5,length=100),
               seq(0.5,10,length=100)) 

cols <- read.csv2(".../colorCode_RNA_61Samples_August2017_metastasisFinal.csv")
cols2 <- read.csv(".../rowside.csv")


pdf(file="../Heatmap.pdf", width=28, height=10)
heatmap3(as.matrix(let2), col=my_palette1, breaks=col_breaks, labR=, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),  main = "Patient ") 
legend("bottomright", inset =.00002, c("Metastasis site","Brain","Liver", "Ovary", "Gastric", "LN", "Meninges", "Lung.pleura", "Pericardium", "Bone", "Breast"), col=c("white","grey","darkblue","black", "snow", "gold",  "brown", "cyan4", "yellow2", "lightblue","deeppink1" ), pch=c(15,15,15,15,15,15,15,15,15,15,15,15), ncol=1, bty="n", lwd = 3)
legend("topleft", inset =.00002, c("Immunogenecity parameter","MHC Class I", "MHC Class II", "MHC - Non Class", "Check Points / immunomodulators" ,"Effector cells", "Supressor cells"), col=c("white","blue", "tomato", "grey", "gold", "green",   "purple" ), pch=c(15,15,15,15,15,15,15), ncol=1, bty="n", lwd = 3)
legend("bottomleft", inset =.00002, c("Weight","Positive", "Negative"), col=c("white","red", "black" ), pch=c(15,15,15), ncol=1, bty="n", lwd = 3)


dev.off()
