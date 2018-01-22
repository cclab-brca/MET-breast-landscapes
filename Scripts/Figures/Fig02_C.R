#### TS: mutation threshold, filtering drivers and non-drivers and Oncoprints from TS #### 


fileA<-".../308.500.sum.variants.txt"
fileB<-".../308.500.sum.variants.zscore.txt"

#convert A to a table where row = ID and column = sample
A <- read.table(fileA, header=T, stringsAsFactors = T, sep="\t",quote = NULL, row.names = NULL, check.names = FALSE)

rownames(A)<-A$ID
A <- A[,c(5:ncol(A))]
A[is.na(A)] <- 0
head(A)
dim(A)
View(A)

#convert B to a table where row = ID and column = sample
B <- read.table(fileB, header=T, stringsAsFactors = T,   sep="\t",quote = NULL, row.names = NULL, check.names = FALSE)
rownames(B)<-B[,1]
B <- B[,c(2:ncol(B))]
B[is.na(B)] <- 0
head(B)
dim(B)
View(B)

#make sure A and B have identical rows and columns
stopifnot(rownames(A)==rownames(B))
stopifnot(colnames(A)==colnames(B))


#keep a mutation if AF >= 0.01 and zscore >=3
af <- apply(A, 2, function(x) x>=.01)
zscore <-  apply(B, 2, function(x) x>=3)
mutationCalled <- af==T & zscore==T

#The table now contains TRUE/FALSE depending on whether a mutation was called or not based on these criteria
mutationCalled
library(reshape)
Z <- melt(mutationCalled)
head(Z)
dim(Z)
Z <- Z[Z$value=="TRUE",]
colnames(Z)[c(1,2 )]=c("ID","mets")
write.table(Z,file="/308_MutationsCalledCriteria.txt", append = FALSE, quote = FALSE, sep="\t", row.names = FALSE, col.names = TRUE)


##Merging for AA position and mutation class
B <- read.table("WESall-filtered.maf.txt", header=T, check.names = T)
B<- B[B$patient=="308",]
head(B)
unique(B$MAF_Variant)
dim(B)
#View(B)
A<- read.delim("/308_MutationsCalledCriteria.txt", header=T, sep="\t")
head(A)
dim(A)
#View(A)
mergedTest <- merge(A, B, all = T, by = c('ID'))
head(mergedTest)
dim(mergedTest)
View(mergedTest)
mergedTest <-mergedTest[mergedTest$mets!="NA",]
unique(mergedTest$mets)
write.table(mergedTest,file="../308_MutationsCalledCriteria_Merged.txt", append = FALSE, quote = FALSE, sep="\t", row.names = FALSE, col.names = TRUE)


# filtering drivers 
set_MB <-read.table("/308_MutationsCalledCriteria_Merged.txt", sep="\t",header=TRUE, fill=TRUE,comment.char="#",quote="",blank.lines.skip=T, row.names = NULL, stringsAsFactors = FALSE) 
set_jco <- read.csv2(".../ALL_637driver_genes.csv", header=TRUE) 
jco_subset <- set_MB[set_MB$Gene %in% set_jco[,,drop=T],] # ALL Drivers

write.table(jco_subset,file=".../308_Drivers.txt",append = FALSE, quote = FALSE, sep="\t", row.names = FALSE, col.names = TRUE)

# filtering non drivers
jco_subsetNOND <- set_MB[!set_MB$Gene %in% set_jco[,,drop=T],] # non Drivers
write.table(jco_subsetNOND,file=".../308_NON_Drivers.txt",append = FALSE, quote = FALSE, sep="\t", row.names = FALSE, col.names = TRUE)


# Casted tables drivers
jco_subset <- jco_subset[,c(1, 2, 4, 7, 5 )] # order  patient amplicon sample  gene zscore
head(jco_subset)
jco_subset1 <- cast(jco_subset, Mutation + ID ~ mets)
write.table(jco_subset1,file=".../308_Drivers_Casted.txt",append = FALSE, quote = FALSE, sep="\t", row.names = FALSE, col.names = TRUE)


# CAsted tables non drivers
jco_subsetNOND <- jco_subsetNOND[,c(1, 2, 4, 7, 5 )] # order  patient amplicon sample  gene zscore
jco_subset1 <- cast(jco_subsetNOND, Mutation + ID ~ mets)
write.table(jco_subset1,file=".../308_NON_Drivers_Casted.txt",append = FALSE, quote = FALSE, sep="\t", row.names = FALSE, col.names = TRUE)


##Oncoprint
#source("http://bioconductor.org/biocLite.R")
#biocLite("ComplexHeatmap")
#library(devtools)
#install_github("jokergoo/ComplexHeatmap")

library(ComplexHeatmap)
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(1, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
  },
  white = function(x, y, w, h) {
    grid.rect(x, y, w-unit(1, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "white", col = NA))
  },
  Missense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(1, "mm"), h*0.66, gp = gpar(fill = "#008000", col = NA))
  },
  TRUNC = function(x, y, w, h) {
    grid.rect(x, y, w-unit(1, "mm"), h*0.33, gp = gpar(fill = "black", col = NA))
  },
  Frame_Shift_Ins = function(x, y, w, h) {
    grid.rect(x, y, w-unit(1, "mm"), h*0.33, gp = gpar(fill = "black", col = NA))
  },
  Frame_Shift_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(1, "mm"), h*0.33, gp = gpar(fill = "black", col = NA))
  },
  Nonsense_Mutation  = function(x, y, w, h) {
    grid.rect(x, y, w-unit(1, "mm"), h*0.33, gp = gpar(fill = "black", col = NA))
  },
  Nonstop_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(1, "mm"), h*0.33, gp = gpar(fill = "black", col = NA))
  },
  Splice_Site = function(x, y, w, h) {
    grid.rect(x, y, w-unit(1, "mm"), h*0.33, gp = gpar(fill = "black", col = NA))
  },
  INFRAME   = function(x, y, w, h) {
    grid.rect(x, y, w-unit(1, "mm"), h-unit(0.4, "mm"), gp = gpar(fill = "brown", col = NA))
  },
  OTHER = function(x, y, w, h) {
    grid.rect(x, y, w-unit(1, "mm"), h*0.2, gp = gpar(fill = "purple", col = NA))
  },
  GERM = function(x, y, w, h) {
    grid.rect(x, y, w-unit(1, "mm"), h-unit(0.4, "mm"), gp = gpar(fill = "darkblue", col = NA))
  }
)
col = c("Missense_Mutation" = "#008000", "TRUNC" = "black", "Frame_Shift_Ins" = "black", "Frame_Shift_Del"= "black","Nonsense_Mutation"= "black", "Nonstop_Mutation"="black","Splice_Site"= "black",
        "INFRAME" = "brown",  "OTHER" = "purple","GERM" = "darkblue", "white"= "white", "Silent"= "purple")




X1 <- read.table("..../308_Drivers_Casted.txt", row.names = 1, header = TRUE, sep = "\t", stringsAsFactors = T,  check.names = FALSE)
head(X1)
X1 <- X1[, -c(1)]
X1 = (as.matrix(X1))
X1[1:4, 1:4]

#PLOT per patient
pdf(".../308_Drivers_final.pdf", width=25, height=20)
oncoPrint(X1, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col, 
          column_title = "308: Driver Mutations",
          show_column_names = TRUE,
          remove_empty_columns = FALSE,
          axis_gp = gpar(fontsize = 14),
          row_order = NULL, column_order = NULL,
          heatmap_legend_param = list(title = "Alterations", 
                                      at = c("Missense_Mutation", "TRUNC", "Frame_Shift_Ins", "Frame_Shift_Del","Nonsense_Mutation", "Nonstop_Mutation","Splice_Site","INFRAME",  "OTHER" ,"GERM"), 
                                      labels = c("Missense", "Truncating", "Truncating","Truncating","Truncating", "Truncating","Truncating","In_Frame", "Other", "Germline")))
dev.off()





X2 <- read.table("308_NON_Drivers_Casted.txt", row.names = 1, header = TRUE, sep = "\t", stringsAsFactors = T, check.names = FALSE)
head(X2)
X2 <- X2[, -c(1)]
X2 = (as.matrix(X2))

pdf(".../308_NON_Drivers_final.pdf", width=25, height=20)
oncoPrint(X2, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col, 
          column_title = " 308: Non Driver Mutations",
          show_column_names = T,
          remove_empty_columns = T,
          axis_gp = gpar(fontsize = 14),
          row_order = NULL, column_order = NULL ,
          column_names_gp = gpar(fontsize = 16),
          heatmap_legend_param = list(title = "Alterations", 
                                      at = c("Missense_Mutation", "TRUNC", "Frame_Shift_Ins", "Frame_Shift_Del","Nonsense_Mutation", "Nonstop_Mutation","Splice_Site","INFRAME",  "OTHER" ,"GERM", "*"), 
                                      labels = c("Missense", "Truncating", "Truncating","Truncating","Truncating", "Truncating","Truncating","In_Frame", "Other", "Germline",  "Manual Inspection"))) 



dev.off()



