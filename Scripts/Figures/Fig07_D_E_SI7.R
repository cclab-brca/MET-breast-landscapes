library(data.table)
library(stringr)
library(dplyr)
library(plyr)
library(scales)
library(vegan)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(dendextend)
library(igraph)
library(phangorn)
library(ape)



##############################################################
# common variables
##############################################################
#change this to the location of the data
data_file <- "./data/All_CDR3_sequences_REMERGED.csv"
meta_data_file <- "./data/Final_names_TCR_LDMA_23.05.2017.csv"


##############################################################
#helper function
##############################################################
ids2inc <- function(ids,  freq, levels, binary = T){
  retVec <- vector(mode = "numeric", length = length(levels))
  names(retVec) <- levels
  len <- length(ids)
  for( i in 1:len){
    if(binary){
      retVec[as.character(ids[i])] <- 1
    }
    else{
      retVec[as.character(ids[i])] <- freq [i]
    }
  }
  retVec
}


##############################################################
# prepare a data frame to hold all the data
# reading and data and some cleaning
##############################################################

cdrs <- read.csv(data_file, sep =",",stringsAsFactors = F)


#redundant header line
cdrs <- cdrs[cdrs$sample != "#sample",]

# remove these samples - not part of the study
# 70H_alpha 70H_beta
# 72H_alpha 72H_beta
# 76_alpha 76_beta
# 77_alpha 77_beta


cdrs <- subset(cdrs, sample != "70H_alpha" & sample != "70H_beta")
cdrs <- subset(cdrs, sample != "72H_alpha" & sample != "72H_beta")

cdrs <- subset(cdrs, sample != "76_alpha" & sample != "76__beta")
cdrs <- subset(cdrs, sample != "77_alpha" & sample != "77__beta")


# match metadata
metadata <- read.csv(meta_data_file, stringsAsFactors = F)

metadata <- subset(metadata, select = c(Old_sample1, chain, patient, Site, FINAL_ID_Metastasis_Code))
metadata <- dplyr::rename(metadata, Metastasis_Code = FINAL_ID_Metastasis_Code)
metadata <- dplyr::rename(metadata, sample = Old_sample1)


cdrs <- left_join(cdrs, metadata, by ="sample")

#before - 342208 sequences total
cdrs <- cdrs[str_detect(cdrs$CDR3_AA,"^C.*F$"),]
#after - 306407 sequences adhere to the C-F rule.

cdrs <- dplyr::rename(cdrs,  aminoAcid = CDR3_AA)
cdrs$count <- as.integer(cdrs$count)
cdrs <- as.data.table(cdrs)

# I also need to remove these samples
# 308-011, 308-012, 308-013

cdrs <- subset(cdrs, Metastasis_Code != "308-011")
cdrs <- subset(cdrs, Metastasis_Code != "308-012")
cdrs <- subset(cdrs, Metastasis_Code != "308-013")
######################################
#finished - constructing data frame
######################################


######################################
#Jaccard heatmap - beta chain
######################################

cdrs_beta <- subset(cdrs, chain == "beta") 

cdrs_beta[, aa_count := sum(count), by = c("patient","sample","aminoAcid")]
cdrs_beta <- unique(cdrs_beta, by =c("patient","sample","aminoAcid"))
cdrs_beta$count <- NULL
cdrs_beta <- dplyr::rename(cdrs_beta, count = aa_count)

metas <- unique(cdrs_beta$Metastasis_Code)

incMat_metas <- ddply(cdrs_beta, .(aminoAcid), function(df){ ids2inc(df$Metastasis_Code, df$count, metas, binary = TRUE) } )
incMat_metas$aminoAcid <- NULL
incMat_metas <- as.matrix(incMat_metas)

######jaccard - beta chain
res_jac_m <- vegan::vegdist(t(incMat_metas), method = "jaccard")
res_jac_m <- as.matrix(res_jac_m)
diag(res_jac_m)<- NA
res_jac_m <- 1-res_jac_m #distance to similarity

pheatmap(res_jac_m, color = brewer.pal(9,"Blues"),fontsize = 5, main = "Jaccard distances between metastases - beta chain")

######################################
#Jaccard heatmap - alpha chain
######################################

cdrs_alpha <- subset(cdrs, chain == "alpha") 
cdrs_alpha[, aa_count := sum(count), by = c("patient","sample","aminoAcid")]
cdrs_alpha <- unique(cdrs_alpha, by =c("patient","sample","aminoAcid"))  
cdrs_alpha$count <- NULL
cdrs_alpha <- dplyr::rename(cdrs_alpha, count = aa_count)

metas <- unique(cdrs_alpha$Metastasis_Code)

incMat_metas_a <- ddply(cdrs_alpha, .(aminoAcid), function(df){ ids2inc(df$Metastasis_Code, df$count, metas, binary = TRUE) } )
incMat_metas_a$aminoAcid <- NULL
incMat_metas_a <- as.matrix(incMat_metas_a)

######jaccard - alpha chain
res_jac_m <- vegan::vegdist(t(incMat_metas_a), method = "jaccard")
res_jac_m <- as.matrix(res_jac_m)
diag(res_jac_m)<- NA
res_jac_m <- 1-res_jac_m #distance to similarity

pheatmap(res_jac_m, color = brewer.pal(9,"Blues"),fontsize = 5, main = "Jaccard distances between metastases - alpha chain")

################################################
#now passing to the per patient analysis
################################################
#first analysis - compare for each patient the correlations 
#between the genetic tree, and the two tcr trees (alpha and beta)
#give a p-value for the correlation

################################################
#helper function - chain is alpha or beta
#returns an hclust object
################################################
clusterPatient <- function(patient_id, chain){
  
  if(chain == "beta"){
    metPat <- subset(cdrs_beta, patient == patient_id)
  }
  if(chain == "alpha"){
    metPat <- subset(cdrs_alpha, patient == patient_id)
  }
  #this deals with the gaps between OncoNem data and TCR data
  #not all samples have the two types of data
  if(patient_id == "288"){
    metPat <- subset(metPat, Metastasis_Code != "288-018")
  }
  
  if(patient_id == "290"){
    metPat <- subset(metPat, Metastasis_Code != "290-022")
  }
  
  if(patient_id == "298"){
    metPat <- subset(metPat, Metastasis_Code != "298-009")
    metPat <- subset(metPat, Metastasis_Code != "298-010")
  }
  
  metPat[,met_count := .N, by = "aminoAcid"]
  metPat <- subset(metPat, met_count >1)
  
  met_codes <- unique(metPat$Metastasis_Code)
  
  incMat_Pat <- ddply(metPat, .(aminoAcid), function(df){ ids2inc(df$Metastasis_Code, df$count, met_codes, binary = TRUE) } )
  incMat_Pat$aminoAcid <- NULL
  incMat_Pat <- as.matrix(incMat_Pat)
  
  ######jaccard
  jacPat <- vegan::vegdist(t(incMat_Pat), method = "jaccard")
  jacPat <- as.matrix(jacPat)
  diag(jacPat)<- NA
  jacPat <- 1-jacPat #distance to similarity
  
 
  #now pass to hclust object
  dPat <- dist(jacPat, method = "euclidean")
  clustsPat <- hclust(dPat, method = "complete")
 
  tree_title <- str_c("beta chain TCR tree patient ", patient_id, sep= "")
  plot(clustsPat, main = tree_title)

  clustsPat
}

#compute the hclust object for each patient - beta chain

clusts288 <- clusterPatient(288, "beta")
clusts290 <- clusterPatient(290, "beta")
clusts298 <- clusterPatient(298, "beta")
clusts308 <- clusterPatient(308, "beta")
clusts315 <- clusterPatient(315, "beta")

plot(clusts308,main =  "beta chain TCR tree - patient 308")

#compute the hclust object for each patient - alpha chain

clusts288_alpha <- clusterPatient(288, "alpha")
clusts290_alpha <- clusterPatient(290, "alpha")
clusts298_alpha <- clusterPatient(298, "alpha")
clusts308_alpha <- clusterPatient(308, "alpha")
clusts315_alpha <- clusterPatient(315, "alpha")

plot(clusts308_alpha,main =  "alpha chain TCR tree - patient 308")

#####################################################
#compare the immunological TCR trees as produced by the 
#clustering procedure to the genetic mutational trees of the
#oncoNEM model - my access to the model is by their RData files
#####################################################

#####################################################
#helper function
#####################################################
#####################################################
#computePvalImmunoGenetic
#function to compute a p value for an already computed correlation
#between two trees of type phylo
#input:
#phylo_tree1 - an object of type phylo
#phylo_tree2 - an object of type phylo
#trees_cor - computed correlation
#####################################################

computePvalPhyloTrees <- function(phylo_tree1, phylo_tree2, trees_cor ){
  
  numOfLabels <- length(phylo_tree1$tip.label)
  numOfRands <- 100
  random_phylos <- vector(mode ="list",length = numOfRands)
  
  for(i in 1:numOfRands){
    phylo_tree1_c <- phylo_tree1
    per <- sample(numOfLabels)
    phylo_tree1_c$tip.label <- phylo_tree1_c$tip.label[per]
    random_phylos[[i]] <- phylo_tree1_c
  }
  
  cop_cors <- vector(mode ="numeric",length = numOfRands)
  
  for(i in 1:numOfRands){
    cop_cors[i] <- cor_cophenetic(phylo_tree2,random_phylos[[i]])
  }
  
  p_val <- max(numOfRands - sum(abs(cop_cors) < trees_cor),1)/numOfRands
  
  p_val
}


#####################################################
#patient 308
#####################################################
load("./onconem_rdata/308_onconemNew.RData", verbose = T)
onconem308 <- onconem
## get samples of onconem tree and the node indices to which these samples belong.
samples <- unlist(onconem308$labels) ## all samples plus normal
samples <- samples[-which(samples == "N")]
samples <- samples[-which(samples == "024")]

nodeIndxPerSample <- sapply(samples,function(s) which(sapply(onconem308$labels, function(x) s%in%x)))

## calculate distance matrix nodes
dist <- shortest.paths(graph = onconem308$tree,v = V(onconem308$tree),to = V(onconem308$tree),weights = onconem308$edgeLength)
## select distances for each sample
dist <- dist[nodeIndxPerSample,nodeIndxPerSample]
rownames(dist) <- colnames(dist) <- str_c("308-", names(nodeIndxPerSample))

## create neighbour joining tree
genetic308 <- NJ(as.dist(dist))
genetic308 <- root(genetic308,outgroup = 1)

plot(genetic308, main = "Genetic tree - patient 308")
 

#compute cophenetic correlations
tcr308_beta <- as.phylo.dendrogram(clusts308)
genetic_beta308_cor <- cor_cophenetic(genetic308,tcr308_beta)
#0.6359549

pval_beta_genetic308 <- computePvalPhyloTrees(genetic308,tcr308_beta,genetic_beta308_cor)
#0.01

tcr308_alpha <- as.phylo.dendrogram(clusts308_alpha)
genetic_alpha308_cor <- cor_cophenetic(genetic308,tcr308_alpha)
#0.4663944

pval_alpha_genetic308 <- computePvalPhyloTrees(genetic308,tcr308_alpha,genetic_alpha308_cor)
#0.01

#correlation between alpha and beta - patient 308
ab_308_cor <- cor_cophenetic(tcr308_beta,tcr308_alpha)
#0.689649

pval_alpha_beta308 <- computePvalPhyloTrees(tcr308_beta,tcr308_alpha,ab_308_cor)
#0.01

###############################################################################
#patient 288
###############################################################################

#patient 288
load("./onconem_rdata/288_onconem.RData", verbose = T)
onconem288 <- onconem

samples <- unlist(onconem288$labels) ## all samples plus normal
samples <- samples[-which(samples == "N")]
samples <- samples[-which(samples == "021")]
samples <- samples[-which(samples == "022")]
samples <- samples[-which(samples == "023")]
samples <- samples[-which(samples == "024")]

nodeIndxPerSample <- sapply(samples,function(s) which(sapply(onconem288$labels, function(x) s%in%x)))

## calculate distance matrix nodes
dist <- shortest.paths(graph = onconem288$tree,v = V(onconem288$tree),to = V(onconem288$tree),weights = onconem288$edgeLength)
## select distances for each sample
dist <- dist[nodeIndxPerSample,nodeIndxPerSample]
rownames(dist) <- colnames(dist) <- str_c("288-", names(nodeIndxPerSample))

## create neighbour joining tree
genetic288 <- NJ(as.dist(dist))
genetic288 <- root(genetic288,outgroup = 1)

#plot(genetic288, main = "Genetic tree - patient 288")

#compute cophenetic correlations
tcr288_beta <- as.phylo.dendrogram(clusts288)
genetic_beta288_cor <- cor_cophenetic(genetic288,tcr288_beta)
#0.4463731

pval_beta_genetic288 <- computePvalPhyloTrees(genetic288,tcr288_beta,genetic_beta288_cor)
#0.04

tcr288_alpha <- as.phylo.dendrogram(clusts288_alpha)
genetic_alpha288_cor <- cor_cophenetic(genetic288,tcr288_alpha)
#0.5293715

pval_alpha_genetic288 <- computePvalPhyloTrees(genetic288,tcr288_alpha,genetic_alpha288_cor)
#0.01

#correlation between alpha and beta - patient 308
ab_288_cor <- cor_cophenetic(tcr288_beta,tcr288_alpha)
#0.9844122

pval_alpha_beta288 <- computePvalPhyloTrees(tcr288_beta,tcr288_alpha,ab_288_cor)
#0.03

######################################
# patient 290
######################################
load("./onconem_rdata/290exOvary_onconem.RData", verbose = T)

onconem290 <- onconem

samples <- unlist(onconem290$labels) ## all samples plus normal
samples <- samples[-which(samples == "N")]
samples <- samples[-which(samples == "016BWT")]
samples <- samples[-which(samples == "016BIDC")]
samples <- samples[-which(samples == "016Bmuc")]
samples <- samples[-which(samples == "004")]
samples[which(samples == "016A")] <- "016"

onconem290$labels[[6]] <- "016"

nodeIndxPerSample <- sapply(samples,function(s) which(sapply(onconem290$labels, function(x) s%in%x)))

## calculate distance matrix nodes
dist <- shortest.paths(graph = onconem290$tree,v = V(onconem290$tree),to = V(onconem290$tree),weights = onconem290$edgeLength)
## select distances for each sample
dist <- dist[nodeIndxPerSample,nodeIndxPerSample]
rownames(dist) <- colnames(dist) <- str_c("290-", names(nodeIndxPerSample))

## create neighbour joining tree
genetic290 <- NJ(as.dist(dist))
genetic290 <- root(genetic290,outgroup = 1)

#plot(genetic290, main = "Genetic tree - patient 290")

#compute cophenetic correlations
tcr290_beta <- as.phylo.dendrogram(clusts290)
genetic_beta290_cor <- cor_cophenetic(genetic290,tcr290_beta)
#-0.01708116

pval_beta_genetic290 <- computePvalPhyloTrees(genetic290,tcr290_beta,genetic_beta290_cor)
#1

tcr290_alpha <- as.phylo.dendrogram(clusts290_alpha)
genetic_alpha290_cor <- cor_cophenetic(genetic290,tcr290_alpha)
#0.006476021

pval_alpha_genetic290 <- computePvalPhyloTrees(genetic290,tcr290_alpha,genetic_alpha290_cor)
#0.97

#correlation between alpha and beta 
ab_290_cor <- cor_cophenetic(tcr290_beta,tcr290_alpha)
#0.247216

pval_alpha_beta290 <- computePvalPhyloTrees(tcr290_beta,tcr290_alpha,ab_290_cor)
#0.01

###################################
#patient 298
###################################

load("./onconem_rdata/298main_onconem.RData", verbose = T)

onconem298 <- onconem

samples <- unlist(onconem298$labels) ## all samples plus normal
samples <- samples[-which(samples == "N")]
samples <- samples[-which(samples == "020")]
samples <- samples[-which(samples == "023")]

nodeIndxPerSample <- sapply(samples,function(s) which(sapply(onconem298$labels, function(x) s%in%x)))

## calculate distance matrix nodes
dist <- shortest.paths(graph = onconem298$tree,v = V(onconem298$tree),to = V(onconem298$tree),weights = onconem298$edgeLength)
## select distances for each sample
dist <- dist[nodeIndxPerSample,nodeIndxPerSample]
rownames(dist) <- colnames(dist) <- str_c("298-", names(nodeIndxPerSample))

## create neighbour joining tree
genetic298 <- NJ(as.dist(dist))
genetic298 <- root(genetic298,outgroup = 1)

#plot(genetic298, main = "Genetic tree - patient 298")

#compute cophenetic correlations
tcr298_beta <- as.phylo.dendrogram(clusts298)
genetic_beta298_cor <- cor_cophenetic(genetic298,tcr298_beta)
#0.154957

pval_beta_genetic298 <- computePvalPhyloTrees(genetic298,tcr298_beta,genetic_beta298_cor)
#0.45

tcr298_alpha <- as.phylo.dendrogram(clusts298_alpha)
genetic_alpha298_cor <- cor_cophenetic(genetic298,tcr298_alpha)
#0.2520648

pval_alpha_genetic298 <- computePvalPhyloTrees(genetic298,tcr298_alpha,genetic_alpha298_cor)
#0.19

#correlation between alpha and beta 
ab_298_cor <- cor_cophenetic(tcr298_beta,tcr298_alpha)
#0.9665051

pval_alpha_beta298 <- computePvalPhyloTrees(tcr298_beta,tcr298_alpha,ab_298_cor)
#0.01

###########################################
#patient 315
###########################################
load("./onconem_rdata/315_onconem.RData", verbose = T)

onconem315 <- onconem

samples <- unlist(onconem315$labels) ## all samples plus normal
samples <- samples[-which(samples == "N")]
samples <- samples[-which(samples == "017")]
samples <- samples[-which(samples == "022")]

nodeIndxPerSample <- sapply(samples,function(s) which(sapply(onconem315$labels, function(x) s%in%x)))

## calculate distance matrix nodes
dist <- shortest.paths(graph = onconem315$tree,v = V(onconem315$tree),to = V(onconem315$tree),weights = onconem315$edgeLength)
## select distances for each sample
dist <- dist[nodeIndxPerSample,nodeIndxPerSample]
rownames(dist) <- colnames(dist) <- str_c("315-", names(nodeIndxPerSample))

## create neighbour joining tree
genetic315 <- NJ(as.dist(dist))
genetic315 <- root(genetic315,outgroup = 1)


#plot(genetic315, main = "Genetic tree - patient 315")

#compute cophenetic correlations
tcr315_beta <- as.phylo.dendrogram(clusts315)
genetic_beta315_cor <- cor_cophenetic(genetic315,tcr315_beta)
#0.202159

pval_beta_genetic315 <- computePvalPhyloTrees(genetic315,tcr315_beta,genetic_beta315_cor)
#0.1

tcr315_alpha <- as.phylo.dendrogram(clusts315_alpha)
genetic_alpha315_cor <- cor_cophenetic(genetic315,tcr315_alpha)
#0.3408417

pval_alpha_genetic315 <- computePvalPhyloTrees(genetic315,tcr315_alpha,genetic_alpha315_cor)
#0.01

#correlation between alpha and beta 
ab_315_cor <- cor_cophenetic(tcr315_beta,tcr315_alpha)
# 0.4360827

pval_alpha_beta315 <- computePvalPhyloTrees(tcr315_beta,tcr315_alpha,ab_315_cor)
#0.01


####################################################################
####################################################################
#second analysis - compare patient 308 (beta) to a random shuffle of its tcrs
#how likely are we to see the biological organization in a random data
####################################################################
####################################################################

####################
# helper function
####################

#clustering of a random repertoire - returns an hclust object
clusterRep <- function (repertorie){
  
  #cluster the sampeled repertoire
  
  mets_labels <- unique(repertorie$Metastasis_Code)
  incMat_per <- ddply(repertorie, .(aminoAcid), function(df){ ids2inc(df$Metastasis_Code, df$count, mets_labels, binary = TRUE) } )
  incMat_per$aminoAcid <- NULL
  incMat_per <- as.matrix(incMat_per)
  
  ######jaccard
  jacPer <- vegan::vegdist(t(incMat_per), method = "jaccard")
  jacPer <- as.matrix(jacPer)
  diag(jacPer)<- NA
  jacPer <- 1-jacPer #distance to similarity
  
  #pass to hclust object
  dPer <- dist(jacPer, method = "euclidean")
  clustsPer <- hclust(dPer)
  clustsPer
}


#this analysis is only for patient 308
met308 <- subset(cdrs_beta, patient == 308)
met308 <- as.data.table(met308)

tcr_total <- sum(met308$count)
all_tcrs <- rep(met308$aminoAcid, met308$count)
all_met_code <- rep(met308$Metastasis_Code, met308$count)
all_site <- rep(met308$Site, met308$count)

number_of_permutations <- 100
clust_beta <- vector(mode = "list", length = number_of_permutations)


for(i in 1:number_of_permutations){
  
  #get a permutaion
  permt <- sample.int(tcr_total)
  
  #shuffle data
  perm_tcr <- all_tcrs[permt]
  shuf_rep <- data.table(Metastasis_Code = all_met_code, Site = all_site, aminoAcid = perm_tcr)
  
  #compute aa_count
  shuf_rep[,count := .N, by = c("Metastasis_Code", "aminoAcid")]
  
  shuf_rep <- unique(shuf_rep, by = c("Metastasis_Code", "aminoAcid"))
  
  clust_beta[[i]] <- clusterRep(shuf_rep)
  
}

dendro_list <- lapply(clust_beta, as.dendrogram)
last <- length(dendro_list)+1
dendro_list[[last]] <- as.dendrogram(clusts308)

dendro_list <- as.dendlist(dendro_list)
names(dendro_list) <- c(paste("per", 1:(length(dendro_list)-1), sep =""),"patient308")

cop_cor <- cor.dendlist(dendro_list)

rand_bio_data <- cop_cor[last,1:(last-1)]

rand_rand_data <- cop_cor[1:(last-1),(1:last-1)]
rand_rand_data <- rand_rand_data[lower.tri(rand_rand_data, diag = FALSE)]

#statistical test
res.ks <- ks.test(rand_rand_data,rand_bio_data)

# Two-sample Kolmogorov-Smirnov test
# 
# data:  rand_rand_data and rand_bio_data
# D = 1, p-value < 2.2e-16
# alternative hypothesis: two-sided


#We also verified the result according to the Robinson-Foulds metric

phylo_list <- lapply(dendro_list, as.phylo.dendrogram)


#distance method is "PH85" (by default) - RF metric
distmat <- matrix(nrow = length(phylo_list),ncol = length(phylo_list))
rownames(distmat) <- names(phylo_list)
colnames(distmat) <- names(phylo_list)
for(i in 1:length(phylo_list)){
  for(j in 1:length(phylo_list)){
    distmat[i,j]<-dist.topo(phylo_list[[i]], phylo_list[[j]])
  }
}

#rf for Robinson and Foulds
rand_bio_data_rf <- distmat[last,1:(last-1)]

rand_rand_data_rf <- distmat[1:(last-1),(1:last-1)]
rand_rand_data_rf <- rand_rand_data_rf[lower.tri(rand_rand_data_rf, diag = FALSE)]

#statistical test
res_rf.ks <- ks.test(rand_rand_data_rf,rand_bio_data_rf)
#Two-sample Kolmogorov-Smirnov test
# data:  rand_rand_data_rf and rand_bio_data_rf
# D = 1, p-value < 2.2e-16
# alternative hypothesis: two-sided


