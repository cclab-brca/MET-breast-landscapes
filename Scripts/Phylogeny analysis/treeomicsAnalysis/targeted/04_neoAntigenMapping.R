library(data.table)
scripdir='~/Documents/Caldas/Data/breastMetRepo/treeomicsAnalysis/'
source(paste0(scripdir,'03_extractTrees.R'))
source(paste0(scripdir,'04_plotting.R'))
source(paste0(scripdir,'05_collapseTree.R'))

datadir_NeoAntigens <- '~/Documents/Caldas/Data/dataRaw/04-GLD-OneHLA/'
datadir2 <- '~/Documents/Caldas/Data/dataProcessed/treeomics/'
## output directory
outdir <- '~/Documents/Caldas/Data/dataResults/treeomics/targeted/'
dir.create(paste0(outdir,'plots'),showWarnings = FALSE)
# list of driver mutations
drivers <- read.csv('~/Documents/Caldas/Data/dataAnnotation/List_Drivers_NZ_Per_Lef.csv',stringsAsFactors = FALSE)$Drivers

cases <- c(288,290,'290exOvary',291,298,'298main',308,315,323,328,330,'DET52') ## TODO 290,308,
SAMPLESETS <- c('All','Biopsies','Fluids','WESsamples')

for (case in cases) {
  if (case=='290exOvary') {
    baseCase=290
  } else if (case=='298Brain'|case=='298main') {
    baseCase=298
  } else if (case=='DET52main') {
    baseCase='DET52'
  } else {
    baseCase=case
  }
  
  ##### Step 1: read in neoantigen mutations and create single list for each case ##### 
  files <- list.files(datadir_NeoAntigens,pattern = as.character(baseCase))
  
  neoAntigens <- fread(paste0(datadir_NeoAntigens,files[1]),stringsAsFactors = FALSE,select = c('Chromosome','Start','Mutation','Gene.Name'),key = c('Chromosome','Start','Mutation','Gene.Name'))
  class(neoAntigens$Chromosome) <- "character"
  class(neoAntigens$Start) <- "numeric"
  
  ## remove duplicates
  neoAntigens <- neoAntigens[!duplicated(neoAntigens)]
  
  for (i in 2:length(files)) {
    dt <- fread(paste0(datadir_NeoAntigens,files[i]),stringsAsFactors = FALSE,select = c('Chromosome','Start','Mutation','Gene.Name'),key = c('Chromosome','Start','Mutation','Gene.Name'))
    class(dt$Chromosome) <- "character"
    class(dt$Start) <- "numeric"
    
    ## remove duplicates
    dt <- dt[!duplicated(dt)]
    
    neoAntigens <- merge(neoAntigens,dt,all=TRUE)
  }
  colnames(neoAntigens)[2] <- "Position"
  neoAntigens$neoAntigen <- 1
  
  for (sampleset in SAMPLESETS) {
    
    if (file.exists(paste0(datadir2,'/output_targeted_selected/extractedTrees/',case,'_',sampleset,'.txt'))) {
      
      ##### Step 2: map mutations to edges in trees #####
      dat <- scan(paste0(datadir2,'/output_targeted_selected/extractedTrees/',case,'_',sampleset,'.txt'),what='character')
      eTree <- extractTree(dat)
      
      edgeAttributes <- eTree$edgeAttributes
      edgeAttributes$Chromosome <- gsub("chr","",edgeAttributes$Chromosome)
      setkey(edgeAttributes,Chromosome,Position)
      
      # read annotation data
      dat <- fread(paste0(datadir2,'/input_targeted_selected/X',case,'_',sampleset,'_covCount.txt'),
                   select=c('Chromosome','Position','Gene'),
                   key = c('Chromosome','Position'),
                   colClasses = c(Chromosome="character",Position="numeric",Gene="character"))
      dat$Chromosome <- gsub("chr","",dat$Chromosome)
      edgeAttributes <- merge(edgeAttributes,dat)
      
      # merge edgeAttributes and neoAntigen data
      edgeAttributes <- merge(edgeAttributes,neoAntigens,all=TRUE)
      
      # # read posterior probabilities
      postFile <- list.files(path=paste0(datadir2,'/output_targeted_selected/X',case,'_',sampleset,'/'),pattern='posterior.txt',full.names = TRUE)
      posterior <- fread(postFile,sep='\t',stringsAsFactors = FALSE)
      ## remove 4 lines of parameters
      posterior <- posterior[-(1:4),]
      colnames(posterior)[2] <- "Position"
      setkey(posterior,Chromosome,Position)
      posterior$Keep <- apply(posterior[,-(1:4)],1,function(x) any(x>0.5))
      posterior <- posterior[,.(Chromosome,Position,Keep)]
      
      edgeAttributes <- merge(edgeAttributes,posterior,all=TRUE)
      edgeAttributes$neoAntigen[is.na(edgeAttributes$neoAntigen)] <- 2
      edgeAttributes$driver <- 2
      edgeAttributes$driver[edgeAttributes$Gene%in%drivers] <- 1
      setkey(edgeAttributes,neoAntigen,driver)
      
      edgeLengths <- data.frame(matrix(0,nrow = ecount(eTree$tree),ncol = 2))#4))
      colnames(edgeLengths) <- c("All","AllNeoAntigens")#,"Filtered","FilteredNeoAntigens")
      ## count edges
      edgeLengths$All <- table(factor(edgeAttributes$EdgeIdx,levels=1:ecount(eTree$tree)))
      edgeLengths$AllNeoAntigens <- table(factor(edgeAttributes[neoAntigen==1,EdgeIdx],levels=1:ecount(eTree$tree)))
      
      ##### check for convergent evolution ##### 
      duplicatedGenes <- edgeAttributes$Gene[duplicated(edgeAttributes$Gene)]
      edgeAttributes$geneAndEdge <- paste0(edgeAttributes$Gene,'_',edgeAttributes$EdgeIdx)
      convergentlyMutatedGenes <- sapply(duplicatedGenes, function(x) {
        if (length(unique(edgeAttributes$geneAndEdge[edgeAttributes$Gene==x]))>1) {
          x  
        } else {
          NA
        }})
      convergentlyMutatedGenes <- convergentlyMutatedGenes[!is.na(convergentlyMutatedGenes)]
      edgeAttributes$convergent <- edgeAttributes$Gene%in%convergentlyMutatedGenes
      
      
      ##### plotting #####
      ## plot with all mutations
      
      ## create edge labels for all mutations
      edgeLabelsNeo <- sapply(1:ecount(eTree$tree),function(x) edgeAttributes[EdgeIdx==x&neoAntigen==1,Gene])
      edgeLabelsAllOther <- sapply(1:ecount(eTree$tree),function(x) edgeAttributes[EdgeIdx==x&neoAntigen!=1,Gene])
      
      ## 
      edgeLabelsNeo_Driver <- lapply(1:length(edgeLabelsNeo), function(i) sapply(edgeLabelsNeo[[i]],function(x) {if (x%in%drivers) {x} else {""}} ))
      edgeLabelsNeo_NonDriver <- lapply(1:length(edgeLabelsNeo), function(i) sapply(edgeLabelsNeo[[i]],function(x) {if (!x%in%drivers) {x} else {""}} ))
      edgeLabelsAllOther_Driver <- lapply(1:length(edgeLabelsAllOther), function(i) sapply(edgeLabelsAllOther[[i]],function(x) {if (x%in%drivers) {x} else {""}} ))
      edgeLabelsAllOther_NonDriver <- lapply(1:length(edgeLabelsAllOther), function(i) sapply(edgeLabelsAllOther[[i]],function(x) {if (!x%in%drivers) {x} else {""}} ))
      
      ## collapse labels with line break
      edgeLabelsAllOther_Driver <- lapply(1:length(edgeLabelsNeo), function(i) paste(c(rep("",length(edgeLabelsNeo[[i]])),edgeLabelsAllOther_Driver[[i]]),collapse='\n'))
      edgeLabelsAllOther_NonDriver <- lapply(1:length(edgeLabelsNeo), function(i) paste(c(rep("",length(edgeLabelsNeo[[i]])),edgeLabelsAllOther_NonDriver[[i]]),collapse='\n'))
      edgeLabelsNeo_Driver <- lapply(edgeLabelsNeo_Driver, function(x) paste(x,collapse='\n'))
      edgeLabelsNeo_NonDriver <- lapply(edgeLabelsNeo_NonDriver, function(x) paste(x,collapse='\n'))
      
      
      ##### creat edge labels for convergent mutations #####
      convGenes <- edgeAttributes[convergent==TRUE]
      convGenes$ID <- paste0(convGenes$Gene,'_',convGenes$Position)
      convMutations <- lapply(1:ecount(eTree$tree), function(x) convGenes[EdgeIdx==x])
      
      convNeo_Driver <- lapply(1:length(convMutations),function(li) if (nrow(convMutations[[li]])>0) {sapply(1:nrow(convMutations[[li]]), function(mi) {
        if (convMutations[[li]][mi,neoAntigen]==1 & convMutations[[li]][mi,driver]==1) {
          convMutations[[li]][mi,ID]
        } else { "" }})} else {""})
      convNeo_NonDriver <- lapply(1:length(convMutations),function(li) if (nrow(convMutations[[li]])>0) {sapply(1:nrow(convMutations[[li]]), function(mi) {
        if (convMutations[[li]][mi,neoAntigen]==1 & convMutations[[li]][mi,driver]==2) {
          convMutations[[li]][mi,ID]
        } else { "" }})} else {""})
      convAllOther_Driver <- lapply(1:length(convMutations),function(li) if (nrow(convMutations[[li]])>0) {sapply(1:nrow(convMutations[[li]]), function(mi) {
        if (convMutations[[li]][mi,neoAntigen]==2 & convMutations[[li]][mi,driver]==1) {
          convMutations[[li]][mi,ID]
        } else { "" }})} else {""})
      convAllOther_NonDriver <- lapply(1:length(convMutations),function(li) if (nrow(convMutations[[li]])>0) {sapply(1:nrow(convMutations[[li]]), function(mi) {
        if (convMutations[[li]][mi,neoAntigen]==2 & convMutations[[li]][mi,driver]==2) {
          convMutations[[li]][mi,ID]
        } else { "" }})} else {""})
      
      convNeo_Driver <- lapply(convNeo_Driver, function(x) paste(x,collapse='\n'))
      convNeo_NonDriver <- lapply(convNeo_NonDriver, function(x) paste(x,collapse='\n'))
      convAllOther_Driver <- lapply(convAllOther_Driver, function(x) paste(x,collapse='\n'))
      convAllOther_NonDriver <- lapply(convAllOther_NonDriver, function(x) paste(x,collapse='\n'))
      
      ## collapse tree if any edges have length 0
      if (any(edgeLengths$All==0)) {
        collapsed <- collapseTree(tree = eTree$tree,
                                  e.length = edgeLengths$All,
                                  labels=eTree$labels)
      } else {
        E(eTree$tree)$idx <- 1:ecount(eTree$tree)
        collapsed <- list(tree=eTree$tree,
                          e.length = edgeLengths$All,
                          labels= eTree$labels)
      }
      
      ###### Plot Tree only #####
      width=8
      height=8
      h.adjust=NULL
      pdf(paste0(outdir,case,'_',sampleset,'_treeomics_treeOnly.pdf'),width = width,height = height)
      par(mar=c(1,4,0,3),xpd=TRUE)
      plotTree(tree = collapsed$tree,
               e.length = collapsed$e.length,
               clones=collapsed$labels,
               edge.arrow.mode='-',
               edge.label.color='black',
               axis=TRUE,
               v.label.cex = 1,
               edge.width=3,
               h.adjust = h.adjust)
      title(ylab="Accumulated mutations")
      dev.off()
      save(treeomicsTS=collapsed,file=paste0(outdir,case,'_',sampleset,'_treeomics_tree.RData'))
      
      ###### Plot Tree with edge annotation #####
      pdf(paste0(outdir,case,'_',sampleset,'_treeomics_tree_scaled.pdf'),width = 9+2*ecount(collapsed$tree),height = max(9,12*max(edgeLengths$All)/60))
      par(mfrow=c(1,2),mar=c(3,4,3,0))
      layout(matrix(c(1,2),ncol=2),widths = c(1,2),heights = c(1,1))
      plotTree(tree = collapsed$tree,
               e.length = collapsed$e.length,
               clones=collapsed$labels,
               edge.label=1:ecount(collapsed$tree),
               edge.arrow.mode='-',
               edge.label.color='black',
               axis=TRUE)
      title(main=paste("Case",case),ylab="Accumulated mutations")
      plot.new()
      n <- ecount(collapsed$tree)
      
      text(x=0,y=1.04,labels="All mutations",cex=0.75,adj=c(0,1),font=2)
      for (i in 1:n) {
        text(x=(i-1)/n,y=1.03,labels=i,cex=0.75,adj=c(0,1),font=2)
        text(x=(i-1)/n,y=1,labels=edgeLabelsNeo_Driver[[E(collapsed$tree)$idx[i]]],cex=0.65,adj=c(0,1),col='deeppink',font = 2)
        text(x=(i-1)/n,y=1,labels=edgeLabelsNeo_NonDriver[[E(collapsed$tree)$idx[i]]],cex=0.65,adj=c(0,1),col='deeppink')
        text(x=(i-1)/n,y=1,labels=edgeLabelsAllOther_Driver[[E(collapsed$tree)$idx[i]]],cex=0.65,adj=c(0,1),font = 2)
        text(x=(i-1)/n,y=1,labels=edgeLabelsAllOther_NonDriver[[E(collapsed$tree)$idx[i]]],cex=0.65,adj=c(0,1))
      }
      
      text(x=0,y=0.18,labels="Genes with multiple hits on different edges",cex=0.75,adj=c(0,1),font=2)
      for (i in 1:n) {
        text(x=(i-1)/n,y=0.16,labels=i,cex=0.75,adj=c(0,1),font=2)
        text(x=(i-1)/n,y=0.14,labels=convNeo_Driver[[E(collapsed$tree)$idx[i]]],cex=0.65,adj=c(0,1),col='deeppink',font = 2)
        text(x=(i-1)/n,y=0.14,labels=convNeo_NonDriver[[E(collapsed$tree)$idx[i]]],cex=0.65,adj=c(0,1),col='deeppink')
        text(x=(i-1)/n,y=0.14,labels=convAllOther_Driver[[E(collapsed$tree)$idx[i]]],cex=0.65,adj=c(0,1),font = 2)
        text(x=(i-1)/n,y=0.14,labels=convAllOther_NonDriver[[E(collapsed$tree)$idx[i]]],cex=0.65,adj=c(0,1))
      }
      dev.off()
      
      if (all(edgeLengths$AllNeoAntigens==0)) {
        print(paste("No neoantigen tree generated for case ",case," ",sampleset,"."))
        next()
      } else if (any(edgeLengths$AllNeoAntigens==0)) {
        collapsed <- collapseTree(tree = eTree$tree,
                                  e.length = edgeLengths$AllNeoAntigens,
                                  labels=eTree$labels)
      } else {
        E(eTree$tree)$idx <- 1:ecount(eTree$tree)
        collapsed <- list(tree=eTree$tree,
                          e.length = edgeLengths$AllNeoAntigens,
                          labels= eTree$labels)
      }
      
      ###### Plot Tree with neo-antigen causing mutation only #####
      if (baseCase%in%c(308,315)) {
        pdf(paste0(outdir,case,'_',sampleset,'_treeomics_tree_scaled_NeoAntigens.pdf'),width = 9+ecount(collapsed$tree)/2, height= 10)
      } else {
        pdf(paste0(outdir,case,'_',sampleset,'_treeomics_tree_scaled_NeoAntigens.pdf'),width = 9+ecount(collapsed$tree)/2, height= max(7,7*max(edgeLengths$AllNeoAntigens)/20))
      }
      par(mfrow=c(1,2),mar=c(3,4,3,0))
      plotTree(tree = collapsed$tree,
               e.length = collapsed$e.length,
               clones=collapsed$labels,
               edge.label=1:ecount(collapsed$tree),
               edge.arrow.mode='-',
               edge.label.color='black',
               axis=TRUE)
      title(main=paste("Case",case),ylab="Accumulated neo-antigen mutations")
      plot.new()
      n <- ecount(collapsed$tree)
      for (i in 1:ceiling(n/2)) {
        text(x=2*(i-1)/n,y=1.03,labels=i,cex=0.75,adj=c(0,1),font=2)
        # text(x=2*(i-1)/n,y=1,labels=edgeLabelsNeo[[E(collapsed$tree)$idx[i]]],cex=0.65,adj=c(0,1),col='deeppink')
        text(x=2*(i-1)/n,y=1,labels=edgeLabelsNeo_Driver[[E(collapsed$tree)$idx[i]]],cex=0.65,adj=c(0,1),col='deeppink',font = 2)
        text(x=2*(i-1)/n,y=1,labels=edgeLabelsNeo_NonDriver[[E(collapsed$tree)$idx[i]]],cex=0.65,adj=c(0,1),col='deeppink')
      }
      if (n>1) {
        for (i in (ceiling(n/2)+1):n) {
          text(x=2*(i-1-ceiling(n/2))/n,y=0.53,labels=i,cex=0.75,adj=c(0,1),font=2)
          # text(x=2*(i-1-ceiling(n/2))/n,y=0.5,labels=edgeLabelsNeo[[E(collapsed$tree)$idx[i]]],cex=0.65,adj=c(0,1),col='deeppink')
          text(x=2*(i-1-ceiling(n/2))/n,y=0.5,labels=edgeLabelsNeo_Driver[[E(collapsed$tree)$idx[i]]],cex=0.65,adj=c(0,1),col='deeppink',font = 2)
          text(x=2*(i-1-ceiling(n/2))/n,y=0.5,labels=edgeLabelsNeo_NonDriver[[E(collapsed$tree)$idx[i]]],cex=0.65,adj=c(0,1),col='deeppink')
        }
      }
      dev.off()
      
      

    }
  }
}
