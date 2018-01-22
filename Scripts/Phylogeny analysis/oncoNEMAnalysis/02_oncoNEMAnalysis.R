## oncoNEM on treeomics posteriors
## scripts
source('~/Documents/Caldas/Data/breastMetRepo/treeomicsAnalysis/04_plotting.R')
source('~/Documents/Caldas/Data/breastMetRepo/treeomicsAnalysis/05_collapseTree.R')
## required input data
# output of neoantigen analysis
datadir_NeoAntigens <- '~/Documents/Caldas/Data/dataRaw/04-GLD-OneHLA/'
# output of treeomics WES analysis (posteriors)
treeomicsDir_WES <- '~/Documents/Caldas/Data/dataProcessed/treeomics/output_WES_all/'
# output of treeomics targeted analysis (posteriors)
treeomicsDir_targeted <- '~/Documents/Caldas/Data/dataProcessed/treeomics/output_targeted_all'
# raw data dir (used for filtering)
rawdatadir <- '~/Documents/Caldas/Data/dataRaw/targetedSeq/merge.calls.sum.500-exomeSelect-normalAsControl/'
# list of driver mutations
drivers <- read.csv('~/Documents/Caldas/Data/dataAnnotation/List_Drivers_NZ_Per_Lef.csv',stringsAsFactors = FALSE)$Drivers
## output directory
outdir <- '~/Documents/Caldas/Data/dataResults/oncoNEM'

## required packages
library(oncoNEM)
library(ggplot2)
library(data.table)
library(igraph)
library(heatmap.plus)
library( 'tikzDevice' )
################################################################################
## helper functions
## Orders edges of tree from top to bottom
orderEdges <- function(tree) {
  nodeOrder <- orderNodes(tree)
  el <- get.edgelist(tree)
  edgeOrder <- rep(0,nrow(el))
  for (i in 1:nrow(el)) {
    edgeOrder[i] <- which(nodeOrder[-1]==el[i,2])
  }
  return(edgeOrder)
}

## orders nodes by branches for nicer heatmap plotting
orderNodes <- function(tree) {
  leaves <- (2:vcount(tree))[degree(tree)[-1]==1]
  nodeOrder <- 1
  currLeaf <- leaves[1]
  while (length(leaves) > 0) {
    leaves <- leaves[leaves!=currLeaf]
    sPath <- get.shortest.paths(tree,from=1,to=currLeaf)$vpath[[1]]
    nodeOrder <- c(nodeOrder,sPath[!sPath%in%nodeOrder])
    currLeaf <- leaves[which.min(shortest.paths(tree,v=currLeaf,to=leaves))]
  }
  return(nodeOrder)
}

## Plot heatmap of mutations
plothm <- function(data,rowGroup,colorsOrigin,rowColors,containsAmbiguous,...) {
  rowColors <- cbind(rowColors,rowColors)
  colnames(rowColors) <- NULL
  heatmap.plus(t(data[order(rowGroup),]),Rowv = NA, Colv = NA, scale = 'none',
               col=colorRampPalette(c('red','white','blue'))(100),
               ColSideColors = cbind(colorsOrigin[sort(rowGroup)],
                                     colorsOrigin[sort(rowGroup)]),
               RowSideColors = rowColors,
               margins = c(max(nchar(rownames(data))) + 8,
                           max(nchar(colnames(data))) + 1),#c(15,15),
               ...)
  if (containsAmbiguous) {
    legend=c(1:(max(rowGroup)-1),'ambiguous')
  } else {
    legend=1:max(rowGroup)
  }
  par(xpd=TRUE)
  ## choose number of columns for legend
  if (max(rowGroup) > 16) {
    ncol=ceiling(max(rowGroup)/3)
  } else if (max(rowGroup) > 8) {
    ncol=ceiling(max(rowGroup)/2)
  } else {
    ncol=max(rowGroup)
  }
  legend(x = 0.05, y = 0,
         title = 'Edge assignment',
         title.adj = 0,
         legend = legend,
         col=colorsOrigin[1:max(rowGroup)],
         pch = 15,
         pt.cex = 1.5,
         bty = 'n',
         ncol = ncol)
}

## color vectors
colorsOrigin <- c('black','darkcyan','magenta4','deeppink2','dodgerblue',
                  'blue4','cyan3','blue','indianred2','cornflowerblue',
                  'yellow3','slategrey','chartreuse3','mediumturquoise',
                  'thistle',rainbow(10))

## plot heatmap of likelihood for parameter estimataion
plotLLH <- function(llh,test.fpr,test.fnr,x=200) {
  ## get estimated parameters
  indx <- which(llh==max(llh), arr.ind=TRUE)
  fpr <- test.fpr[indx[1]]
  fnr <- test.fnr[indx[2]]
  
  ## reformat data for plotting
  xy <- expand.grid(1:length(test.fpr),1:length(test.fnr))
  llhdf <- data.frame(fpr=test.fpr[xy[,1]],
                      fnr=test.fnr[xy[,2]],
                      llh=apply(xy,1,function(x) llh[x[1],x[2]]))
  
  ## find range of llh values to plot
  x <- min(x,max(llh)+quantile(-llh,probs = 0.8))
  
  ## generate plot
  gg <- ggplot(llhdf,aes(x=fpr,y=fnr)) +
    geom_tile(aes(fill = llh),colour = "white") +
    scale_fill_gradientn(colors=rainbow(7),
                         limits=c(max(llhdf$llh)-x,max(llhdf$llh))) +
    geom_point(data=data.frame(fpr=fpr,fnr=fnr),
               aes(x=fpr,y=fnr),shape=24,
               colour="black",
               bg="white",
               size=2)
  print(gg)
}
################################################################################

# TYPES <- c('WES','targeted') #TODO
TYPES <- 'WES'
for (type in TYPES) {
  
  if (type=='WES') {
    SAMPLESETS <- 'All'
    CASES <- c(288,'290exOvary',291,'298Brain','298main','298All',308,315,'315_018',323,328,330,'DET52','DET52main','DET52main_mt3','DET52main_mt3mt4','DET52_mt3','DET52_mt3mt4') 
  } else {
    SAMPLESETS <- c('All','Biopsies','Fluids','WESsamples')
    CASES <- c(288,290,'290exOvary',291,298,'298main',308,315,323,328,330,'DET52')
    
  }
  
  dir.create(paste0(outdir,'/',type),showWarnings = FALSE, recursive = TRUE)
  
  for (case in CASES) {
    
    if (case=='290exOvary') {
      baseCase=290
    } else if (case=='298Brain'|case=='298main'|case=='298All') {
      baseCase=298
    } else if (case=='DET52main'|case=='DET52main_mt3'|case=='DET52main_mt3mt4'|case=='DET52_mt3'|case=='DET52_mt3mt4') {
      baseCase='DET52'
    } else if (case=='315_018') {
      baseCase=315
    } else {
      baseCase=case
    }
    
    for (sampleSet in SAMPLESETS) {
      if (type=='WES') {
        outName <- paste0(outdir,'/WES/',case,'_oncoNEM_')
      } else {
        outName <- paste0(outdir,'/targeted/',case,'_',sampleSet,'_oncoNEM_')
      }
      
      ############################ data preprocessing  #############################
      ######## read treeomics posterior
      if (type=='WES') {
        if (case=='315_018'|case=='DET52_mt3'|case=='DET52_mt3mt4') {
          postFile <- list.files(path=paste0(treeomicsDir_WES,'/X',baseCase),pattern='posterior.txt',full.names = TRUE)
        } else if (case=='DET52main_mt3' | case == 'DET52main_mt3mt4') {
          postFile <- list.files(path=paste0(treeomicsDir_WES,'/XDET52main'),pattern='posterior.txt',full.names = TRUE)
        } else {
          postFile <- list.files(path=paste0(treeomicsDir_WES,'/X',case),pattern='posterior.txt',full.names = TRUE)
        }
      }
      if (type=='targeted') {
        postFile <- list.files(path=paste0(treeomicsDir_targeted,'/X',case,'_',sampleSet),pattern='posterior.txt',full.names = TRUE)
      }
      
      if (length(postFile)==0) {
        print(paste("Post file does not exist for case", case, type ,sampleSet))
        next
      } else if (length(postFile)>1) {
        print(paste("More than one post files detected for case", case, type ,sampleSet))
        next
      }
      posterior <- read.csv(postFile,sep='\t',stringsAsFactors = FALSE)
      if (case=='315_018') {
        posterior <- posterior[,-grep('X018',colnames(posterior))]
      } else if (case=='DET52_mt3'|case=='DET52main_mt3') {
        posterior <- posterior[,-grep('Xmt3',colnames(posterior))]
      } else if (case=='DET52_mt3mt4'|case=='DET52main_mt3mt4') {
        posterior <- posterior[,-c(grep('Xmt3',colnames(posterior)),
                                   grep('Xmt4',colnames(posterior)))]
      }
      genes <- posterior[-c(1:4),c(1,2,4)]
      colnames(genes) <- c('Chromosome','Position','Gene.Name')
      
      ## remove rows (purities, VAFs and priors) and columns (Chr, Start, End, Gene) 
      posterior <- posterior[-(1:4),-(1:4)]
      posterior <- apply(posterior,2,as.numeric)
      
      ## binarize data
      binary <- posterior
      binary[binary>=0.5] <- 1 
      binary[binary<0.5] <- 0
      colnames(binary) <- colnames(posterior)
      
      if (type=='targeted') {
        ## only use sites that are in read depth filtered file.
        
        ## get index of these sites
        filteredFile <- list.files(path=rawdatadir,pattern=paste0(baseCase,'.500..*sum.variants.txt'),full.names = TRUE) 
        filtered <- fread(filteredFile)
        genes$ID <- paste0(genes$Chromosome,":",genes$Position)
        filterIdx <- genes$ID%in%filtered$ID
        
        ## remove all other sites from data
        genes <- genes[filterIdx,]
        binary <- binary[filterIdx,]
        posterior <- posterior[filterIdx,]
      }      
      
      ## filter: remove sites without mutation in at least one sample
      idx <- apply(binary,1,function(x) sum(x==1))!=0 & !is.na(apply(binary,1,sum))
      binary <- binary[idx,]
      genes <- genes[idx,]
      genes$order <- 1:nrow(genes)
      if (type=='targeted') {
        genes$Gene.Name <- gsub("([^\\_]*)\\_.*","\\1",genes$Gene.Name)
      }
      genes$ID <- paste0(genes$Gene.Name,'_',genes$Position)
      posterior <- posterior[idx,]
      rownames(posterior) <- genes$Gene.Name
      colnames(posterior) <- gsub("X\\_?","",colnames(posterior))
      
      ######## read neoantigen data
      files <- list.files(datadir_NeoAntigens,pattern = as.character(baseCase))
      
      neoAntigens <- fread(paste0(datadir_NeoAntigens,files[1]),stringsAsFactors = FALSE,select = c('Chromosome','Start','Gene.Name'),key = c('Chromosome','Start','Gene.Name'))
      class(neoAntigens$Chromosome) <- "character"
      class(neoAntigens$Start) <- "numeric"
      
      ## remove duplicates
      neoAntigens <- neoAntigens[!duplicated(neoAntigens)]
      
      for (i in 2:length(files)) {
        dt <- fread(paste0(datadir_NeoAntigens,files[i]),stringsAsFactors = FALSE,select = c('Chromosome','Start','Gene.Name'),key = c('Chromosome','Start','Gene.Name'))
        class(dt$Chromosome) <- "character"
        class(dt$Start) <- "numeric"
        
        ## remove duplicates
        dt <- dt[!duplicated(dt)]
        
        neoAntigens <- merge(neoAntigens,dt,all=TRUE)
      }
      colnames(neoAntigens)[2] <- "Position"
      neoAntigens$neoAntigen <- 1
      
      ############################# oncoNEM Analysis  ##############################
      
      #### optimize error rates #####
      test.fpr <- seq(from = 0.0001, to = 0.1, length.out = 10)
      test.fnr <- seq(from = 0.0001, to = 0.1, length.out = 10)
      if (!file.exists(paste0(outName,'llh.RData'))) {
        llh <- matrix(0,nrow = length(test.fpr),ncol = length(test.fnr))
        for (i.fpr in 1:length(test.fpr)) {
          for (i.fnr in 1:length(test.fnr)) {
            ## initialize oncoNEM
            oNEM <- oncoNEM$new(Data = binary,
                                FPR = test.fpr[i.fpr],
                                FNR = test.fnr[i.fnr])

            ## run initial search
            oNEM$search(delta = 100)
            oNEM.expanded <- expandOncoNEM(oNEM,epsilon = 1,delta = 100,
                                           checkMax = 10000,app = TRUE)

            ## save log-likelihood of highest-scoring tree
            llh[i.fpr,i.fnr] <- oNEM.expanded$best$llh
          }
        }

        save(llh,file=paste0(outName,'llh.RData'))
      } else {
        load(paste0(outName,'llh.RData'))
      }
      pdf(paste0(outName,'llh.pdf'),
          width=5,height=4)
      plotLLH(llh = llh,test.fpr = test.fpr,test.fnr = test.fnr)
      dev.off()

      indx <- which(llh==max(llh), arr.ind=TRUE)
      fpr.est <- test.fpr[indx[1]]
      fnr.est <- test.fnr[indx[2]]
      
      # fpr.est=0.025
      # fnr.est=0.025
      
      ##### run oncoNEM ##### 
      oNEM <- oncoNEM$new(Data = binary,
                          FPR = fpr.est,
                          FNR = fnr.est)
      # run initial search until best tree has not changed for 10 steps
      oNEM$search(delta = 100)
      oNEM.expanded <- expandOncoNEM(oNEM,epsilon = 2, delta = 200,checkMax = 40000)
      oNEMclust <- clusterOncoNEM(oNEM.expanded, epsilon = 2)
      
      clones <- lapply(oNEMclust$clones, function(x) colnames(posterior)[x])
      clones[[1]] <- c("N",clones[[1]])
      
      ## calculate posterior and edge lengths
      post <- oncoNEMposteriors(tree = oNEMclust$g,
                                clones = oNEMclust$clones,
                                Data = binary,
                                FPR = fpr.est,
                                FNR = fnr.est)
      theta <- post$p_theta
      theta <- apply(theta,1,function(x) which(x==max(x)))
      idx <- sapply(theta,length)
      eLength <- table(factor(apply(post$p_theta[idx==1,-1,drop=FALSE],1,which.max),levels=1:ecount(oNEMclust$g)))
      
      if (any(eLength==0)) {
        ## collapse Tree (if mutations are assigned to multiple nodes it can happen that some branches are not assigned a single mutation on their own)
        collapsed <- collapseTree(oNEMclust$g,e.length = eLength,labels=1:vcount(oNEMclust$g))
        collapsed$clones <- lapply(collapsed$labelsSplit, function(i) unlist(clones[i]))
        collapsedClonesNumeric <- lapply(collapsed$labelsSplit, function(i) unlist(oNEMclust$clones[i]))
        
        ## reorder edges, otherwise there are problems with the index of theta
        el <- get.edgelist(collapsed$tree)
        el <- el[order(el[,2]),]
        collapsed$tree <- graph.edgelist(el)
        
        ## recalculate posterior and edgeLengths
        post <- oncoNEMposteriors(tree = collapsed$tree,
                                  clones = collapsedClonesNumeric,
                                  Data = binary,
                                  FPR = fpr.est,
                                  FNR = fnr.est)
        theta <- post$p_theta
        theta <- apply(theta,1,function(x) which(x==max(x)))
        idx <- sapply(theta,length)
        eLength <- table(factor(apply(post$p_theta[idx==1,-1,drop=FALSE],1,which.max),levels=1:ecount(collapsed$tree)))
        collapsed$e.length <- eLength
      } else {
        collapsed <- list(tree=oNEMclust$g,
                          e.length=eLength,
                          clones=clones)
      }
      save(collapsed,file=paste0(outName,'tree.RData'))
      
      ##### assign mutations to edges ##### 
      assignedMutations <- theta
      assignedMutations[idx!=1] <- max(unlist(theta)) + 1 ## ambiguously assigned mutations 
      assignedMutations <- unlist(assignedMutations)
      genes$edgeAssignment <- assignedMutations - 1
      genes$driver <- genes$Gene.Name%in%drivers
      
      ##### check for convergent evolution ##### 
      duplicatedGenes <- genes$Gene.Name[duplicated(genes$Gene.Name)]
      genes$geneAndEdge <- paste0(genes$Gene.Name,'_',genes$edgeAssignment)
      convergentlyMutatedGenes <- sapply(duplicatedGenes, function(x) {
        if (length(unique(genes$geneAndEdge[genes$Gene.Name==x]))>1) {
          x  
        } else {
          NA
        }})
      convergentlyMutatedGenes <- convergentlyMutatedGenes[!is.na(convergentlyMutatedGenes)]
      genes$convergent <- genes$Gene.Name%in%convergentlyMutatedGenes
      
      ## get order of edges from top to bottom of tree (used for plotting and annotation)
      o <- orderEdges(tree = collapsed$tree)
      ## add index for ambiguous mutations to order vector
      o <- c(o,length(o)+1)
      ## use this to assign new edge labels
      genes$edgeAssignmentNew <- o[genes$edgeAssignment]
      
      ## add possible edges to ambiguously assigned mutations
      genes$possibleEdges <- NA
      for (i in which(idx!=1)) {
        edgeLabelOld <- theta[[i]]-1
        genes$possibleEdges[i] <- paste(sort(o[edgeLabelOld]),collapse=" ")
      }
      
      ## merge gene annotation and neoAntigen annotation
      genes <- as.data.table(genes)
      setkey(genes,'Chromosome','Position','Gene.Name')
      genes <- merge(genes,neoAntigens,all.x=TRUE)
      
      ## for mutations that are not assigned unambiguously add possible edge indices to name
      genes[!is.na(possibleEdges)]$Gene.Name <- paste(genes[!is.na(possibleEdges),Gene.Name],genes[!is.na(possibleEdges),possibleEdges])
      
      ##### create edge labels ##### 
      genes$type <- sapply(1:nrow(genes),function(i) {
        if (!is.na(genes[i,neoAntigen]) & genes[i,driver]) {
          1 ## neo antigen causing and driver mutation
        } else if (!is.na(genes[i,neoAntigen])) {
          2 ## neo antigen causing
        } else if (genes[i,driver]) {
          3 ## driver
        } else {
          4 ## neither neo antigen causing nor driver
        }
      })
      
      setkey(genes,edgeAssignment,type)
      assignedMutations <- lapply(1:(ecount(collapsed$tree)+1), function(x) genes[edgeAssignment==x])
      
      edgeLabelsNeo_Driver <- lapply(1:length(assignedMutations),function(li) if (nrow(assignedMutations[[li]])>0) {sapply(1:nrow(assignedMutations[[li]]), function(mi) {
        if (assignedMutations[[li]][mi,type]==1) {
          assignedMutations[[li]][mi,Gene.Name]
        } else { "" }})} else {""})
      edgeLabelsNeo_NonDriver <- lapply(1:length(assignedMutations),function(li) if (nrow(assignedMutations[[li]])>0) {sapply(1:nrow(assignedMutations[[li]]), function(mi) {
        if (assignedMutations[[li]][mi,type]==2) {
          assignedMutations[[li]][mi,Gene.Name]
        } else { "" }})} else {""})
      edgeLabelsAllOther_Driver <- lapply(1:length(assignedMutations),function(li) if (nrow(assignedMutations[[li]])>0) {sapply(1:nrow(assignedMutations[[li]]), function(mi) {
        if (assignedMutations[[li]][mi,type]==3) {
          assignedMutations[[li]][mi,Gene.Name]
        } else { "" }})} else {""})
      edgeLabelsAllOther_NonDriver <- lapply(1:length(assignedMutations),function(li) if (nrow(assignedMutations[[li]])>0) {sapply(1:nrow(assignedMutations[[li]]), function(mi) {
        if (assignedMutations[[li]][mi,type]==4) {
          assignedMutations[[li]][mi,Gene.Name]
        } else { "" }})} else {""})
      
      edgeLabelsNeo_Driver <- lapply(edgeLabelsNeo_Driver, function(x) paste(x,collapse='\n'))
      edgeLabelsNeo_NonDriver <- lapply(edgeLabelsNeo_NonDriver, function(x) paste(x,collapse='\n'))
      edgeLabelsAllOther_Driver <- lapply(edgeLabelsAllOther_Driver, function(x) paste(x,collapse='\n'))
      edgeLabelsAllOther_NonDriver <- lapply(edgeLabelsAllOther_NonDriver, function(x) paste(x,collapse='\n'))
      
      ##### creat edge labels for convergent mutations #####
      convGenes <- genes[convergent==TRUE]
      if (nrow(convGenes)>0) {
        convGenesTmp <- convGenes[,.(Chromosome,Position,Gene.Name,edgeAssignmentNew,driver)]
        fwrite(convGenesTmp,file=paste0(outName,"convergentGenes.txt"),quote = FALSE, sep="\t",row.names = FALSE)
      }
      convMutations <- lapply(1:(ecount(collapsed$tree)+1), function(x) convGenes[edgeAssignment==x])
      
      convNeo_Driver <- lapply(1:length(convMutations),function(li) if (nrow(convMutations[[li]])>0) {sapply(1:nrow(convMutations[[li]]), function(mi) {
        if (convMutations[[li]][mi,type]==1) {
          convMutations[[li]][mi,ID]
        } else { "" }})} else {""})
      convNeo_NonDriver <- lapply(1:length(convMutations),function(li) if (nrow(convMutations[[li]])>0) {sapply(1:nrow(convMutations[[li]]), function(mi) {
        if (convMutations[[li]][mi,type]==2) {
          convMutations[[li]][mi,ID]
        } else { "" }})} else {""})
      convAllOther_Driver <- lapply(1:length(convMutations),function(li) if (nrow(convMutations[[li]])>0) {sapply(1:nrow(convMutations[[li]]), function(mi) {
        if (convMutations[[li]][mi,type]==3) {
          convMutations[[li]][mi,ID]
        } else { "" }})} else {""})
      convAllOther_NonDriver <- lapply(1:length(convMutations),function(li) if (nrow(convMutations[[li]])>0) {sapply(1:nrow(convMutations[[li]]), function(mi) {
        if (convMutations[[li]][mi,type]==4) {
          convMutations[[li]][mi,ID]
        } else { "" }})} else {""})
      
      convNeo_Driver <- lapply(convNeo_Driver, function(x) paste(x,collapse='\n'))
      convNeo_NonDriver <- lapply(convNeo_NonDriver, function(x) paste(x,collapse='\n'))
      convAllOther_Driver <- lapply(convAllOther_Driver, function(x) paste(x,collapse='\n'))
      convAllOther_NonDriver <- lapply(convAllOther_NonDriver, function(x) paste(x,collapse='\n'))
      
      ##### plot tree only ##### 
      if (case == '290exOvary' & type=='WES') {
        width=8
        height=8
      } else if (case==308 & type=='WES') {
        width=10
        height=8
      } else if (case==315 & type=='WES') {
        width=8
        height=5
      } else if (type=='WES') {
        width=5
        height=5
      } else if (type=='targeted' & paste0(case,'_',sampleSet) %in% c('288_Fluids','288_WESsamples','290_Fluids','290exOvary_Fluids',
                                                                      '291_WESsamples','298main_Biopsies','298main_Fluids',
                                                                      '298main_WESsamples','308_Fluids','315_Fluids',
                                                                      '323_Fluids','323_WESsamples','328_All','328_Biopsies',
                                                                      '328_Fluids','328_WESsamples','330_Fluids')) {
        ## small
        width=5.5
        height=5.5
      } else if (type=='targeted' & paste0(case,'_',sampleSet) %in% c(
                                                                      '315_Biopsies')) {
        ## wide
        width=10
        height=7
      } else if (type=='targeted' & paste0(case,'_',sampleSet) %in% c(
                                                                      '290exOvary_Biopsies',
                                                                      '291_Fluids')) {
        ## high
        width=5
        height=8
      } else if (type=='targeted' & paste0(case,'_',sampleSet) %in% c('290_Biopsies','290_WESsamples','290exOvary_All',
                                                                      '290exOvary_WESsamples','298_Fluids',
                                                                      '298_WESsamples','298main_All','315_WESsamples','323_All','323_Biopsies','330_WESsamples')) {
        ## wide and high
        width=8
        height=8
      } else {
        c('288_All','288_Biopsies','290_All','291_All','291_Biopsies','298_All','298_Biopsies','308_All','308_Biopsies','308_WESsamples','315_All','330_All',
          '330_Biopsies')
        ## large
        ## wide and high
        width=10
        height=10
      }
      
      ## avoid overlaps in plots
      if (type=='targeted'&case=='290exOvary'&sampleSet=='All') {
        h.adjust <- rbind(c(4,-0.75),
                          c(5,-0.3))
      } else if (type=='targeted'&case==308&sampleSet=='All') { 
        h.adjust <- rbind(c(2,-0.75),
                          c(3,-1.5),
                          c(4,-1.5),
                          c(5,-1.75),
                          c(6,-1.25),
                          c(7,-1.25),
                          c(8,-1.25))
      } else if (type=='targeted'&case==308&sampleSet%in%c('Biopsies','WESsamples')) {
        h.adjust <- rbind(c(4,-0.5),
                          c(5,-0.25),
                          c(7,-0.9),
                          c(8,-0.9),
                          c(9,-0.9))
      } else if (type=='targeted'&case==315&sampleSet=='All') { 
        h.adjust <- rbind(c(6,1.25),
                          c(7,-0.5))
      }else if (type=='targeted'&case==315&sampleSet%in%c('Biopsies','WESsamples')) {
        h.adjust <- rbind(c(6,0.5),
                          c(7,0.5))
      } else if (type=='targeted'&case==298&sampleSet%in%c('All')) {
        h.adjust <- rbind(c(3,0.5))
      } else {
        h.adjust <- NULL
      }
      
      pdf(paste0(outName,'treeOnly.pdf'),width = width,height = height)
      par(mar=c(1,4,0,3),xpd=TRUE)
      plotTree(collapsed$tree,
               clones=collapsed$clones,
               e.length = collapsed$e.length,
               edge.arrow.mode='-',
               axis=TRUE,
               v.label.cex = 1,
               edge.width=3,
               h.adjust = h.adjust)
      title(ylab="Accumulated mutations")
      dev.off()
      
      
      tikz(paste0(outName,'treeOnly.tex'),width = width,height = height)
      par(mar=c(1,4,0,3),xpd=TRUE)
      plotTree(collapsed$tree,
               clones=collapsed$clones,
               e.length = collapsed$e.length,
               edge.arrow.mode='-',
               axis=TRUE,
               v.label.cex = 1,
               edge.width=3,
               h.adjust = h.adjust)
      title(ylab="Accumulated mutations")
      dev.off()
      
      ##### plot tree ##### 
      
      pdf(paste0(outName,'tree.pdf'),width = 9+2*length(o),height = max(9,11*max(collapsed$e.length)/60))
      par(mfrow=c(1,2),mar=c(3,4,3,0))
      layout(matrix(c(1,2),ncol=2),widths = c(1,2),heights = c(1,1))
      plotTree(collapsed$tree,
               clones=collapsed$clones,
               e.length = collapsed$e.length,
               edge.arrow.mode='-',
               axis=TRUE,
               v.label.cex = 1,
               edge.label=o,
               edge.label.color='black',
               h.adjust=h.adjust)
      title(main=paste("Case",case),ylab="Accumulated mutations")
      plot.new()
      n <- ecount(collapsed$tree)+1
      
      text(x=0,y=1.04,labels="All mutations",cex=0.75,adj=c(0,1),font=2)
      for (i in 1:n) {
        if (i==n) {
          text(x=(i-1)/n,y=1.02,labels="ambiguous",cex=0.75,adj=c(0,1),font=2)
        } else {
          text(x=(i-1)/n,y=1.02,labels=i,cex=0.75,adj=c(0,1),font=2)
        }
        text(x=(i-1)/n,y=1,labels=edgeLabelsNeo_Driver[[which(o==i)]],cex=0.65,adj=c(0,1),col='deeppink',font = 2)
        text(x=(i-1)/n,y=1,labels=edgeLabelsNeo_NonDriver[[which(o==i)]],cex=0.65,adj=c(0,1),col='deeppink')
        text(x=(i-1)/n,y=1,labels=edgeLabelsAllOther_Driver[[which(o==i)]],cex=0.65,adj=c(0,1),font = 2)
        text(x=(i-1)/n,y=1,labels=edgeLabelsAllOther_NonDriver[[which(o==i)]],cex=0.65,adj=c(0,1))
      }
      
      text(x=0,y=0.1,labels="Genes with multiple hits on different edges",cex=0.75,adj=c(0,1),font=2)
      for (i in 1:n) {
        if (i==n) {
          text(x=(i-1)/n,y=0.08,labels="ambiguous",cex=0.75,adj=c(0,1),font=2)
        } else {
          text(x=(i-1)/n,y=0.08,labels=i,cex=0.75,adj=c(0,1),font=2)
        }
        text(x=(i-1)/n,y=0.06,labels=convNeo_Driver[[which(o==i)]],cex=0.65,adj=c(0,1),col='deeppink',font = 2)
        text(x=(i-1)/n,y=0.06,labels=convNeo_NonDriver[[which(o==i)]],cex=0.65,adj=c(0,1),col='deeppink')
        text(x=(i-1)/n,y=0.06,labels=convAllOther_Driver[[which(o==i)]],cex=0.65,adj=c(0,1),font = 2)
        text(x=(i-1)/n,y=0.06,labels=convAllOther_NonDriver[[which(o==i)]],cex=0.65,adj=c(0,1))
      }
      
      dev.off()
      
      ##### plot mutation heatmap #####
      setkey(genes,order)
      ## order nodes along branches for nicer plotting
      oN <- orderNodes(collapsed$tree)
      oSamples <- unlist(collapsed$clones[oN])
      oSamples <- oSamples[!oSamples=="N"]
      nodeSizes <- sapply(1:length(oN),function(x) length(collapsed$clones[oN][[x]]))
      ## remove normal
      nodeSizes[1] <- nodeSizes[1]-1
      rowColors <- rep(c("white",colorsOrigin)[1:length(nodeSizes)],times=nodeSizes)
      ## invert order so that least mutated samples are at top of heatmap
      oSamples <- oSamples[length(oSamples):1]
      rowColors <- rowColors[length(rowColors):1]
      pdf(paste0(outName,'_heatmap.pdf'),width = max(10,nrow(posterior)/10+5),height = max(4,ncol(posterior)/3+1))
      plothm(data=posterior[,oSamples],
             rowGroup=genes$edgeAssignmentNew,
             colorsOrigin=colorsOrigin,
             rowColors=rowColors,
             containsAmbiguous=!all(is.na(genes$possibleEdges)),
             cexCol=1.15,
             cexRow=1.5)
      dev.off()
    }
  }
}
