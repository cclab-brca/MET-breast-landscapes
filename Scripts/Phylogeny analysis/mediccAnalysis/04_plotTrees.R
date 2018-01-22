library(phangorn)
library(igraph)
library(oncoNEM)
library( 'tikzDevice' )

source('~/Documents/Caldas/Data/breastMetRepo/treeomicsAnalysis/04_plotting.R')
source('~/Documents/Caldas/Data/breastMetRepo/treeomicsAnalysis/05_collapseTree.R')
resdir <- "~/Documents/Caldas/Data/dataResults/medicc_final/"
dir.create(paste0(resdir,'plots'),showWarnings = FALSE)

CASES <- c('288','290exOvary','298main','308','315','323','328','330','DET52')

for (case in CASES) {
  outName <- paste0(resdir,'plots/',case)
  
  if (case=='290exOvary') {
    baseCase=290
  } else if (case=='298main'|case=='298mainSingleD') {
    baseCase=298
  } else if (case=='308SingleFiltered500k') {
    baseCase=308
  } else {
    baseCase=case
  }

  noFINAL <- FALSE
  if (file.exists(paste0(resdir,case,'/tree_final.new'))) {
    M <- read.tree(paste0(resdir,case,'/tree_final.new'))
  } else {
    M <- read.tree(paste0(resdir,case,'/tree_nj.new'))
    M <- root(M,'diploid')
  }
  
  M$tip.label <- gsub(pattern = paste0(baseCase,'-'),replacement = "",x = M$tip.label)
  pdf(paste0(outName,'_medicc_phylogram.pdf'),width = 7,height = 7)
  plot(M,"phylogram",main="")
  dev.off()
  
  g <- as.igraph(M)
  el <- get.edgelist(g)
  idx <- which(el[,2]=='diploid')
  el[idx,] <- el[idx,2:1]
  internalVertices <- unique(el[,1])
  internalVerticesOrdered <- internalVertices[order(shortest.paths(g,v='diploid',to=internalVertices))]
  idx <- unlist(sapply(internalVerticesOrdered, function(x) which(el[,1]==x)))
  el <- el[idx,]
  edge.length <- M$edge.length[idx]
  g <- graph.edgelist(el,directed = TRUE)

  ## replace each node by number (needed for plotting)
  Nodes <- unique(as.vector(el))
  labs <- Nodes
  labs[grep(pattern = 'internal',x=labs)] = ''
  el <- matrix(as.numeric(factor(as.vector(el),levels=Nodes)),ncol=2,byrow = FALSE)
  # el <- matrix(as.numeric(factor(el,levels=Nodes)),ncol=2,byrow = TRUE)
  # nNodes <- length(Nodes)
  # for (i in 1:nrow(el)) {
  #   for (j in 1:ncol(el)) {
  #     el[i,j] <- (1:nNodes)[which(Nodes==el[i,j])]
  #   }
  # }
  g <- graph.edgelist(el)
  
  if (any(M$edge.length==0)) {
    cTree <- collapseTree(tree = g,e.length = edge.length,labels = labs)
  } else {
    cTree <- list(tree=g,
                  e.length=edge.length,
                  labels=labs)
  }
  
  ## plot
  if (case%in%c('290exOvary',308)) {
    pdf(paste0(outName,'_medicc_tree.pdf'),width = 11,height = 7)
  } else {
    pdf(paste0(outName,'_medicc_tree.pdf'),width = 7,height = 7)
  }
  par(mar=c(0.1,4,0.1,0.1),xpd=TRUE)
  plotTree(cTree$tree,
           clones=cTree$labels,
           e.length = cTree$e.length,
           edge.arrow.mode='-',
           axis=TRUE,
           v.label.cex = 1)
  title(ylab="Evolutionary distance")
  dev.off()
  
  if (case%in%c('290exOvary',308)) {
    tikz(paste0(outName,'_medicc_tree.tex'),width = 11,height = 7)
  } else {
    tikz(paste0(outName,'_medicc_tree.tex'),width = 7,height = 7)
  }
  par(mar=c(0.1,4,0.1,0.1),xpd=TRUE)
  plotTree(cTree$tree,
           clones=cTree$labels,
           e.length = cTree$e.length,
           edge.arrow.mode='-',
           axis=TRUE,
           v.label.cex = 1)
  title(ylab="Evolutionary distance")
  dev.off()
  
  save(cTree,file=paste0(outName,'_medicc_tree.RData'))
}
