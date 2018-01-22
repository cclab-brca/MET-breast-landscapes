collapseTree <- function(tree,e.length,labels) {
  ## assign IDs to edges
  E(tree)$idx <- 1:ecount(tree)
  
  ## find indices of nodes we want to join
  DistWeighted <- shortest.paths(tree,weights = e.length,mode = 'out')
  DistSteps <- shortest.paths(tree,mode='out')
  collapseIdx <- which(DistWeighted==0&DistSteps==1,arr.ind = TRUE)
  
  collapsedNodes <- as.list(1:vcount(tree))
  for (i in nrow(collapseIdx):1) {
    ## add edges to children of the node that is about to be removed
    el <- get.edgelist(tree,names = FALSE)
    edgeIndices <- E(tree)$idx
    idxAdd <- which(el[,1] == collapseIdx[i,2])
    if (length(idxAdd)>0) {
      for (j in 1:length(idxAdd)) {
        tree <- add.edges(tree,c(collapseIdx[i,1],el[,2][idxAdd[j]]),attr=list(idx=edgeIndices[idxAdd[j]]))
        tree[from=el[idxAdd[j],1],to=el[idxAdd[j],2]] <- 0
      }
    }
    
    ## delete edge to node that is about to be removed
    tree[from=collapseIdx[i,1],to=collapseIdx[i,2]] <- 0
    collapsedNodes[[collapseIdx[i,1]]] <- c(collapsedNodes[[collapseIdx[i,1]]],
                                            collapsedNodes[[collapseIdx[i,2]]])
    collapsedNodes[[collapseIdx[i,2]]] <- numeric()
  }
  
  ## delete unconnected vertices
  tree <- delete.vertices(tree, which(degree(tree)==0) )
  
  el <- get.edgelist(tree,names = FALSE)
  treeNew <- graph.edgelist(el[order(el[,1]),,drop=FALSE])
  E(treeNew)$idx=E(tree)$idx[order(el[,1])]
  
  ## make new node Labels
  collapsedNodes <- collapsedNodes[sapply(collapsedNodes,function(x) length(x)!=0)]
  
  labelsNew <- sapply(collapsedNodes,function(x) {
    labs=labels[x]
    paste0(labs[labs!=""],collapse=", ")}
  )

  return(list(tree=treeNew,e.length=e.length[E(treeNew)$idx],labels=labelsNew,labelsSplit=collapsedNodes))
}