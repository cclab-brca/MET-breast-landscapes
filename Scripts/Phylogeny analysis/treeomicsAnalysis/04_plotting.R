plotTree <- function(tree,clones=NULL,e.length=NULL,scale.width=NULL,axis=FALSE,
                     label.length="max",f.x=0.1,f.y=0.1,v.label.cex=1,h.adjust=NULL,layout.help=FALSE,...) {
  if (is.vector(tree)) {
    g <- igraph::graph.edgelist(el=cbind(tree+1,(1:length(tree))+1))
  } else if (igraph::is.igraph(tree)) {
    g <- tree
  } else {
    stop("tree has unknown input format. Should be in vector format or an igraph
         object")
  }
  if (!is.null(clones)) {
    if(length(clones)!=igraph::vcount(g)) {
      stop(paste0("List of clones does not match tree.
                  Has to have length vcount(g)=",igraph::vcount(g),'.'))
    }
    if (label.length=="max") {
      igraph::V(g)$label <- sapply(clones,function(x) toString({
        if ("N"%in%x) {
          c("N",sort(x[x!="N"]))
        } else {
          sort(x)
        }}))
    } else {
      igraph::V(g)$label <- sapply(clones,function(x) {
        x <- sort(x)
        if ("N" %in% x) {
         x <- c("N",x[x!="N"]) 
        }
        y <- vector()
        while (length(x)>label.length) {
          y <- c(y,x[1:label.length],'\n')
          x <- x[-(1:label.length)]
        }
        if (length(x)>0) {
          y <- c(y,x)
        }
        y <- toString(y)
        gsub(' \n,',replacement = '\n',x = y)
      } 
      )
    }
  } else {
    igraph::V(g)$label <- c(1:(igraph::vcount(g)))
  }
  
  g.orig <- g
  
  nv <- igraph::vcount(g)
  layout <- igraph::layout.reingold.tilford(g)
  tree.depth <- max(layout[,2])
  leaf.nodes <- (1:nv)[igraph::degree(g,mode="out")==0]
  leaf.depth <- igraph::shortest.paths(g,v = 1,to = leaf.nodes)
  ## remove leaves with maximum depth
  leaf.nodes <- leaf.nodes[leaf.depth<tree.depth]
  leaf.depth <- leaf.depth[leaf.depth<tree.depth]
  counter <- nv
  if (length(leaf.nodes)>=1) {
    for (i in 1:length(leaf.nodes)) {
      ## add nodes and edges so that every leaf node has maximum depth
      delta <- tree.depth-leaf.depth[i]
      g <- g+delta
      #     g <- g+path(c(leaf.nodes[i],counter+1:delta))
      if (delta==1) {
        g[from=c(leaf.nodes[i]),
          to=counter+1] <- TRUE
      } else {
        g[from=c(leaf.nodes[i],counter+(1:(delta-1))),
          to=counter+(1:delta)] <- TRUE
      }
      counter <- counter+delta
    }
  }
  ## get new layout matrix.
  layout.new <- igraph::layout.reingold.tilford(g)[1:nv,]
  if (!is.null(e.length)) {
    el <- igraph::get.edgelist(g.orig)  
    el <- matrix(as.numeric(el),ncol=2)
    Pred <- transitiveClosure(el, returnAdjMat = TRUE)
    for (i in 1:nrow(el)) {
      ## update y coordinates
      layout.new[Pred[el[i,2],]==1,2] <- layout.new[Pred[el[i,2],]==1,2]-e.length[i]+1
    }
  }
  if (!is.null(scale.width)) {
    layout.new[,1] <- layout.new[,1]*scale.width
  }
  ## set y value of normal to 0
  layout.new[,2] <- layout.new[,2]-max(layout.new[,2])
  if (layout.help) {
    print("Original Layout matrix")
    print(layout.new)
  }
  if (!is.null(h.adjust)) {
    h.adjust.vector <- rep(0,nrow(layout.new))
    h.adjust.vector[h.adjust[,1]] <- h.adjust[,2]
    layout.new[,1] <- layout.new[,1]+h.adjust.vector
  }
  if (is.null(clones)) {
    plot(0, type="n", ann=FALSE, axes=FALSE, xlim=extendrange(layout.new[,1]), 
         ylim=extendrange(layout.new[,2]))
    if (axis) {
      x <- axis(2,labels = FALSE)
      axis(2,labels = as.character(-x),at=x)
    }
    igraph::plot.igraph(g.orig, layout=layout.new, rescale=FALSE, add=TRUE,...)
  } else {
    plot(0, type="n", ann=FALSE, axes=FALSE, xlim=extendrange(layout.new[,1],f=f.x), 
         ylim=extendrange(layout.new[,2],f=f.y))
    if (axis) {
      x <- axis(2,labels = FALSE)
      axis(2,labels = as.character(-x),at=x)
    }
    
    vshape <- rep('rectangle',length(igraph::V(g.orig)$label))
    vshape[strwidth(igraph::V(g.orig)$label)==0] <- 'none'
    
    igraph::plot.igraph(g.orig, layout=layout.new, rescale=FALSE, add=TRUE,
                        vertex.shape=vshape,
                        vertex.size=(strwidth(igraph::V(g.orig)$label) + strwidth("o")) * 90 * v.label.cex,
                        vertex.size2=(strheight(igraph::V(g.orig)$label) + strheight("o")) * 100 * v.label.cex,
                        vertex.label.cex=v.label.cex,
                        vertex.color="white",
                        ...)
  }  
}
