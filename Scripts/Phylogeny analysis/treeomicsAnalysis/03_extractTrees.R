
library(igraph)
library(oncoNEM)
library(data.table)


extractTree <- function(dat) {
  
  tree <- data.frame(from=rep(0,length(dat)),
                     to=rep(0,length(dat)),
                     length=rep(0,length(dat)),
                     AcqNames=rep(0,length(dat)),
                     AcqPos=rep(0,length(dat)),
                     AcqIdx=rep(0,length(dat)))
  
  ## initialize
  currIdx <- 1
  fromIdx <- 1
  tree$from[currIdx] <- dat[1]
  if (length(grep("\\+",dat[2]))>0) {
    tree$length[currIdx] <- dat[2]
    tree$AcqNames[currIdx] <- dat[3]
    tree$AcqPos[currIdx] <- dat[4]
    tree$AcqIdx[currIdx] <- dat[5]
    tree$to[currIdx] <- dat[6]
    i=6
  } else {
    tree$to[currIdx] <- dat[2]
    i=2
  }
  parentNodes <- dat[1]
  
  while (i<=length(dat)) {
    
    if (length(grep("\\[",dat[i]))>0) {
      #tree$to[currIdx] <- dat[i]
      currIdx=currIdx+1
      tree$from[currIdx] <- dat[i]
      parentNodes <- c(parentNodes,dat[i])
      if (length(grep("\\+",dat[i+1]))>0) {
        tree$length[currIdx] <- dat[i+1]
        tree$AcqNames[currIdx] <- dat[i+2]
        tree$AcqPos[currIdx] <- dat[i+3]
        tree$AcqIdx[currIdx] <- dat[i+4]
        tree$to[currIdx] <- dat[i+5]
        i=i+5
      } else {
        tree$length[currIdx] <- 0.1
        tree$to[currIdx] <- dat[i+1]
        i=i+1
      }
    } else if (length(grep("\\]",dat[i]))>0) {
      while (i<length(dat) && length(grep("\\]",dat[i+1]))>0) {
        i=i+1
        parentNodes <- parentNodes[-length(parentNodes)]
      }
      if (i==length(dat)) {
        break()
      }
      currIdx=currIdx+1
      parentNodes <- parentNodes[-length(parentNodes)]
      tree$from[currIdx] <- parentNodes[length(parentNodes)]
      if (length(grep("\\+",dat[i+1]))>0) {
        tree$length[currIdx] <- dat[i+1]
        tree$AcqNames[currIdx] <- dat[i+2]
        tree$AcqPos[currIdx] <- dat[i+3]
        tree$AcqIdx[currIdx] <- dat[i+4]
        tree$to[currIdx] <- dat[i+5]
        i=i+5
      } else {
        tree$length[currIdx] <- 0.1
        tree$to[currIdx] <- dat[i+1]
        i=i+1
      }
    } else if (length(grep("\\+",dat[i+1]))>0) {
      currIdx=currIdx+1
      tree$from[currIdx] <- parentNodes[length(parentNodes)]
      tree$length[currIdx] <- dat[i+1]
      tree$AcqNames[currIdx] <- dat[i+2]
      tree$AcqPos[currIdx] <- dat[i+3]
      tree$AcqIdx[currIdx] <- dat[i+4]
      tree$to[currIdx] <- dat[i+5]
      i=i+5
    } else if (length(grep("\\+",dat[i]))>0) {
      currIdx=currIdx+1
      tree$from[currIdx] <- parentNodes[length(parentNodes)]
      tree$length[currIdx] <- dat[i]
      tree$AcqNames[currIdx] <- dat[i+1]
      tree$AcqPos[currIdx] <- dat[i+2]
      tree$AcqIdx[currIdx] <- dat[i+3]
      tree$to[currIdx] <- dat[i+4]
      i=i+4
    } else if (length(grep("\\]",dat[i+1]))>0) {
      i=i+1
    } else {
      currIdx=currIdx+1
      tree$from[currIdx] <- parentNodes[length(parentNodes)]
      tree$length[currIdx] <- 0.1
      tree$to[currIdx] <- dat[i+1]
      i=i+1
    }
  }
  
  idx <- apply(tree,1,function(x) all(x=="0"))
  tree <- tree[!idx,]
  tree$from <- gsub(pattern = "\\[","",tree$from)
  tree$to <- gsub(pattern = "\\[","",tree$to)
  tree$length <- gsub(pattern = "\\+","",tree$length)
  
  ## make igraph object
  treeGraph <- graph.edgelist(el=as.matrix(tree[,1:2]))
  
  ## set labels
  labels <- gsub("SC\\~.*","",V(treeGraph)$name)
  labels <- gsub("X","",labels)
  labels[1] <- "N"
  V(treeGraph)$name <- 1:length(V(treeGraph))
  edgeLength <- as.numeric(tree$length)
  edgeLength[is.na(edgeLength)] <- 0
  
  tree$AcqNames <- gsub("%Acquiredmutations:","",tree$AcqNames)
  tree$AcqPos <- gsub("%Acquiredmutations:","",tree$AcqPos)
  tree$AcqIdx <- gsub("%Acquiredmutations:","",tree$AcqIdx)
  
  edgeAttributes <- data.table(matrix(0,nrow = sum(floor(as.numeric(tree$length))),ncol=4)) #6
  colnames(edgeAttributes)  <- c('EdgeIdx','Chromosome','Position','Mut')#,'GeneName','MutIdx')
  class(edgeAttributes$EdgeIdx) <- "numeric"
  class(edgeAttributes$Position) <- "numeric"
  #class(edgeAttributes$MutIdx) <- "numeric"
  edgeAttributes$EdgeIdx <- rep(1:nrow(tree),times=floor(as.numeric(tree$length)))
  trueEdges <- which(tree$length>=1)
  tmp <- matrix(strsplit(paste0(tree$AcqPos[trueEdges],collapse=","),split = ',|\\_\\_')[[1]],ncol=3,byrow = TRUE)
  edgeAttributes$Chromosome <- tmp[,1]
  edgeAttributes$Position <- as.numeric(tmp[,2])
  edgeAttributes$Mut <- tmp[,3]
  #edgeAttributes$GeneName <- strsplit(paste0(tree$AcqNames,collapse=","),',')[[1]]
  #edgeAttributes$MutIdx <- as.numeric(strsplit(paste0(tree$AcqIdx,collapse=","),',')[[1]])+1
  
  return(list(tree=treeGraph,labels=labels,edgeLength=edgeLength,edgeAttributes=edgeAttributes))
}  
