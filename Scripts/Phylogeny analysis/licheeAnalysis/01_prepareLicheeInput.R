library(data.table)
pycloneData <- '~/Cancer Research UK, Cambridge Institute/Stephen-John Sammut - autopsy-project/stephen-analysis/pyclone/'
dataProcessed0 <- '~/Documents/Caldas/Data/dataProcessed/lichee/'
dataResults0 <- '~/Documents/Caldas/Data/dataResults/lichee/'

TYPES <- c("WES","targeted")

for (type in TYPES) {
  dataProcessed <- paste0(dataProcessed0,type)
  dataResults <- paste0(dataResults0,type,'/CPspectrum/')
  dir.create(paste0(dataProcessed,'input'),recursive = TRUE,showWarnings = FALSE)
  dir.create(dataResults,recursive = TRUE,showWarnings = FALSE)
  
  if (type=="WES") {
    CASES <- c(288,290,'298All',298,'298_en',308,315,323,328,330,'DET52')
  } else {
    CASES <- c(288,290,298,308,315,323,328,330,'DET52')
  }
  
  for (case in CASES) {
    if (case=='298_en'|case=='298All') {
      baseCase='298'
    } else {
      baseCase=case
    }
    
    ## read data
    if (type=="WES") {
      loci <- fread(paste0(pycloneData,'/exome/data/',case,'/tables/loci.tsv'))
    } else {
      ## if targeted analysis
      loci <- fread(paste0(pycloneData,'/ts/tables/',case,'/tables/loci.tsv'))
      ## read z-score for binarization
      zScore <- fread(paste0(pycloneData,'../murtaza-stephen-bind/Multi-regional_20170627/seq_results/merge.freq/merge.calls.sum.500-exomeSelect-normalAsControl/',case,'.500.sum.variants.zscore.txt'))
    }
    nSamples <- length(unique(loci$sample_id))
    clusterIdx <- loci$cluster_id[seq(1,nrow(loci),nSamples)]
    mutationIdx <- lapply(1:max(clusterIdx),function(x) which(clusterIdx==x))
    mutationIdx <- sapply(mutationIdx,function(x) paste0(x,collapse=','))
    ## create matrices with mutations in rows and samples in columns where values are vafs/cellular prevalences
    vaf <- dcast(loci,mutation_id ~ sample_id, value.var='variant_allele_frequency')
    vaf$cluster <- clusterIdx
    setkey(vaf,'cluster')
    loci <- dcast(loci,mutation_id ~ sample_id,value.var = 'cellular_prevalence')
    ## add info for each mutation
    loci[,c("#chr","position","description") := tstrsplit(loci$mutation_id,split=":|_")]
    ## delete mutation_id column
    loci[,mutation_id:=NULL]
    ## add normal sample
    loci[,"N" := rep(0,nrow(loci))]
    ## reorder columns as required by LICHeE
    n <- ncol(loci)
    setcolorder(loci, c((n-3):n,1:(n-4)))
    
    if (type=='targeted') {
      ## assign cluster ID to each mutation z Score entry and remove mutation if not in vaf data set
      clusterAssigment <- data.table(V1=gsub("(.*)\\_.*","\\1",vaf$mutation_id),
                                     cluster=vaf$cluster,
                                     order=1:nrow(vaf))
      zScore <- merge(x = clusterAssigment,y = zScore,by='V1' ,all.x = TRUE,all.y = FALSE)
      setkey(zScore,order)
      zScore$mutation_id <- vaf$mutation_id  
      
      if (case==298) {
        ## remove brain sample (different tumour)
        zScore[,`298-010`:=NULL]
      }
      
      ## select only columns that are present in both vaf and zScore
      vaf <- vaf[,colnames(vaf)%in%colnames(zScore),with=FALSE]
      loci <- loci[,c(rep(TRUE,4),colnames(loci)[-(1:4)]%in%colnames(vaf)),with=FALSE]
      zScore <- zScore[,match(colnames(vaf)[-1],colnames(zScore)),with=FALSE]
      ## remove na values
      zScoreTmp <- as.data.frame(zScore)
      zScoreTmp[is.na(zScoreTmp)] <- 0
      zScore <- data.table(zScoreTmp)
    }
  
    ## change sample names
    colnames(loci) <- gsub(paste0(baseCase,"-"),"",colnames(loci))
    
    ## read data
    if (type=="WES") {
      cluster <- fread(paste0(pycloneData,'/exome/data/',case,'/tables/cluster.tsv'))
    } else {
      cluster <- fread(paste0(pycloneData,'/ts/tables/',case,'/tables/cluster.tsv'))
    }
    clusterID <- 1:max(cluster$cluster_id)
    
    clusterSize <- cluster$size[seq(1,nrow(cluster),nSamples)]
    if (type=="WES") {
      ## remove clusters with fewer than x mutations
      x=3
      clusterID <- clusterID[clusterSize>=x]
      mutationIdx <- mutationIdx[clusterSize>=x]
      clusterSize <- clusterSize[clusterSize>=x]
      cluster <- cluster[size>=x]
      if (case=='298All') {
        ## remove ubiquitous noise cluster (containing 15 mutations)
        ## noise cluster filter does not work here because of independent tumours
        x=15
        clusterID <- clusterID[clusterSize!=x]
        mutationIdx <- mutationIdx[clusterSize!=x]
        clusterSize <- clusterSize[clusterSize!=x]
        cluster <- cluster[size!=x]
      }
    }
    
    ## create matrix with mutations in rows and samples in columns where values are cellular prevalences
    cluster <- dcast(cluster,cluster_id ~ sample_id,value.var = 'mean')
    if (type=='targeted') {
      ## only select samples for which zScore values were calculated 
      cluster <- cluster[,c(TRUE,colnames(cluster)[-1] %in% colnames(zScore)),with=FALSE]
    }
    
    ## plot cellular prevalence spectrum
    pdf(paste0(dataResults,case,'_CPspectrum.pdf'),width = 10,height = 6)
    par(mar=c(10,4,1,1),xpd=FALSE)
    plot(1:(ncol(cluster)-1),cluster[1,-1],type='l',col=rainbow(nrow(cluster))[1],ylim=c(0,1),xaxt='n',xlab='',ylab='cellular prevalence',lwd=2)
    grid()
    axis(1, at=1:(ncol(cluster)-1), labels=colnames(cluster)[-1],cex=0.1,las=3)
    for (i in 2:(nrow(cluster))) {
      lines(1:(ncol(cluster)-1),
            cluster[i,-1],
            col=rainbow(nrow(cluster))[i],
            lwd=2)
    }
    par(xpd=TRUE)
    legend(x=(ncol(cluster)-1)/2,y=-0.4,horiz=TRUE,xjust = 0.5,x.intersp = 0.25,
           legend=clusterSize,
           col=rainbow(length(clusterSize)),
           lty=1,lwd=2,
           bty='n',title = "Cluster size",title.adj = 0)
    dev.off()
    
    # calculate binary profile for each cluster
    vaf[,mutation_id:=NULL]
    if (type=="WES") {
      ## vaf greater than 0.01 in at least 40% of samples
      bin <- t(sapply(clusterID,function(x) (apply(vaf[cluster==x]>0.01,2,mean)>=0.4)*1))
    } else {
      ## zscore greater than 3 in at least 40% of samples
      bin <- t(sapply(clusterID,function(x) (apply(zScore[cluster==x]>3,2,mean)>=0.4)*1))
    }   
    ## delete cluster column
    bin <- bin[,-ncol(bin),drop=FALSE]
    
    if (type=="WES") {
      ## remove clusters that don't exceed 0.1 in any of the samples (likely to be noise)
      idx <- which(apply(cluster[,-1],1,max)<0.1)
    } else {
      ## if targeted analysis, remove clusters that have been classified as absent in all samples
      idx <- which(apply(bin,1,sum)==0)
    }
    ## remove clusters that are present in all samples but not at a high frequency (only keep cluster with highest average prevalence)
    removePresentInAll <- !apply(bin,1,function(x) any(x==0))
    meanCP <- apply(cluster[,-1],1,mean)
    if (any(removePresentInAll)) {
      removePresentInAll[removePresentInAll][which(meanCP[removePresentInAll]==max(meanCP[removePresentInAll]))] <- FALSE
    }
    
    if (case==288) {
      idx <- c(idx,4)
    }
    idx <- unique(sort(c(idx,which(removePresentInAll))))
    print(paste("Removing cluster(s):",idx))
    
    if (length(idx)>0) {
      clusterID <- clusterID[-idx]
      bin <- bin[-idx,,drop=FALSE]
      cluster <- cluster[-idx,,drop=FALSE]
      mutationIdx <- mutationIdx[-idx]
      clusterSize <- clusterSize[-idx]
    }
    
    ## select mutations that are in one of the selected clusters
    loci[clusterIdx%in%clusterID,]
    
    ## write data to file
    fwrite(loci,file=paste0(dataProcessed,'input/',case,'_cp.txt'),sep = '\t')
    
    
    ## plot cellular prevalence spectrum
    pdf(paste0(dataResults,case,'_CPspectrum_filtered.pdf'),width = 12,height = 6)
    par(mar=c(10,4,1,1),xpd=FALSE)
    plot(1:(ncol(cluster)-1),cluster[1,-1],type='l',col=rainbow(nrow(cluster))[1],ylim=c(0,1),xaxt='n',xlab='',ylab='cellular prevalence',lwd=2)
    grid()
    axis(1, at=1:(ncol(cluster)-1), labels=colnames(cluster)[-1],cex=0.1,las=3)
    for (i in 2:(nrow(cluster))) {
      lines(1:(ncol(cluster)-1),
            cluster[i,-1],
            col=rainbow(nrow(cluster))[i],
            lwd=2)
    }
    par(xpd=TRUE)
    legend(x=(ncol(cluster)-1)/2,y=-0.4,horiz=TRUE,xjust = 0.5,x.intersp = 0.25,
           legend=clusterSize,
           col=rainbow(length(clusterSize)),
           lty=1,lwd=2,
           bty='n',title = "Cluster size",title.adj = 0)
    dev.off()
    
    clusterBin <- as.matrix(cluster[,-1])
    clusterBin[bin==0] <- 0
    clusterBin <- as.data.table(data.frame(cluster_id=cluster$cluster_id,clusterBin))
    pdf(paste0(dataResults,case,'_CPspectrum_filtered_binarized.pdf'),width = 12,height = 6)
    par(mar=c(10,4,1,1),xpd=FALSE)
    plot(1:(ncol(cluster)-1),clusterBin[1,-1],type='l',col=rainbow(nrow(cluster))[1],ylim=c(0,1),xaxt='n',xlab='',ylab='cellular prevalence',lwd=2)
    grid()
    axis(1, at=1:(ncol(cluster)-1), labels=colnames(cluster)[-1],cex=0.1,las=3)
    for (i in 2:(nrow(cluster))) {
      lines(1:(ncol(cluster)-1),
            clusterBin[i,-1],
            col=rainbow(nrow(cluster))[i],
            lwd=2)
    }
    par(xpd=TRUE)
    legend(x=(ncol(cluster)-1)/2,y=-0.4,horiz=TRUE,xjust = 0.5,x.intersp = 0.25,
           legend=clusterSize,
           col=rainbow(length(clusterSize)),
           lty=1,lwd=2,
           bty='n',title = "Cluster size",title.adj = 0)
    dev.off()
    
    profile <- apply(bin,1,function(x) paste0(c(0,x),collapse =""))
    N <- rep(0,nrow(cluster))
    datFinal <- cbind(profile,N,cluster[,-1],mutationIdx)
    fwrite(datFinal,file=paste0(dataProcessed,'input/',case,'_clusters.txt'),sep = '\t',col.names = FALSE)
    
    
  }
}