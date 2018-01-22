## Modified from ascatNGS

library(ASCAT) ## needs to be ASCAT 2.5 or higher

##First read in the arguments listed at the command line
args=(commandArgs(TRUE))

##args is now a list of character vectors
## First check to see if arguments are passed.
if(length(args)==0){
  stop("No arguments supplied.")
} else{
	tumour_sample = args[1]
	gender = args[2]
	chrCount = as.numeric(args[3])
  rdat_out = args[4]
  purity = as.numeric(args[5])
  ploidy = as.numeric(args[6])

}

nonGenderChrs = chrCount - 2

load(rdat_out)

ascat.output = ascat.runAscat(ascat.bc, gamma = 1, rho_manual = purity, psi_manual = ploidy)

## make output files
if(!is.null(ascat.output$nA)) {

  ## make output for Catherine's visualization tool
  gCN = matrix(nrow = dim(ascat.bc$Tumor_LogR)[1], ncol = 9)
  rownames(gCN) = rownames(ascat.bc$Tumor_LogR)
  colnames(gCN) = c("Chromosome","Position","Log R", "segmented LogR", "BAF", "segmented BAF", "Copy number", "Minor allele", "Raw copy number")
  gCN[,1]=as.vector(ascat.bc$SNPpos[,1])
  gCN[,2]=ascat.bc$SNPpos[,2]
  # X chr is the first one after the main, all includes X+Y
  gCN[gCN[,1]=="X",1]=chrCount-1
  gCN[gCN[,1]=="Y",1]=chrCount
  gCN[,3]=ascat.bc$Tumor_LogR[,1]
  gCN[,4]=ascat.bc$Tumor_LogR_segmented[,1]
  gCN[,5]=ascat.bc$Tumor_BAF[,1]
  gCN[rownames(ascat.bc$Tumor_BAF_segmented[[1]]),6]=ascat.bc$Tumor_BAF_segmented[[1]][,1]
  seg = ascat.output$segments

  majorallele = NULL
  minorallele = NULL
  for (i in 1:dim(seg)[1]) {
    start = which(ascat.bc$SNPpos[,1]==seg[i,2] & ascat.bc$SNPpos[,2]==seg[i,3])[1]
    end = which(ascat.bc$SNPpos[,1]==seg[i,2] & ascat.bc$SNPpos[,2]==seg[i,4])[1]
    majorallele=c(majorallele,rep(seg[i,5],end-start+1))
    minorallele=c(minorallele,rep(seg[i,6],end-start+1))
  }

  gCN[,7]=majorallele+minorallele
  gCN[,8]=minorallele

  rho = ascat.output$aberrantcellfraction[1]
  psi = ascat.output$psi[1]

  gCN[,9]=(2*rho - 2 + (2*(1-rho)+rho*psi)*2^(ascat.bc$Tumor_LogR[,1]/0.55))/rho

  write.table(gCN,paste(tumour_sample, ".copynumber.txt",sep=""),sep="\t",row.names=T,col.names=NA,quote=F)

  ## make input for Caveman:
  cavemanSegs = cbind(seg[,2],
          seg[,3],
          seg[,4],2,1,
          seg[,5]+seg[,6],
          seg[,6])

  # make Major germline 1 when MALE
  if(gender=="XY") {
    cavemanSegs[cavemanSegs[,1] %in% c('X','Y'),4] = 1
    cavemanSegs[cavemanSegs[,1] %in% c('X','Y'),5] = 0
  }

  rownames(cavemanSegs) = 1:dim(cavemanSegs)[1]

  write.table(cavemanSegs,paste(tumour_sample,".copynumber.caveman.csv",sep=""),row.names=T,col.names=F,sep=",",quote=F)

  normalContamination = 2*(1-rho)/(2*(1-rho)+rho*ascat.output$ploidy[1])

  ss = matrix(ncol=1,nrow=5)
  rownames(ss) = c("NormalContamination","Ploidy","rho","psi", "goodnessOfFit")
  ss[,1] = c(normalContamination,ascat.output$ploidy,rho,psi,ascat.output$goodnessOfFit)

  write.table(ss,paste(tumour_sample,".samplestatistics.txt",sep=""),row.names=T,col.names=F,quote=F)

}

## remove purity, ploidy and diploidchromosomes from RData file, otherwise we can't reiterate
rm(purity)
rm(ploidy)

## save file
save.image(rdat_out)

q(save="no")


