library(ASCAT) ## TODO ASCATmod
# 
ascatOutDir="~/Documents/Caldas/Data/dataProcessed/ascat/ascatOutput/288/"
ascatInDir="~/Documents/Caldas/Data/dataProcessed/ascat/ascatInput/288/"
pid='298'
samplesFile="~/Documents/Caldas/Data/dataRaw/samples_288.txt"
normal_sample="288-024"
wsample=NULL
gender="XX"
chrCount=24

# # 
# ascatOutDir="~/Documents/Caldas/Data/dataProcessed/ascat/ascatOutput/291"
# ascatInDir="~/Documents/Caldas/Data/X_finalAnalysis/intdata/ascatInput/Patient_291"
# pid=291
# samplesFile="~/Documents/Caldas/Data/X_finalAnalysis/intdata/sampleFiles/samples_291.txt"
# normal_sample="291-005"
# wsample=NULL
# gender="XX"
# chrCount=24


##First read in the arguments listed at the command line
args=(commandArgs(TRUE))
##args is now a list of character vectors
## First check to see if arguments are passed.
if(length(args)==0){
  stop("No arguments supplied.")
} else{
  ascatInDir = args[1]
  ascatOutDir = args[2]
  pid = args[3]
  samplesFile = args[4]
  normal_sample = args[5]
  wsample = args[6]
  gender = args[7]
  chrCount = as.numeric(args[8])
}

# print(paste('ascatInDir',ascatInDir))
# print(paste('ascatOutDir',ascatOutDir))
# print(paste('pid',pid))
# print(paste('samplesFile',samplesFile))
# print(paste('normal_sample',normal_sample))
# print(paste('wsample',wsample))
# print(paste('gender',gender))
# print(paste('chrCount',chrCount))
# 
# stop()


nonGenderChrs = chrCount - 2

dir.create(ascatOutDir,recursive=TRUE,showWarnings = FALSE)
setwd(ascatOutDir)

## read in sample information
samples <- read.table(samplesFile,stringsAsFactors = FALSE)$V1
nsamples <- length(samples)

if (is.null(wsample)||wsample=="NULL"||is.na(wsample)) {
  wsample = NULL
} else {
  wsample = strsplitAsListOfIntegerVectors(wsample)[[1]]
  if (nsamples!=length(wsample)) {
    stop("Sample weights should be a vector of length nsamples.")
  }
}

################################################################################
#### allele specific multi sample segmentation ####

## load data
ascat.bcAll = ascat.loadData(Tumor_LogR_file = file.path(ascatInDir,'/merged.tumour.LogR.txt'),
                             Tumor_BAF_file = file.path(ascatInDir,'/merged.tumour.BAF.txt'),
                             Germline_LogR_file = file.path(ascatInDir,paste0(normal_sample,'.LogR.txt')),
                             Germline_BAF_file = file.path(ascatInDir,paste0(normal_sample,'.BAF.txt')),
                             gender = rep(gender,nsamples), 
                             chrs=c(1:nonGenderChrs,"X"), 
                             sexchromosomes = c("X")) 

## run multi sample segmentation
ascat.bcAll = ascat.asmultipcf(ASCATobj = ascat.bcAll,
                               ascat.gg = NULL,
                               wsample=wsample,
                               selectAlg = "exact",
                               refine = TRUE)

################################################################################
#### run ascat separately on each segmented sample ####

## select single sample from ascat output
for (sample in 1:length(ascat.bcAll$samples)) {
  tumour_sample=ascat.bcAll$samples[sample]
  ascat.bc = list(Tumor_LogR = ascat.bcAll$Tumor_LogR[,sample,drop=FALSE],
                  Tumor_BAF = ascat.bcAll$Tumor_BAF[,sample,drop=FALSE],
                  Tumor_LogR_segmented = ascat.bcAll$Tumor_LogR_segmented[,sample,drop=FALSE],
                  Tumor_BAF_segmented = ascat.bcAll$Tumor_BAF_segmented[sample],
                  Germline_LogR = ascat.bcAll$Germline_LogR,
                  Germline_BAF = ascat.bcAll$Germline_BAF,
                  SNPpos = ascat.bcAll$SNPpos,
                  ch = ascat.bcAll$ch,
                  chr = ascat.bcAll$chr,
                  chrs = ascat.bcAll$chrs,
                  samples = colnames(ascat.bcAll$Tumor_LogR)[sample], gender = ascat.bcAll$gender[sample],
                  sexchromosomes = ascat.bcAll$sexchromosomes, failedarrays = ascat.bcAll$failedarrays)
  
  ascat.output = ascat.runAscat(ASCATobj = ascat.bc, gamma = 1)
  
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
  
  ## save file
  save.image(file.path(ascatOutDir,paste0(tumour_sample,'.RData')))
  
}
