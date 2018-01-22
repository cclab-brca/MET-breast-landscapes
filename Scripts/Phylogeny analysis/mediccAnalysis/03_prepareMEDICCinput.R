library(data.table)

##First read in the arguments listed at the command line
args=(commandArgs(TRUE))

##args is now a list of character vectors
## First check to see if arguments are passed.
if(length(args)==0){
  stop("No arguments supplied.")
} else{
  samplesFile = args[1]
  ascatOutDir = args[2]
  mediccInDir = args[3]
}

################################################################################
#### create dirs and load sample information ####

dir.create(mediccInDir,recursive = TRUE,showWarnings = FALSE)

## read in sample information
samples <- read.table(samplesFile,stringsAsFactors = FALSE)$V1

################################################################################
#### 

## Read CN data for first sample
s <- samples[1]
dat <- fread(paste0(ascatOutDir,'/',s,'.copynumber.txt'),
             colClasses = c('character','character',rep('numeric',8)))
setnames(dat,c("V1","Chr","Pos","LogR","segLogR","BAF","segBAF","CN","MinAllele","RawCN"))
## calculate major allele
dat[,MajAllele:=(CN-MinAllele)]
## test if total copynumber is 0 anywhere
dat[,CN0:=(CN==0)]
## remove unneeded columns
dat[,c("V1","LogR","segLogR","BAF","segBAF","CN","RawCN"):=NULL]
setnames(dat,c("MinAllele","MajAllele"),c(paste0(s,"MinAllele"),paste0(s,"MajAllele")))


## Read CN data for other samples one by one
for (s in samples[-1]) {
  
  dat2 <- fread(paste0(ascatOutDir,'/',s,'.copynumber.txt'),
               colClasses = c('character','character',rep('numeric',8)))
  setnames(dat2,c("V1","Chr","Pos","LogR","segLogR","BAF","segBAF","CN","MinAllele","RawCN"))
  ## calculate major allele
  dat2[,MajAllele:=(CN-MinAllele)]
  ## test if total copynumber is 0 anywhere
  dat2[,CN0:=(CN==0)]
  ## combine with data from other samples
  dat[,c(paste0(s,"MinAllele"),paste0(s,"MajAllele")):=dat2[,.(MinAllele,MajAllele)]]
  dat[,CN0:= (CN0 | dat2[,CN0])]

}

## remove rows where total CN is 0 in any of the samples
dat <- dat[CN0==FALSE]
## delete CN0 column
dat[,CN0:=NULL]
## save SNP positions
pos <- dat[,Pos]
## delete position column
dat[,Pos:=NULL]

## duplicate data.table and shift by one
datShift <- dat[-nrow(dat),]
datShift <- rbind(as.list(rep(999,ncol(datShift))),datShift)

## find rowIndex where at least one copy number or chromosome changes
idx <- which(rowSums(dat!=datShift)>0)

dat <- dat[idx,]
segments <- data.frame(chr=dat[,Chr],
                       start=pos[idx],
                       end=pos[c(idx[-1]-1,length(pos))])
save(segments,file=file.path(mediccInDir,'segments.RData'))

## add maj and min CN for diploid
dat[,diploid:=1]

## write to fasta file for each chromosome and maj/minor allele separately
if (file.exists(file.path(mediccInDir,'desc.txt'))) {
  file.remove(file.path(mediccInDir,'desc.txt'))
}

replaceValuesGreater9 = function(DT) {
  for (i in names(DT))
    DT[get(i)>9,(i):=9]
}


for (chr in as.numeric(unique(dat$Chr))) {
  datChr <- dat[Chr==chr]
  datChr[,Chr:=NULL]
  
  ## replace values greater than 9 by 9 (otherwise length of string changes wich causes MEDICC to break down)
  replaceValuesGreater9(datChr)
  
  fileMaj <- file.path(mediccInDir,paste0("major_chr",chr,'.fasta'))
  fileMin <- file.path(mediccInDir,paste0("minor_chr",chr,'.fasta'))
  
  write(">diploid",file=fileMaj)
  write(">diploid",file=fileMin)
  
  write(datChr[,diploid],file=fileMaj,sep = "",append = TRUE,ncol=80)
  write(datChr[,diploid],file=fileMin,sep = "",append = TRUE,ncol=80)
  
  for (s in samples) {
    write(paste0(">",s),file=fileMaj,append = TRUE)
    write(paste0(">",s),file=fileMin,append = TRUE)
    
    write(datChr[,get(paste0(s,'MajAllele'))],file=fileMaj,sep = "",append = TRUE,ncol=80)
    write(datChr[,get(paste0(s,'MinAllele'))],file=fileMin,sep = "",append = TRUE,ncol=80)
  }
  
  write(paste(paste0('chrom',chr),paste0('major_chr',chr,'.fasta'),paste0('minor_chr',chr,'.fasta')),
          file=file.path(mediccInDir,'desc.txt'), sep ='\t',
          append=TRUE)
  
}



