args=(commandArgs(TRUE))

##args is now a list of character vectors
## First check to see if arguments are passed.
if(length(args)==0){
  stop("No arguments supplied.")
} else{
  metaDataFile = args[1]
  normalDirectory = args[2]
  Rdirectory = args[3]
  plotDirectory = args[4]
  captureRegionsFile = args[5]
  reference = args[6]
  vepCall = args[7]
}

library(superFreq)

data <- superFreq(metaDataFile=metaDataFile,
                  captureRegions=captureRegionsFile,
                  normalDirectory=normalDirectory,
                  Rdirectory=Rdirectory,
                  plotDirectory=plotDirectory,
                  reference=reference,
                  genome = "hg19",
                  BQoffset = 33,
                  cpus = 8,
                  outputToTerminalAsWell = T,
                  forceRedo = forceRedoNothing(),
                  normalCoverageDirectory = "",
                  systematicVariance = 0.02,
                  maxCov = 150,
                  cloneDistanceCut = -qnorm(0.01),
                  dbSNPdirectory = "superFreqDbSNP",
                  cosmicDirectory = "superFreqCOSMIC",
                  mode = "exome",
                  splitRun = F,
                  participants = "all",
                  manualStoryMerge = F,
                  vepCall = vepCall)
