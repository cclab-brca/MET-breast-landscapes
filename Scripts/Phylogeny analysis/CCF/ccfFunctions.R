library(GenomicRanges)
library(Hmisc)
library(sequenza)
library(bootstrap)
library(boot)

###Function to compute CCF 
###Code based on McGranahan et al. 2015 (doi: 10.1126/scitranslmed.aaa1408 )

computeCCF <- function(vaf, depth, purity, tumCN){
	
	#VAF - confidence intervals of VAF estimate
	vafCI_lower <- round(binom.test(round(vaf*depth,0), depth)$conf[1],4)
	vafCI_upper <- round(binom.test(round(vaf*depth,0), depth)$conf[2],4)
	
	#Multiplicity - average number of mutant copies per cell.  Can be above 1
	#Assume no subclonal copy number events
	mult_point <- round((((binconf(x=vaf*depth, n=depth))/purity) * ( purity*tumCN + (1-purity)*2))[1],4)
	mult_lower <- round((((binconf(x=vafCI_lower*depth, n=depth))/purity) * ( purity*tumCN + (1-purity)*2))[1],4)
	mult_upper <- round((((binconf(x=vafCI_upper*depth, n=depth))/purity) * ( purity*tumCN + (1-purity)*2))[1],4)
	
	#Absolute CCF
	possibleCCFs <- seq(0.01, 1, 0.01)
	possibleVafs <- (purity*possibleCCFs)/((2*(1-purity)) + (purity*tumCN) ) #Expected VAF for each CCF
	probs <- dbinom(x=round(vaf*depth), size=depth, prob=possibleVafs) #Prob of observed VAF
	names(probs) <- possibleCCFs
	probs_norm <- probs/sum(probs) #Normalise to get posterior distribution
	probs_sort <- sort(probs_norm, decreasing=T)
	probs_cum <- cumsum(probs_sort)
	n <- sum(probs_cum < 0.95) + 1 #Get 95% confidence interval (95% of probability)
	threshold <- probs_sort[n]
	cint  <- probs[probs_norm >= threshold]
	
	#Prob clonal/subclonal (one version - uses 0.9 as cutoff) 
	prob_clonal <- round(sum(probs_norm[91:100]),4)
	prob_subclonal <- round(sum(probs_norm[1:90]),4)
	
	ccf_point <- as.numeric(names(which.max(probs_norm)))
	ccf_lower <- as.numeric(names(cint)[1])
	ccf_upper <- as.numeric(names(cint)[length(cint)])
	
	res <- c(vaf_CI_lower=vafCI_lower, vaf_CI_upper=vafCI_upper, mult_point=mult_point, mult_lower=mult_lower, mult_upper=mult_upper, ccf_point=ccf_point, ccf_lower=ccf_lower, ccf_upper=ccf_upper, prob_clonal=prob_clonal, prob_subclonal=prob_subclonal)
	return(res)
}





