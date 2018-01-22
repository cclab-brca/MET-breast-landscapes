#############
## Cluster filtering criteria
## - At least 3 mutations
## - Cellular prevalence > 0.1 in at least 1 sample
## - Binary profile: VAF > 0.01 in at leat 40% of samples
## - Keep 1 cluster that is present in all samples (pick the one with highest cellular prevalence)
## NOTE THAT LICHeE is NOT completely DETERMINISTIC! (can produce different clone selection on different runs!!!)
## -maxVAFAbsent and -minVAFPresent required but do not seem to influence result in case of clustered input

indir=~/Documents/Caldas/Data/dataProcessed/lichee/WESinput/
resultsdir=~/Documents/Caldas/Data/dataResults/lichee/WES/
licheedir=~/Documents/Caldas/Analysis/lichee-master/LICHeE/release

#mkdir -p ${resultsdir}

cd ${licheedir}

for case in 288 290 298_en 298All 308 315 323 328 330 DET52; do

	./lichee -build -i ${indir}/${case}_cp.txt -cp -clustersFile ${indir}/${case}_clusters.txt -maxVAFAbsent 0.01 -minVAFPresent 0.05 -n 0 -color -e 0.3 -dot -dotFile ${resultsdir}/${case}.dot

	dot -Tpdf ${resultsdir}/${case}.dot -O

done
