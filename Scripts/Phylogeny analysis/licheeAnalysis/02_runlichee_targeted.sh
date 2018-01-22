#############

indir=~/Documents/Caldas/Data/dataProcessed/lichee/targetedinput/
resultsdir=~/Documents/Caldas/Data/dataResults/lichee/targeted/
licheedir=~/Documents/Caldas/Analysis/lichee-master/LICHeE/release

#mkdir -p ${resultsdir}

cd ${licheedir}

for case in 288 290 298 308 315 323 328 330 DET52; do

	./lichee -build -i ${indir}/${case}_cp.txt -cp -clustersFile ${indir}/${case}_clusters.txt -maxVAFAbsent 0.01 -minVAFPresent 0.05 -n 0 -color -e 0.3 -dot -dotFile ${resultsdir}/${case}.dot

	dot -Tpdf ${resultsdir}/${case}.dot -O

done
