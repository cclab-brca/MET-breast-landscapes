#!/bin/bash

################################################################################
##### input files needed #####
## Output from 01_prepareASCATinput.sh

##### programs and scripts needed #####
## rerunASCAT.R

################################################################################
##### --- Adjust as needed --- #####
## Data directories
rawdatadir=~/Caldas/Data/dataRaw/
intdatadir=~/Caldas/Data/dataProcessed/ascat/
referencedir=~/Caldas/Data/dataAnnotation/Reference/

## Program locations
ascatdir=~/Caldas/Programs/ascatNGS/
scriptdir=~/Caldas/Data/breastMetRepo/mediccAnalysis

################################################################################

#cp -R ${intdatadir}/ascatOutput/* ${intdatadir}/ascatRefitted/

### 288 ###
case=288
baseCase=288
cd ${intdatadir}/ascatRefitted/${case}

sample=288-016
rho=0.53
psi=3.5
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}

sample=288-017
rho=0.27 #0.35
psi=3.2 #2.1
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}

sample=288-020
rho=0.6
psi=3.1
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}



### 290exOvary ###
case=290exOvary
cd ${intdatadir}/ascatRefitted/${case}

## exclude 005 (something wrong, total CN flat)
rm ${intdatadir}/ascatRefitted/${case}/290-005*

sample=290-007
rho=0.85
psi=4.1
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}

sample=290-008
rho=0.72
psi=4.1
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}

sample=290-014
rho=0.68
psi=4
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}

sample=290-016-B-IDC
rho=0.7
psi=4.1
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}

sample=290-016-B-muc
rho=0.9
psi=4.4
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}

sample=290-016-B-WT
rho=0.9
psi=4
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}

sample=290-018
rho=0.85
psi=4.1
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}

sample=290-021
rho=0.75
psi=4.2
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}


### 298main ###
case=298main
cd ${intdatadir}/ascatRefitted/${case}
 
sample=298-004
rho=0.98
psi=4.2
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}

sample=298-012
rho=0.96
psi=4.1
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}
 
sample=298-020
rho=0.99
psi=4.35
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}

sample=298-022
rho=0.99
psi=4.2
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}


### 308 ###
case=308
cd ${intdatadir}/ascatRefitted/${case}

sample=308-003
rho=0.75
psi=3.7
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}

sample=308-004
rho=0.75
psi=3.2
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}

sample=308-006
rho=0.75
psi=3.1
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}

sample=308-008
rho=0.75
psi=3.6
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}

sample=308-010
rho=0.72
psi=3.3
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}

sample=308-014
rho=0.63
psi=3.15
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}

sample=308-015
rho=0.75
psi=3.2
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}

sample=308-017
rho=0.75
psi=3.15
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}

sample=308-018
rho=0.43
psi=3.1
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}

sample=308-019
rho=0.75
psi=3.1
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}

sample=308-020
rho=0.8
psi=3.1
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}

sample=308-021
rho=0.64
psi=3.2
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}

sample=308-022
rho=0.6
psi=3.1
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}


### 315 ###
case=315
cd ${intdatadir}/ascatRefitted/${case}

sample=315-001
rho=0.5
psi=3.9
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}

sample=315-002
rho=0.85
psi=3.9
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}

sample=315-003
rho=0.8
psi=3.9
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}

sample=315-004
rho=0.77
psi=3.9
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}

sample=315-005
rho=0.8
psi=4.0
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}

sample=315-007
rho=0.68
psi=4.0
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}

sample=315-009
rho=0.7
psi=4.0
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}

sample=315-012
rho=0.65
psi=3.9
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}

sample=315-014
rho=0.55
psi=3.9
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}

sample=315-015
rho=0.46
psi=3.9
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}

rm ${intdatadir}/ascatRefitted/${case}/315-017*
rm ${intdatadir}/ascatRefitted/${case}/315-018*



### 323 ###
case=323
cd ${intdatadir}/ascatRefitted/${case}

sample=323-003
rho=0.3
psi=4.1
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}


### 328 ###
case=328
cd ${intdatadir}/ascatRefitted/${case}

sample=328-004
rho=0.99
psi=3.1
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}

sample=328-005
rho=0.94
psi=3.1
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}


### 330 ###
case=330
cd ${intdatadir}/ascatRefitted/${case}

sample=330-002
rho=0.77
psi=3.6
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}

sample=330-003
rho=0.75
psi=4.0
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}

sample=330-004
rho=0.28
psi=2.5
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}

sample=330-005
rho=0.55
psi=3.8
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}


### DET52 ###
case=DET52
cd ${intdatadir}/ascatRefitted/${case}

sample=DET52-mt3
rho=0.44
psi=4.2
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}

sample=DET52-mt5
rho=0.4
psi=4.0
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}

sample=DET52-mt7
rho=0.25
psi=4.15
/usr/bin/Rscript ${scriptdir}/02_rerunASCAT.R ${sample}.merged.realn.RG.bam XX 24 ${intdatadir}/ascatRefitted/${case}/${sample}.RData ${rho} ${psi}


## exclude mt8 (something wrong with fit, can find good parameters)
rm ${intdatadir}/ascatRefitted/${case}/DET52-mt8*
