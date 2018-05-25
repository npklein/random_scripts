#!/bin/bash
#SBATCH --job-name=1_retrieve_intensities
#SBATCH --output=1_retrieve_intensities.out
#SBATCH --error=1_retrieve_intensities.err
#SBATCH --time=5:59:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 80gb
#SBATCH --nodes 1
#SBATCH --export=NONE
#SBATCH --get-user-env=30L
#SBATCH --qos=dev
set -u
set -x
set -e
# Retrieves signal-intensities from Illumina idat-files. Returns detection p-values and signal intensities for genomic CpGs and control probes.
# Written by Benjamin Lehne (Imperial College London), Alexander Drong (Oxford University) and Reiner Schulz (King's College London)
# Modified by Niek de Klein
OUTDIR=/groups/umcg-lld/tmp03/umcg-ndeklein/EWAS/output/CPACOR
mkdir -p $OUTDIR

module load R/3.3.3-foss-2015b
module list

R --vanilla << "EOF"

setwd('/groups/umcg-lld/tmp03/umcg-ndeklein/EWAS/output/CPACOR')

require(minfi)
require(IlluminaHumanMethylation450kmanifest)

filenames <- Sys.glob("/groups/umcg-lld/tmp03/rawdata/IlluminaHumanMethylation450k/*/*_Grn.idat")


interval=2
for(i in 1:floor(length(filenames)/interval)){
    if((length(filenames)%%interval!=0) && (length(filenames)-(i*interval)<interval)){
        filenamesSubset <- filenames[(((i-1)*interval)+1):length(filenames)] 
        fileBaseNames <- unlist(strsplit(filenamesSubset,'_Grn.idat'))
    }else{
        filenamesSubset <- filenames[(((i-1)*interval)+1):(i*interval)]
        fileBaseNames <- unlist(strsplit(filenamesSubset,'_Grn.idat'))
    }
    print('read data')
    RGset <- read.metharray(fileBaseNames, verbose=TRUE)
    print('Correct background')
	RGset <- bgcorrect.illumina(RGset)  # Illumina background subtraction

	# Type II probes
    print('Type II probes')
    probeInfo <- getProbeInfo(RGset, type = "II")
	TypeII.Name <- probeInfo$Name
    TypeII.Adress <- probeInfo$Address
    green <- getGreen(RGset)
	TypeII.Green <- green[TypeII.Adress,]
    red <- getRed(RGset)
	TypeII.Red <- red[TypeII.Adress,]

    sampleNames <- sampleNames(RGset)
	rownames(TypeII.Red) <- TypeII.Name
	colnames(TypeII.Red) <- sampleNames
	rownames(TypeII.Green) <- TypeII.Name
	colnames(TypeII.Green) <- sampleNames

	# Type I probes, split into green and red channels
    print('Type I probes')
    print('green')
	TypeI.Green.Name <- getProbeInfo(RGset, type = "I-Green")$Name
	TypeI.Green.M <- getGreen(RGset)[getProbeInfo(RGset, type = "I-Green")$AddressB,]
	rownames(TypeI.Green.M) <- TypeI.Green.Name
	colnames(TypeI.Green.M) <- sampleNames(RGset)
	TypeI.Green.U <- getGreen(RGset)[getProbeInfo(RGset, type = "I-Green")$AddressA,]
	rownames(TypeI.Green.U) <- TypeI.Green.Name
	colnames(TypeI.Green.U) <- sampleNames(RGset)
    print('red')
	TypeI.Red.Name <- getProbeInfo(RGset, type = "I-Red")$Name
	TypeI.Red.M <- getRed(RGset)[getProbeInfo(RGset, type = "I-Red")$AddressB,]
	rownames(TypeI.Red.M) <- TypeI.Red.Name
	colnames(TypeI.Red.M) <- sampleNames(RGset)
	TypeI.Red.U <- getRed(RGset)[getProbeInfo(RGset, type = "I-Red")$AddressA,]
	rownames(TypeI.Red.U) <- TypeI.Red.Name
	colnames(TypeI.Red.U) <- sampleNames(RGset)

	#Control probes
    print('Control probes')
	controls=getProbeInfo(RGset, type = "Control")
	types=unique(controls$Type)
	types=types[types!='NEGATIVE']
	ctrl.names=controls[controls$Type %in% types,'ExtendedType']
	ctrl.address=controls[controls$Type %in% types,'Address']
	ctrl.Green <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(ctrl.names), dimnames = list(ctrl.names, sampleNames(RGset)))
	ctrl.Green[ctrl.names,] <- getGreen(RGset)[ctrl.address,]
	ctrl.Red <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(ctrl.names), dimnames = list(ctrl.names, sampleNames(RGset)))
	ctrl.Red[ctrl.names,] <- getRed(RGset)[ctrl.address,]
	rownames(ctrl.Green)=paste(rownames(ctrl.Green), '.grn', sep='')
	rownames(ctrl.Red)=paste(rownames(ctrl.Red), '.red', sep='')
    ctrl = rbind(ctrl.Green, ctrl.Red)

	#detection p-values
    print('detect p-values')
	dp = detectionP(RGset, type = "m+u")

	#add data for the new samples
	if(exists("TypeII.Red.All")) {
    	TypeII.Red.All <- cbind(TypeII.Red.All,TypeII.Red)
    	TypeII.Green.All <- cbind(TypeII.Green.All,TypeII.Green)
	    TypeI.Red.M.All <- cbind(TypeI.Red.M.All,TypeI.Red.M)
    	TypeI.Red.U.All <- cbind(TypeI.Red.U.All,TypeI.Red.U)
	    TypeI.Green.M.All <- cbind(TypeI.Green.M.All,TypeI.Green.M)
    	TypeI.Green.U.All <- cbind(TypeI.Green.U.All,TypeI.Green.U)
	    ctrl.all <- rbind(ctrl.all, t(ctrl))
    	dp.all <- cbind(dp.all, dp)
	}
	else {
	    TypeII.Red.All <- TypeII.Red
    	TypeII.Green.All <- TypeII.Green
	    TypeI.Red.M.All <- TypeI.Red.M
    	TypeI.Red.U.All <- TypeI.Red.U
	    TypeI.Green.M.All <- TypeI.Green.M
    	TypeI.Green.U.All <- TypeI.Green.U
	    ctrl.all <- t(ctrl)
        dp.all <- dp
	}

}

#PCA of control-probe intensities
print('pca')
pca <- prcomp(na.omit(ctrl.all))
ctrlprobes.scores = pca$x
colnames(ctrlprobes.scores) = paste(colnames(ctrlprobes.scores), '_cp', sep='')
print("red.all")
print(head(TypeII.Red.All[1:2]))
print(dim(TypeII.Red.All))
print("green.all")
print(head(TypeII.Green.All[1:2]))
print(dim(TypeII.Green.All))
save(TypeII.Red.All ,TypeII.Green.All ,TypeI.Red.M.All ,TypeI.Red.U.All ,TypeI.Green.M.All ,TypeI.Green.U.All , file="intensities.RData")
save(ctrl.all,ctrlprobes.scores, file = "ctrlprobes.RData")
save(dp.all, file = "detectionPvalue.RData")

print(paste('Saved to',getwd()))
EOF

