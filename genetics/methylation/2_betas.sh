#!/bin/bash
#SBATCH --job-name=2_betas
#SBATCH --output=2_betas.out
#SBATCH --error=2_betas.err
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
# Sets NAs based on the detection p-value and calculates marker and sample call-rates. Filters samples based on sample calle-rate
# Performs Quantile Normalisation and Calculates Beta-values.
# Written by Benjamin Lehne (Imperial College London) and Alexander Drong (Oxford University)
# Modified by Niek de Klein (niekdeklein@gmail.com)

INPUTDIR=/groups/umcg-lld/tmp03/umcg-ndeklein/EWAS/output/CPACOR

module load R/3.3.3-foss-2015b

time R --vanilla << "EOF"

require(limma)
load("anno450K.RData")
rownames(anno)=as.character(anno$Name)
cas=anno[substr(as.character(anno$Name), 1,3)=='ch.' & !(anno$CHR %in% c('X','Y')),]
cgs=anno[substr(as.character(anno$Name), 1,2)=='cg'& !(anno$CHR %in% c('X','Y')),]
auto = c(as.character(cgs$Name), as.character(cas$Name))
#auto=as.matrix(auto)

#detection p-value
thres=1E-16
load('/groups/umcg-lld/tmp03/umcg-ndeklein/EWAS/output/CPACOR/intensities.RData')
load('/groups/umcg-lld/tmp03/umcg-ndeklein/EWAS/output/CPACOR/detectionPvalue.RData')
d=dp.all[rownames(TypeII.Green.All),colnames(TypeII.Green.All)]
TypeII.Green.All.d = ifelse(d<thres,TypeII.Green.All,NA)
TypeII.Red.All.d = ifelse(d<thres,TypeII.Red.All,NA)
d=dp.all[rownames(TypeI.Green.M.All),colnames(TypeI.Green.M.All)]
TypeI.Green.M.All.d = ifelse(d<thres,TypeI.Green.M.All,NA)
TypeI.Green.U.All.d = ifelse(d<thres,TypeI.Green.U.All,NA)
d=dp.all[rownames(TypeI.Red.M.All),colnames(TypeI.Red.M.All)]
TypeI.Red.M.All.d = ifelse(d<thres,TypeI.Red.M.All,NA)
TypeI.Red.U.All.d = ifelse(d<thres,TypeI.Red.U.All,NA)
rm(dp.all,d)
#####beta incl sex chrom
print(head(TypeII.Green.All[1:3]))
print(head(TypeII.Red.All[1:3]))

TypeII.betas.d = TypeII.Green.All.d/(TypeII.Red.All.d+TypeII.Green.All.d+100)
TypeI.Green.betas.d = TypeI.Green.M.All.d/(TypeI.Green.M.All.d+TypeI.Green.U.All.d+100)
TypeI.Red.betas.d = TypeI.Red.M.All.d/(TypeI.Red.M.All.d+TypeI.Red.U.All.d+100)
beta = as.matrix(rbind(TypeII.betas.d,TypeI.Green.betas.d,TypeI.Red.betas.d))
save(beta, file="/groups/umcg-lld/tmp03/umcg-ndeklein/EWAS/output/CPACOR/beta_sex.RData")


# autosomes ------------------------------------------------------------------
samples=colnames(TypeI.Red.M.All)
category=auto
markers=as.matrix(intersect(rownames(TypeII.Green.All.d), category))
TypeII.Green = TypeII.Green.All.d[markers,samples]
TypeII.Red = TypeII.Red.All.d[markers,samples]
markers=intersect(rownames(TypeI.Green.M.All.d), category)
TypeI.Green.M = TypeI.Green.M.All.d[markers,samples]
TypeI.Green.U = TypeI.Green.U.All.d[markers,samples]
markers=intersect(rownames(TypeI.Red.M.All.d), category)
TypeI.Red.M = TypeI.Red.M.All.d[markers,samples]
TypeI.Red.U = TypeI.Red.U.All.d[markers,samples]
#save(TypeI.Red.U.All.d,TypeI.Red.M.All.d,markers,samples,auto,cas,file="DEBUG.RData") 

#raw betas
print(typeof(TypeII.Green))
print(head(TypeII.Green))
print(typeof(TypeII.Red))
print(head(TypeII.Red))
TypeII.betas = TypeII.Green/(TypeII.Red+TypeII.Green+100)
TypeI.Green.betas = TypeI.Green.M/(TypeI.Green.M+TypeI.Green.U+100)
TypeI.Red.betas = TypeI.Red.M/(TypeI.Red.M+TypeI.Red.U+100)
beta = as.matrix(rbind(TypeII.betas,TypeI.Green.betas,TypeI.Red.betas))
sample.call=colSums(!is.na(beta))/nrow(beta)
marker.call=rowSums(!is.na(beta))/ncol(beta)
save(sample.call, marker.call, file='/groups/umcg-lld/tmp03/umcg-ndeklein/EWAS/output/CPACOR/callRates.RData')
save(beta, file='/groups/umcg-lld/tmp03/umcg-ndeklein/EWAS/output/CPACOR/beta_raw.RData')

#call-rate filtering
callrate.thres=0.95
samples=names(sample.call[sample.call>callrate.thres])
markers=as.matrix(intersect(rownames(TypeII.Green.All.d), category))
TypeII.Green = TypeII.Green.All.d[markers,samples]
TypeII.Red = TypeII.Red.All.d[markers,samples]
markers=intersect(rownames(TypeI.Green.M.All.d), category)
TypeI.Green.M = TypeI.Green.M.All.d[markers,samples]
TypeI.Green.U = TypeI.Green.U.All.d[markers,samples]
markers=intersect(rownames(TypeI.Red.M.All.d), category)
TypeI.Red.M = TypeI.Red.M.All.d[markers,samples]
TypeI.Red.U = TypeI.Red.U.All.d[markers,samples]
rm(TypeII.Green.All.d,TypeII.Red.All.d,TypeI.Green.M.All.d,TypeI.Green.U.All.d,TypeI.Red.M.All.d,TypeI.Red.U.All.d)

#QN
TypeII.Green=normalizeQuantiles(TypeII.Green)
TypeII.Red = normalizeQuantiles(TypeII.Red)
TypeI.Green.M = normalizeQuantiles(TypeI.Green.M)
TypeI.Green.U = normalizeQuantiles(TypeI.Green.U)
TypeI.Red.M = normalizeQuantiles(TypeI.Red.M)
TypeI.Red.U = normalizeQuantiles(TypeI.Red.U)
TypeII.betas = TypeII.Green/(TypeII.Red+TypeII.Green+100)
TypeI.Green.betas = TypeI.Green.M/(TypeI.Green.M+TypeI.Green.U+100)
TypeI.Red.betas = TypeI.Red.M/(TypeI.Red.M+TypeI.Red.U+100)
beta = as.matrix(rbind(TypeII.betas,TypeI.Green.betas,TypeI.Red.betas))
rm(TypeII.Green,TypeII.Red,TypeI.Green.M,TypeI.Green.U,TypeI.Red.M,TypeI.Red.U,TypeII.betas,TypeI.Green.betas,TypeI.Red.betas)
save(beta, file="/groups/umcg-lld/tmp03/umcg-ndeklein/EWAS/output/CPACOR/beta_QN.RData")

EOF

