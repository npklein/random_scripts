#!/bin/bash
# Written by ?, modified by Niek
module load R/3.3.3-foss-2015b


time R --vanilla << "EOF"
path.out="/groups/umcg-lld/tmp03/umcg-ndeklein/EWAS/output/"
print(path.out)
phe <- read.table("/groups/umcg-lld/tmp03/phenotype_data/LLD_data.txt",sep="\t",header=T)
phe2 <- read.table("/groups/umcg-lld/tmp03/umcg-ndeklein/EWAS/LLD_1135subjects_ImmuneMarkers.txt",
                    sep="\t",header=T,row.names=1)
cellcounts <- read.table("/groups/umcg-lld/tmp03/umcg-ndeklein/EWAS/output/CPACOR/cellcounts.txt",header=T)
LLD_name_match <- read.table("/groups/umcg-lld/tmp03/LLDeep_Methylation_GTE.txt")
cellcounts$LLD <- LLD_name_match[match(row.names(cellcounts),LLD_name_match$V2),]$V1
cellcounts <- cellcounts[! is.na(cellcounts$LLD),]
rownames(cellcounts) <- cellcounts$LLD
cellcounts$LLD <- NULL
pheMerged <- merge(phe, phe2, by='row.names')
pheMerged$antrop_gender.F1M2 <- pheMerged$antrop_gender.F1M2-1
rownames(pheMerged) <- pheMerged$Row.names
pheMerged$Row.names <- NULL
pheMerged <- merge(pheMerged,cellcounts, by="row.names")
pheMerged$logCRP <- log(pheMerged$hs.CRP+0.01)
load('/groups/umcg-lld/tmp03/umcg-ndeklein/EWAS/output/CPACOR/ctrlprobes.RData')
ctrlprobes.scores <- as.data.frame(ctrlprobes.scores)
ctrlprobes.scores$LLD <- LLD_name_match[match(rownames(ctrlprobes.scores),LLD_name_match$V2),]$V1
ctrlprobes.scores <- ctrlprobes.scores[!is.na(ctrlprobes.scores$LLD),]
rownames(ctrlprobes.scores) <- ctrlprobes.scores$LLD
ctrlprobes.scores$LLD <- NULL
load('/groups/umcg-lld/tmp03/umcg-ndeklein/EWAS/output/CPACOR/beta_QN.RData')
nvar = nrow(beta)
sampleNames <- as.data.frame(colnames(beta))
colnames(sampleNames) = c('V1')
sampleNames$LLD <- LLD_name_match[match(sampleNames$V1,LLD_name_match$V2),]$V1
sampleNames <- sampleNames[! is.na(sampleNames$LLD) ,]
beta <- beta[,sampleNames$V1]
colnames(beta) <- sampleNames[match(colnames(beta), sampleNames$V1),]$LLD
#ctrlprobes.scores
#----------------------------------------------------------
rownames(pheMerged) <- pheMerged$Row.names
pheMerged$Row.names <- NULL
pheMerged <- merge(pheMerged, ctrlprobes.scores, by='row.names')
pheMerged <- pheMerged[match(colnames(beta), pheMerged$Row.names),]
model=c(paste('pheMerged[,"logCRP"] ~ beta[i, ] +',
                          'pheMerged[,"antrop_gender.F1M2"] +',
                          'pheMerged[,"antrop_age"] +',
                          'pheMerged[,"CD8T"] +',
                          'pheMerged[,"CD4T"] +',
                          'pheMerged[,"NK"] +',
                          'pheMerged[,"Bcell"] +',
                          'pheMerged[,"Mono"] +',
                          'pheMerged[,"Gran"] +',
                          'pheMerged[,"PC1_cp"] + pheMerged[,"PC2_cp"] + pheMerged[,"PC3_cp"] + pheMerged[,"PC4_cp"] + pheMerged[,"PC5_cp"] + pheMerged[,"PC6_cp"] + pheMerged[,"PC7_cp"] + pheMerged[,"PC8_cp"] + pheMerged[,"PC9_cp"] + pheMerged[,"PC10_cp"] + pheMerged[,"PC11_cp"] + pheMerged[,"PC12_cp"] + pheMerged[,"PC13_cp"] + pheMerged[,"PC14_cp"] + pheMerged[,"PC15_cp"] + pheMerged[,"PC16_cp"] + pheMerged[,"PC17_cp"] + pheMerged[,"PC18_cp"] + pheMerged[,"PC19_cp"] + pheMerged[,"PC20_cp"] + pheMerged[,"PC21_cp"] + pheMerged[,"PC22_cp"] + pheMerged[,"PC23_cp"] + pheMerged[,"PC24_cp"] + pheMerged[,"PC25_cp"] + pheMerged[,"PC26_cp"] + pheMerged[,"PC27_cp"] + pheMerged[,"PC28_cp"] + pheMerged[,"PC29_cp"] + pheMerged[,"PC30_cp"]'))
lfla=as.formula(model)

# regression
res=matrix(ncol=7, nrow=nvar)

for(i in 1:nvar) {
  if (i %% 10000 == 0){
      print(paste0(i,'/',nvar))
  }
  fit= lm(lfla)
  nsamp=nobs(fit)
  m = mean(beta[i,])
  s = sd(beta[i,])
  res[i,] = c(summary(fit)$coefficients[2,],nsamp,m,s)
  rm(fit)
}
res <- as.data.frame(res)
res$probeID <- rownames(beta)[1:nvar]
colnames(res) <- c('BETA','SE','t-val','P_VAL','N_samp','mean','SD','probeID')
res <- res[c('probeID','BETA','SE','P_VAL','N_samp','mean','SD')]
d <- format(Sys.Date(),"%d%m%Y")
out <- paste0(path.out,'CRP_model1_LLD_',d,'.txt')
write.table(res,out,sep="\t",quote=F,row.names=FALSE)
print(paste0('written to ',out))
rm(res)
EOF

#-------------------------------------------------------------------------------
#echo -e "ID\tEstimate\tSE\tt-val\tP-val\tN" | cat - ${path}result_trait.txt > $TMPDIR/tmpfile; mv $TMPDIR/tmpfile ${path}result_trait.txt




