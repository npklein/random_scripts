load('/groups/umcg-lld/tmp03/umcg-ndeklein/EWAS/output/CPACOR/beta_QN.RData')
nvar = nrow(beta)
sampleNames <- as.data.frame(colnames(beta))
colnames(sampleNames) = c('V1')
sampleNames$LLD <- LLD_name_match[match(sampleNames$V1,LLD_name_match$V2),]$V1
sampleNames <- sampleNames[! is.na(sampleNames$LLD) ,]
beta <- beta[,sampleNames$V1]
colnames(beta) <- sampleNames[match(colnames(beta), sampleNames$V1),]$LLD


