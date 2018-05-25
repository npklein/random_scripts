phe <- read.table("/groups/umcg-lld/tmp03/phenotype_data/LLD_data.txt",sep="\t",header=T)
phe2 <- read.table("/groups/umcg-lld/tmp03/umcg-ndeklein/EWAS/LLD_1135subjects_ImmuneMarkers.txt",
                    sep="\t",header=T,row.names=1)
phe3 <- read.table('LLD_waistCircumference.txt', sep='\t', header=T,
                    row.names=1)
pheMerged <- merge(phe, phe2, by='row.names')
rownames(pheMerged) <- pheMerged$Row.names
pheMerged$Row.names <- NULL
pheMerged <- merge(pheMerged, phe3, by='row.names')
rownames(pheMerged) <- pheMerged$Row.names
pheMerged$Row.names <- NULL


print(colnames(pheMerged))
print(summary(pheMerged$antrop_age))
print(table(pheMerged$antrop_gender.F1M2))
print(dim(pheMerged))
print(table(pheMerged$Smoking))


print(mean(pheMerged$antrop_SBP))
print(sd(mean(pheMerged$antrop_SBP))
print(mean(pheMerged$antrop_DBP))
print(sd(mean(pheMerged$antrop_DBP))
