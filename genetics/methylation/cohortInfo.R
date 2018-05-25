
pheno <- read.table('/groups/umcg-lld/tmp03/umcg-ndeklein/EWAS/LLD_1135subjects_ImmuneMarkers.txt',sep="\t",header=T)
print(head(pheno))
pheno2 <- read.table('/groups/umcg-lld/tmp03/phenotype_data/LLD_data.txt',sep="\t",header=T)

print(paste('sample size:',nrow(pheno2)))
print(paste('Age(yrs):',median(pheno2$antrop_age)))
s <- table(pheno2$antrop_gender.F1M2)
percent_F <- (s[[1]]/nrow(pheno2))*100
print(paste('Sex (% female):',percent_F))
hsCRP <- pheno$hs.CRP[!is.na(pheno$hs.CRP)]
print(paste0('c-reactive protein: ',median(hsCRP),' (',sd(hsCRP),')'))
hsCRP_logged <- log(pheno$hs.CRP[!is.na(pheno$hs.CRP)]+0.01)
print(paste0('log(c-reactive protein+0.01): ',median(hsCRP_logged),' (',sd(hsCRP_logged),')'))
print(paste0('BMI (median): ',median(pheno2$antrop_BMI),' (',sd(pheno2$antrop_BMI),')'))
print(paste0('BMI (mean): ',mean(pheno2$antrop_BMI),' (',sd(pheno2$antrop_BMI),')'))
