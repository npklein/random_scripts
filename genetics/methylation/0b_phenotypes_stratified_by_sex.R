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

females <- pheMerged[pheMerged$antrop_gender.F1M2==1,]
males <- pheMerged[pheMerged$antrop_gender.F1M2==2,]


measurements <- c('antrop_BMI',
                  'waist_circumference',
                  'Biochem_TG_log',
                  'Biochem_HDL',
                  'Biochem_Glucose_log',
                  'antrop_SBP',
                  'antrop_DBP')

for(m in measurements){
    print(paste0(m,' females: ',signif(mean(females[[m]]),4),
                                ' (',signif(sd(females[[m]]),4),')'))

    print(paste0(m,' males: ',signif(mean(males[[m]]),4),
                                ' (',signif(sd(males[[m]]),4),')'))

    print(paste('pvalue:', signif(t.test(females[[m]], males[[m]])$p.val,4)))
    print('----------------')
}


# CG sites
load('/groups/umcg-lld/tmp03/umcg-ndeklein/EWAS/output/CPACOR/beta_QN.RData')
LLD_name_match <- read.table("/groups/umcg-lld/tmp03/LLDeep_Methylation_GTE.txt")
nvar = nrow(beta)
sampleNames <- as.data.frame(colnames(beta))
colnames(sampleNames) = c('V1')
sampleNames$LLD <- LLD_name_match[match(sampleNames$V1,LLD_name_match$V2),]$V1
sampleNames <- sampleNames[! is.na(sampleNames$LLD) ,]
beta <- beta[,sampleNames$V1]
colnames(beta) <- sampleNames[match(colnames(beta), sampleNames$V1),]$LLD

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

females <- rownames(pheMerged[pheMerged$antrop_gender.F1M2==1,])
males <- rownames(pheMerged[pheMerged$antrop_gender.F1M2==2,])


for(cpg in c('cg14179389','cg09935388','cg12876356','cg18316974',
             'cg09662411','cg04535902','cg18146737','cg10399789')){
    cpg_data <- beta[cpg,]
    female_data <- cpg_data[names(cpg_data) %in% females]
     male_data <- cpg_data[names(cpg_data) %in% males] 
     
    print(paste0(cpg,' females: ',signif(mean(female_data),4),
                                ' (',signif(sd(female_data),4),')'))

    print(paste0(cpg,' males: ',signif(mean(male_data),4),
                                ' (',signif(sd(male_data),4),')'))

    print(paste('pvalue:', signif(t.test(female_data, male_data)$p.val,4)))
    print('----------------')        
}
