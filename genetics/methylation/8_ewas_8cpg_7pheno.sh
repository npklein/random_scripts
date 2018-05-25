# Written by ?, modified by Niek
#!/bin/bash
module load R/3.3.3-foss-2015b


time R --vanilla << "EOF"
path.out="/groups/umcg-lld/tmp03/umcg-ndeklein/EWAS/output/"
print(path.out)
phe <- read.table("/groups/umcg-lld/tmp03/phenotype_data/LLD_data.txt",sep="\t",header=T)
phe2 <- read.table("/groups/umcg-lld/tmp03/umcg-ndeklein/EWAS/LLD_1135subjects_ImmuneMarkers.txt",
                    sep="\t",header=T,row.names=1)
phe3 <- read.table('LLD_waistCircumference.txt', sep='\t', header=T,
                    row.names=1)
cellcounts <- read.table("/groups/umcg-lld/tmp03/umcg-ndeklein/EWAS/output/CPACOR/cellcounts.txt",header=T)
LLD_name_match <- read.table("/groups/umcg-lld/tmp03/LLDeep_Methylation_GTE.txt")
cellcounts$LLD <- LLD_name_match[match(row.names(cellcounts),LLD_name_match$V2),]$V1
cellcounts <- cellcounts[! is.na(cellcounts$LLD),]
rownames(cellcounts) <- cellcounts$LLD
cellcounts$LLD <- NULL
pheMerged <- merge(phe, phe2, by='row.names')
rownames(pheMerged) <- pheMerged$Row.names
pheMerged$Row.names <- NULL
pheMerged <- merge(pheMerged, phe3, by='row.names')
pheMerged$antrop_gender.F1M2 <- pheMerged$antrop_gender.F1M2-1
rownames(pheMerged) <- pheMerged$Row.names
pheMerged$Row.names <- NULL
pheMerged <- merge(pheMerged,cellcounts, by="row.names")
load('/groups/umcg-lld/tmp03/umcg-ndeklein/EWAS/output/CPACOR/ctrlprobes.RData')
ctrlprobes.scores <- as.data.frame(ctrlprobes.scores)
ctrlprobes.scores$LLD <- LLD_name_match[match(rownames(ctrlprobes.scores),LLD_name_match$V2),]$V1
ctrlprobes.scores <- ctrlprobes.scores[!is.na(ctrlprobes.scores$LLD),]
rownames(ctrlprobes.scores) <- ctrlprobes.scores$LLD
ctrlprobes.scores$LLD <- NULL
load('/groups/umcg-lld/tmp03/umcg-ndeklein/EWAS/output/CPACOR/beta_QN_8cpg.RData')
nvar = nrow(beta_subset)
sampleNames <- as.data.frame(colnames(beta_subset))
colnames(sampleNames) = c('V1')
sampleNames$LLD <- LLD_name_match[match(sampleNames$V1,LLD_name_match$V2),]$V1
sampleNames <- sampleNames[! is.na(sampleNames$LLD) ,]
beta_subset <- beta_subset[,sampleNames$V1]
colnames(beta_subset) <- sampleNames[match(colnames(beta_subset), sampleNames$V1),]$LLD
#ctrlprobes.scores
#----------------------------------------------------------
rownames(pheMerged) <- pheMerged$Row.names
pheMerged$Row.names <- NULL
pheMerged <- merge(pheMerged, ctrlprobes.scores, by='row.names')
pheMerged <- pheMerged[match(colnames(beta_subset), pheMerged$Row.names),]

phenotypesToCalculate <- c('antrop_BMI',
                           'waist_circumference',
                           'Biochem_TG_log',
                           'Biochem_Glucose_log',
                           'Biochem_HDL',
                           'antrop_SBP',
                           'antrop_DBP')
print(dim(pheMerged))
pheno_res1 <- data.frame()
pheno_res2 <- data.frame()
for(phenotypeName in phenotypesToCalculate){
    print(phenotypeName)
    phenotype = scale(pheMerged[,phenotypeName])
    print(paste0('mean: ',mean(phenotype, na.rm=T)))
    print(paste0('sd: ',sd(phenotype, na.rm=T)))
    model1 = c(paste('phenotype ~ beta_subset[i, ] +',
                              'pheMerged[,"CD8T"] +',
                              'pheMerged[,"CD4T"] +',
                              'pheMerged[,"NK"] +',
                              'pheMerged[,"Bcell"] +',
                              'pheMerged[,"Mono"] +',
                              'pheMerged[,"Gran"] +',
                              'pheMerged[,"PC1_cp"] + pheMerged[,"PC2_cp"] + pheMerged[,"PC3_cp"] + pheMerged[,"PC4_cp"] + pheMerged[,"PC5_cp"] + pheMerged[,"PC6_cp"] + pheMerged[,"PC7_cp"] + pheMerged[,"PC8_cp"] + pheMerged[,"PC9_cp"] + pheMerged[,"PC10_cp"] + pheMerged[,"PC11_cp"] + pheMerged[,"PC12_cp"] + pheMerged[,"PC13_cp"] + pheMerged[,"PC14_cp"] + pheMerged[,"PC15_cp"] + pheMerged[,"PC16_cp"] + pheMerged[,"PC17_cp"] + pheMerged[,"PC18_cp"] + pheMerged[,"PC19_cp"] + pheMerged[,"PC20_cp"] + pheMerged[,"PC21_cp"] + pheMerged[,"PC22_cp"] + pheMerged[,"PC23_cp"] + pheMerged[,"PC24_cp"] + pheMerged[,"PC25_cp"] + pheMerged[,"PC26_cp"] + pheMerged[,"PC27_cp"] + pheMerged[,"PC28_cp"] + pheMerged[,"PC29_cp"] + pheMerged[,"PC30_cp"]'))
    model2 = c(paste('phenotype ~ beta_subset[i, ] +',
                              'pheMerged[,"antrop_gender.F1M2"] +',
                              'pheMerged[,"antrop_age"] +',
                              'pheMerged[,"Smoking"] + ',
                              'pheMerged[,"CD8T"] +',
                              'pheMerged[,"CD4T"] +',
                              'pheMerged[,"NK"] +',
                              'pheMerged[,"Bcell"] +',
                              'pheMerged[,"Mono"] +',
                              'pheMerged[,"Gran"] +',
                              'pheMerged[,"PC1_cp"] + pheMerged[,"PC2_cp"] + pheMerged[,"PC3_cp"] + pheMerged[,"PC4_cp"] + pheMerged[,"PC5_cp"] + pheMerged[,"PC6_cp"] + pheMerged[,"PC7_cp"] + pheMerged[,"PC8_cp"] + pheMerged[,"PC9_cp"] + pheMerged[,"PC10_cp"] + pheMerged[,"PC11_cp"] + pheMerged[,"PC12_cp"] + pheMerged[,"PC13_cp"] + pheMerged[,"PC14_cp"] + pheMerged[,"PC15_cp"] + pheMerged[,"PC16_cp"] + pheMerged[,"PC17_cp"] + pheMerged[,"PC18_cp"] + pheMerged[,"PC19_cp"] + pheMerged[,"PC20_cp"] + pheMerged[,"PC21_cp"] + pheMerged[,"PC22_cp"] + pheMerged[,"PC23_cp"] + pheMerged[,"PC24_cp"] + pheMerged[,"PC25_cp"] + pheMerged[,"PC26_cp"] + pheMerged[,"PC27_cp"] + pheMerged[,"PC28_cp"] + pheMerged[,"PC29_cp"] + pheMerged[,"PC30_cp"]'))
    lfla1=as.formula(model1)
    lfla2=as.formula(model2)

    # regression
    res1=matrix(ncol=7, nrow=nvar)
    res2=matrix(ncol=7, nrow=nvar)

    for(i in 1:nvar) {
        print(paste0(i,'/',nvar))
        fit1 = lm(lfla1)
        fit2 = lm(lfla2)
        nsamp1=nobs(fit1)
        nsamp2=nobs(fit2)
        m = mean(beta_subset[1,])
        s = sd(beta_subset[1,])
        x <- tryCatch({
    	    res1[i,] = c(summary(fit1)$coefficients[2,],nsamp1,m,s)
    	    res2[i,] = c(summary(fit2)$coefficients[2,],nsamp2,m,s)
        }, warning=function(w) {
    	    res1[i,] = c(summary(fit1)$coefficients[2,],nsamp1,m,s)
    	    res2[i,] = c(summary(fit2)$coefficients[2,],nsamp2,m,s)
            message("handling warning: ", conditionMessage(w))
        })
        rm(fit1)
        rm(fit2)
    }
    res1 <- as.data.frame(res1)
    res2 <- as.data.frame(res2)
    res1$probeID <- rownames(beta_subset)[1:nvar]
    res2$probeID <- rownames(beta_subset)[1:nvar]
    colnames(res1) <- c('BETA','SE','t-val','P_VAL','N_samp','mean','SD','probeID')
    colnames(res2) <- c('BETA','SE','t-val','P_VAL','N_samp','mean','SD','probeID')
    res1 <- res1[c('probeID','BETA','SE','P_VAL','N_samp','mean','SD')]
    res2 <- res2[c('probeID','BETA','SE','P_VAL','N_samp','mean','SD')]
    res1$pheno <- phenotypeName
    res2$pheno <- phenotypeName
    pheno_res1 <- rbind(pheno_res1, res1)
    pheno_res2 <- rbind(pheno_res2, res2)
    rm(res1)
    rm(res2)
}
d <- format(Sys.Date(),"%d%m%Y")
out1 <- paste0(path.out,'7pheno_8cg_ewas_model1_cpg_technic_cov_',d,'.txt')
out2 <- paste0(path.out,'7pheno_8cg_ewas_model2_cpg_technic_cov_age_sex_smoking_',d,'.txt')
write.table(pheno_res1, out1, sep="\t",quote=F,row.names=FALSE)
write.table(pheno_res2, out2, sep="\t",quote=F,row.names=FALSE)
print(paste0('written to ',out1))
print(paste0('written to ',out2))

EOF

#-------------------------------------------------------------------------------
#echo -e "ID\tEstimate\tSE\tt-val\tP-val\tN" | cat - ${path}result_trait.txt > $TMPDIR/tmpfile; mv $TMPDIR/tmpfile ${path}result_trait.txt




