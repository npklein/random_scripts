load('/groups/umcg-lld/tmp03/umcg-ndeklein/EWAS/output/CPACOR/beta_QN.RData')
cpgs <- c('cg04535902', 'cg09662411', 'cg09935388', 'cg10399789', 'cg12876356','cg18146737','cg14179389','cg18316974')
beta_subset <- beta[cpgs,]
save(beta_subset, file = '/groups/umcg-lld/tmp03/umcg-ndeklein/EWAS/output/CPACOR/beta_QN_8cpg.RData')
