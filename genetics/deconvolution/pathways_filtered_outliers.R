library('data.table')
library('WGCNA')
library("org.Hs.eg.db")
library(ggplot2)
##### set constants ######
significance <- 0.05
nTests <- 6555886
#######

###### read in deconGW filtered set #####
path <- "DeconGW_filtered_for\ Outliers_eQTLs"
files <- list.files(path)
decon_eQTLs <- vector(mode="list", lengt=length(files))
celltypes <- c()
x <- 0
for (file in files){
  x = x + 1
  decon_eQTLs_tmp <- read.csv(paste0(path,'/',file))
  genename <- sapply(decon_eQTLs_tmp$V1, function(x) strsplit(as.character(x), "_")[[1]][1])
  rownames(decon_eQTLs_tmp) <- genename
  pvalue <- decon_eQTLs_tmp[grep('pvalue',colnames(decon_eQTLs_tmp))]
  decon_eQTLs_tmp$pvalue_fdr <- p.adjust(pvalue[,1], method="fdr", n=nTests)
  decon_eQTLs_tmp_filtered <- decon_eQTLs_tmp[decon_eQTLs_tmp$pvalue_fdr < significance,]
  decon_eQTLs[[x]] <- rownames(decon_eQTLs_tmp_filtered)
  celltype <- strsplit(file, "_")[[1]][1]
  celltypes <- c(celltypes, celltype)
  print(celltype)
  print(nrow(decon_eQTLs_tmp_filtered))
  print('-----')
}
names(decon_eQTLs) <- celltypes
#####

##### Read in expression data #####
expression_data <- data.frame(fread('/Users/NPK/UMCG/projects/deconvolution/expData/expression_all_genes.txt'))
rownames(expression_data) <- expression_data$V1
expression_data$V1 <- NULL
#####

##### Read in full cell name data #####
itInfo <- read.csv('itDecon_extra_colors.csv')
#####

##### PLot number of genes per celltype #####
celltypeLength <- as.data.frame(do.call(rbind, lapply(decon_eQTLs, function(x) length(!is.na(x)))))
colnames(celltypeLength) <- c('numberOfGenes')
celltypeLength$cellCode <- rownames(celltypeLength)
celltypeLength$celltype <- itInfo[match(celltypeLength$cellCode,itInfo$code),]$finalName
celltypeLength$colorCat <- factor(itInfo[match(celltypeLength$cellCode,itInfo$code),]$ColorCat,levels= c("Myeloid", "Lymphocytes" , "T cells","B cells", "NK cells"))

p <- ggplot(celltypeLength, aes(x=celltype,y=numberOfGenes))+
  geom_bar(stat = "identity")+
  facet_grid(~colorCat, space = 'free', scales='free')+
  theme_bw()+ ylab("")+xlab("")+
  theme(legend.position ="right",
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 9.5, family = "Helvetica"),
        strip.text.x = element_blank(),
        strip.text.y = element_text(angle = 0, hjust = 0, size= 9.5),
        strip.background = element_blank())
ggsave('figures_filteredSamples/unique_genes_per_celltype.png',plot=p, width=14, height =8)
#####


##### WGNCA soft threshold function #####
soft_threshold <- function(expression_data_subset,name){
  options(stringsAsFactors = FALSE);
  # Choose a set of soft-thresholding powers
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  # Call the network topology analysis function
  sft = pickSoftThreshold(expression_data_subset, powerVector = powers, verbose = 5)
  # Plot the results:
  sizeGrWindow(9, 5)
  pdf(paste0('figures_filteredSamples/clustering/',name,'_sof_threshold.pdf'), width=8, height=5)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.90,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  dev.off()
}
#####

##### WGNCA network #####
wgnca_network <- function(expression_data_subset, name, power){
  net = blockwiseModules(expression_data_subset, power = power,
                         TOMType = "unsigned", minModuleSize = 10,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = paste0('data/',name,'_TOM'),
                         verbose = 3)
  # open a graphics window
  sizeGrWindow(12, 9)
  pdf(paste0('figures_filteredSamples/clustering/',name,'_dendogram.pdf'), width=8, height=5)
  # Convert labels to colors for plotting
  mergedColors = labels2colors(net$colors)
  # Plot the dendrogram and the module colors underneath
  plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
  moduleLabels = net$colors
  moduleColors = labels2colors(net$colors)
  MEs = net$MEs;
  geneTree = net$dendrograms[[1]];
  save(MEs, moduleLabels, moduleColors, geneTree,
       file = paste0('data/',name,'_networkConstructionAuto.RData'))
  return(net)
}
#####

##### WGNCA go enrichment function #####
go_enrich <- function(net,genes,name){
  moduleColors = labels2colors(net$colors)
  genes_unlisted <- sapply(genes, `[[`, 1)
  GOenr = GOenrichmentAnalysis(moduleColors, genes_unlisted, organism = "human", nBestP = 10);
  tab = GOenr$bestPTerms[[4]]$enrichment
  tab$name <- name
  write.table(tab, file = paste0('goTable_filteredSamples/',name,"_GOEnrichmentTable.csv"), sep = "\t", quote = FALSE, row.names = FALSE)
  for (module in unique(tab$module)){
    tab_subset <- tab[tab$module==module,]
    tab_subset <- tab_subset[c('termID','BonferoniP')]
    tab_subset <- tab_subset[tab_subset$BonferoniP < 0.05,]
    if(nrow(tab_subset) > 0){
      write.table(tab_subset, file = paste0('modules_filteredSamples/',name,"_",module,"_GOtermsPvalue.csv"), sep = "\t", quote = FALSE, row.names = FALSE)
    }
  }
}
#####
converter <- as.list(org.Hs.egENSEMBL2EG)
tablist = list()
i = 0 
expression_data_t <- as.data.frame(t(expression_data))
for (name in names(celltypes)){
  i = i + 1
  print(name)
  print(length(celltypes[[name]]))
  matchedGenes <- match(celltypes[[name]],colnames(expression_data_t))
  matchedGenes <- matchedGenes[!is.na(matchedGenes)]
  expression_data_subset <- expression_data_t[matchedGenes]
  soft_threshold(expression_data_subset,name)
}

# these powers are decided from soft threshold figures made above, needs to be inspected and changed by hand
##### powers #####
power <- data.frame('IT83'=3,
                    'IT73'=3,
                    'IT70'=2,
                    'IT66'=1,
                    'IT65'=2,
                    'IT63'=3,
                    'IT61'=2,
                    'IT60'=2,
                    'IT59'=4,
                    'IT56'=3,
                    'IT51'=3,
                    'IT50'=3,
                    'IT34'=4,
                    'IT33'=5,
                    'IT21'=2,
                    'IT17'=2,
                    'IT15'=1,
                    'IT13'=2,
                    'IT12'=2,
                    'IT11'=2,
                    'IT8'=2,
                    'IT7'=2,
                    'IT5'=2,
                    'IT2'=3,
                    'IT1'=4)
######

for (name in names(celltypes)){
  i = i + 1
  print(name)
  print(length(celltypes[[name]]))
  matchedGenes <- match(celltypes[[name]],colnames(expression_data_t))
  matchedGenes <- matchedGenes[!is.na(matchedGenes)]
  expression_data_subset <- expression_data_t[matchedGenes]
  # power is chosen by inspecting soft_threshold graphs in figures/
  net <- wgnca_network(expression_data_subset, name, power[[toupper(name)]])
  genes <- converter[colnames(expression_data_subset)]
  tab <- go_enrich(net, genes,name)
}

filenames <- list.files("goTables/", pattern="*.csv", full.names=TRUE)
go_info <- lapply(filenames, fread)

go_info_combined <- do.call(rbind,lapply(go_info,data.frame))

go_info_relevant <- go_info_combined[c('termName','BonferoniP','name', 'termOntology','modSize','rank')]
go_info_relevant <- go_info_relevant[go_info_relevant$BonferoniP<0.05,]
go_info_relevant <- go_info_relevant[go_info_relevant$termOntology == 'BP',]




go_info_relevant$celltype <- itInfo[match(go_info_relevant$name,tolower(itInfo$code)),]$finalName
go_info_relevant$colorCat <- factor(itInfo[match(go_info_relevant$name,tolower(itInfo$code)),]$ColorCat,levels= c("Myeloid", "Lymphocytes" , "T cells","B cells", "NK cells"))

go_info_relevant_data <- go_info_relevant[c('termName','celltype','BonferoniP')]


orderTrait <- aggregate(-log10(go_info_relevant_data$BonferoniP), by=list(go_info_relevant_data$termName), FUN=sum)
go_info_relevant$termName <- factor(go_info_relevant$termName, levels=orderTrait[order(orderTrait$x, decreasing = FALSE),1])

go_info_relevant <- go_info_relevant[go_info_relevant$rank <= 5,]

# filled in by hand the general term
#write.table(unique(go_info_relevant$termName), 'terms2.csv', quote=F, row.names=F, col.names=F)
termGeneral <- read.table('terms.csv',sep="\t", header=T)
go_info_relevant$generalTerm <- termGeneral[match(go_info_relevant$termName,termGeneral$term),]$general
go_info_relevant$generalTerm <- factor(go_info_relevant$generalTerm, levels = c('immune', 'sensory', 'signaling', 'metabolic', 'cell cycle', 'differentiation', 'other'))
#geom_tile(aes(fill = -log10(BonferoniP))) +
library(viridis)
p <- ggplot(go_info_relevant, aes(x=celltype, y=termName)) +
  geom_tile(aes(fill = -log10(BonferoniP)), color='white') +
  facet_grid(generalTerm~colorCat, space = 'free', scales='free')+
  scale_fill_viridis()+
  theme_bw()+ ylab("")+xlab("")+
  theme(legend.position ="right",
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 9.5, family = "Helvetica"),
        strip.text.x = element_blank(),
        strip.text.y = element_text(angle = 0, hjust = 0, size= 9.5),
        strip.background = element_blank())

save(go_info_relevant,p, file="go_plot_and_data.RData")
ggsave('figures/goFunctions.pdf',height=12, width=10)

