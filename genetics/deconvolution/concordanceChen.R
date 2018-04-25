# Author: Raúl Aguire
# Modified by: Niek

library(data.table)
library(DT)
library(caret)
library(ggplot2)
library(pheatmap)
library(reshape2)
library(ggrepel)
library(ggExtra)
library(gridExtra)
library(viridis)
library(RColorBrewer)
library(ggsci)
library(UpSetR)
####### allelic concordance plot function ######
allelic.conconcordance.plot <- function(eqtls, 
                                        concordance, 
                                        decon.interaction.beta,
                                        chens.col, 
                                        chen.sig.rows, 
                                        x.label.plot="Decon interaction term", 
                                        decon.color= "red", 
                                        y.label.plot= "Chens z-score")
  {
  title <- paste0(
    paste0("Concordance= ", round(concordance, digits = 2)*100, "%"), "\n",
    paste0("n.decon.sig= ", nrow(eqtls), "\n"),
    paste0("n.chen= ", sum(!is.na(eqtls[,chens.col])))
  )
  
  eqtls$color <- as.factor(sign(eqtls[,"zScoreSwapped"]) == sign(eqtls[,chens.col]))
  
  # grey color for non concordant eqtls
  plot.colors <- c("#696969", decon.color)
  names(plot.colors) <- levels(eqtls$color)
  
  plot <- ggplot(eqtls[chen.sig.rows,], aes_string(x= "zScoreSwapped", y= chens.col))+
    geom_point(size=0.7,alpha=0.4, aes(color=color))+
    geom_hline(alpha=0.7, yintercept = 0)+
    geom_vline(alpha=0.7, xintercept = 0)+
    ggtitle(title)+
    xlab(x.label.plot)+
    ylab(y.label.plot)+
    scale_color_manual(values = plot.colors)+
    theme_minimal()+ 
    theme(legend.position = "none", 
          text= element_text(family = "Helvetica", size=8),
          title = element_text(size=8))
  return(plot)
}
#########


###### calculating chens concordance #######
chens.decon.concordance.proxyEqtls <- function(eqtls, 
                                               chens, 
                                               decon.CT, 
                                               decon.CT.name,
                                               chen.sig.threshold= 0.05){

   # load data from Chen's
  eqtls$chens.neut <- chens[["neut"]][eqtls$qtl, "beta"] / chens[["neut"]][eqtls$qtl, "std.error_of_beta"]
  eqtls$chens.neut.p <- chens[["neut"]][eqtls$qtl, "p.value"]
  eqtls$chens.mono <- chens[["mono"]][eqtls$qtl, "beta"] / chens[["mono"]][eqtls$qtl, "std.error_of_beta"]
  eqtls$chens.mono.p <- chens[["mono"]][eqtls$qtl, "p.value"]
  eqtls$chens.tcell <- chens[["tcel"]][eqtls$qtl, "beta"] / chens[["tcel"]][eqtls$qtl, "std.error_of_beta"]
  eqtls$chens.tcell.p <- chens[["tcel"]][eqtls$qtl, "p.value"]
  
  # hardcoded neut since all of them are aligned in the same way
  eqtls$chen.ref <- chen.var.info[eqtls$qtl, "ref_allele"]
  eqtls$chen.alt <- chen.var.info[eqtls$qtl, "alt_allele"]
  
  eqtls$zScoreSwapped <- ifelse(eqtls$alt_allele == eqtls$chen.alt, eqtls$zScore, eqtls$zScore*-1)

  ## calcultae overall allelic concordance on: 
  # calculate concordance before log modulus transformation
  # eQTLs that pass sig.threhold 
  # that are also present in Chens eQTLs. 
  neut_specific_eqtls <- eqtls[eqtls$Neutrophils=='+',]
  neut.chen.sig <- which(neut_specific_eqtls$chens.neut.p <= chen.sig.threshold)
  neut.concordance <- (sum(sign(neut_specific_eqtls$zScoreSwapped) == sign(neut_specific_eqtls$chens.neut), na.rm = TRUE)) / sum(!is.na(neut_specific_eqtls$chens.neut))

  mono_specific_eqtls <- eqtls[eqtls$Monocytes=='+',]
  mono.chen.sig <- which(mono_specific_eqtls$chens.mono.p <= chen.sig.threshold)
  mono.concordance <- (sum(sign(mono_specific_eqtls$zScoreSwapped) == sign(mono_specific_eqtls$chens.mono), na.rm = TRUE)) / sum(!is.na(mono_specific_eqtls$chens.mono))
  
  cd4_specific_eqtls <- eqtls[eqtls$CD4=='+',]
  tcel.chen.sig <- which(cd4_specific_eqtls$chens.tcell.p <= chen.sig.threshold)
  tcel.concordance <- (sum(sign(cd4_specific_eqtls$zScoreSwapped) == sign(cd4_specific_eqtls$chens.tcell), na.rm = TRUE)) / sum(!is.na(cd4_specific_eqtls$chens.tcell))
  
  ##################
  ## if logModulus is selected. then transform interaction.beta and make the label for the plot. 
  x.label.plot <- paste0(decon.CT.name, "xGenotype")
  x.label.plot <- paste0("log.modulus(", x.label.plot, ")")
  
  
  neut.plot <- allelic.conconcordance.plot(eqtls =neut_specific_eqtls, x.label.plot=x.label.plot, 
                                           concordance = neut.concordance, 
                                           chens.col = "chens.neut", 
                                           chen.sig.rows = neut.chen.sig,
                                           decon.color=deconColors["IT1"], 
                                           y.label.plot = "Neutrophil eQTLs z-score"
  )
  
  
    mono.plot <- allelic.conconcordance.plot(eqtls =mono_specific_eqtls, x.label.plot=x.label.plot,
                                           concordance = mono.concordance, 
                                           chens.col = "chens.mono", 
                                           chen.sig.rows = mono.chen.sig, 
                                           decon.color=deconColors["IT2"], 
                                           y.label.plot = "Monocytes eQTLs z-score")
  
    tcel.plot <- allelic.conconcordance.plot(eqtls =cd4_specific_eqtls, x.label.plot=x.label.plot, 
                                           concordance = tcel.concordance, 
                                           chens.col = "chens.tcell", 
                                           chen.sig.rows = tcel.chen.sig, 
                                           decon.color=deconColors["IT12"], 
                                           y.label.plot = "T cell eQTLs z-score")
  
  
   concordance.plot <- arrangeGrob(neut.plot,mono.plot, tcel.plot , nrow=1)
   
  ## adding back to eQTL the other decon results for further comparisson
  all.decon.colums <- colnames(decon)[grep(colnames(decon), pattern = "IT")]
  decon <- decon[rownames(eqtls),all.decon.colums]
  
  return(list(eqtls=eqtls,
              decon.filtered=decon,
              concordance= c(neut=neut.concordance, mono=mono.concordance, tcel=tcel.concordance),
              plot=concordance.plot))
}




################################
proxy_based_eqtls <- read.table('proxy_celltype_specific_effects.txt', header=T,sep='\t',check.names = F,stringsAsFactors = FALSE)
celltype_effects <- data.frame(celltypes=c('Neutrophil','CD4','CD8','Monocytes','Bcells'),
                               counts=c(table(proxy_based_eqtls$Neutrophils)[['+']],
                                        table(proxy_based_eqtls$CD4)[['+']],
                                        table(proxy_based_eqtls$CD8)[['+']],
                                        table(proxy_based_eqtls$Monocytes)[['+']],
                                        table(proxy_based_eqtls$Bcells)[['+']]))


ggplot(celltype_effects, aes(x=celltypes, y = counts, fill=celltypes))+
  geom_bar(stat='identity')+
  theme_classic()+
  ggtitle(paste0('Number of eQTL per celltype module'))+
  geom_text(aes(label=counts), position=position_dodge(width=0.9), vjust=-0.25)+ 
  scale_fill_manual(values=c("grey10","grey25","grey40","grey55","grey70","grey85"))+
  guides(fill=F)
ggsave('figures/number_of_eqtl_per_celltype_module.png',width=8, height=8)


mono.chen <- read.delim("/Users/NPK/UMCG/projects/deconvolution/chen_filteredStat/bios_freeze2_eqtls_Chens_mono.txt", 
                        stringsAsFactors = FALSE, as.is = TRUE, sep="\t")

neut.chen <- read.delim("/Users/NPK/UMCG/projects/deconvolution/chen_filteredStat/bios_freeze2_eqtls_Chens_neut.txt", 
                        stringsAsFactors = FALSE, as.is = TRUE, sep="\t")

tcel.chen <- read.delim("/Users/NPK/UMCG/projects/deconvolution/chen_filteredStat/bios_freeze2_eqtls_Chens_tcel.txt", 
                        stringsAsFactors = FALSE, as.is = TRUE, sep="\t")

tcel.chen$ref_allele <- sapply(tcel.chen$pos, function(x){unlist(strsplit(x, split="_"))[2]})
tcel.chen$alt_allele <- sapply(tcel.chen$pos, function(x){unlist(strsplit(x, split="_"))[3]})
#tcel.chen <- tcel.chen[-which(is.na(tcel.chen$eQTL)),]
rownames(tcel.chen) <- tcel.chen$eQTL

neut.chen$ref_allele <- sapply(neut.chen$pos, function(x){unlist(strsplit(x, split="_"))[2]})
neut.chen$alt_allele <- sapply(neut.chen$pos, function(x){unlist(strsplit(x, split="_"))[3]})
#3neut.chen <- neut.chen[-which(is.na(neut.chen$eQTL)),]
rownames(neut.chen) <- neut.chen$eQTL

mono.chen$ref_allele <- sapply(mono.chen$pos, function(x){unlist(strsplit(x, split="_"))[2]})
mono.chen$alt_allele <- sapply(mono.chen$pos, function(x){unlist(strsplit(x, split="_"))[3]})
#mono.chen <- mono.chen[-which(is.na(mono.chen$eQTL)),]
rownames(mono.chen) <- mono.chen$eQTL


chen.index <- Reduce(intersect, list(rownames(mono.chen), rownames(neut.chen), rownames(tcel.chen)))

all(mono.chen[chen.index, "ref_allele"] == neut.chen[chen.index, "ref_allele"])
all(mono.chen[chen.index, "ref_allele"] == tcel.chen[chen.index, "ref_allele"])
all(neut.chen[chen.index, "ref_allele"] == tcel.chen[chen.index, "ref_allele"])

all.chen.eQTLs <- unique(c(tcel.chen$eQTL, mono.chen$eQTL, neut.chen$eQTL))
chen.var.info <- data.frame(eqtl_pair = all.chen.eQTLs, row.names = all.chen.eQTLs)
chen.var.info[rownames(tcel.chen),"ref_allele"] <- tcel.chen[,"ref_allele"]
chen.var.info[rownames(mono.chen),"ref_allele"] <- mono.chen[,"ref_allele"]
chen.var.info[rownames(neut.chen),"ref_allele"] <- neut.chen[,"ref_allele"]

chen.var.info[rownames(tcel.chen),"alt_allele"] <- tcel.chen[,"alt_allele"]
chen.var.info[rownames(mono.chen),"alt_allele"] <- mono.chen[,"alt_allele"]
chen.var.info[rownames(neut.chen),"alt_allele"] <- neut.chen[,"alt_allele"]


chens <- list(neut=neut.chen, mono=mono.chen, tcel=tcel.chen)

SNP.info <- read.delim("/Users/NPK/UMCG/projects/deconvolution/chen_filteredStat/CODAM_bios_freeze2_ciseqtls_variantInfo.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE, as.is = TRUE)
row.names(SNP.info) <- SNP.info[,1]
colnames(SNP.info)[1:3] <- c("SNP", "ALT","REF")


itInfo <- read.csv('/Users/NPK/UMCG/projects/deconvolution/pathway/itDecon_extra_colors.csv')

split_path <- function(x) if (dirname(x)==x) x else c(basename(x),split_path(dirname(x)))

logModulus <- function(x){
  signs <- sign(x)
  l <- log(abs(x)+1)
  l <- l*signs
  return(l)
}


deconColors <- c(IT1="#E41A1C",  IT2="#377EB8", IT12="#4DAF4A",IT11="#984EA3",  IT13="#FF7F00", IT8="#FFFF33")
  

concordance.gran <- chens.decon.concordance.proxyEqtls(eqtls = proxy_based_eqtls, 
                                              chens = chens, 
                                              decon.CT="IT1", 
                                              decon.CT.name = "Granulocytes")

  
  
pdf(paste0('figures/conc_proxy_eqtl.pdf'), width=7, height=4)
grid.arrange(concordance.gran$plot)
dev.off()
