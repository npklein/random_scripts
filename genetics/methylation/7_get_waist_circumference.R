#######################################################
# Script to parse LLDeep pheno data for BIOS freeze 2 #
#######################################################

library(memisc)
library(data.table)
library(dplyr)


setwd("/groups/umcg-lld/tmp03/umcg-aclaringbould/phenotypes/DAG1_DAG2_Phenodata_Release_2015_01_extended")

savFiles <- dir(path = ".", pattern='\\.sav$', full.names = TRUE)

datas <- list()
descriptions <- list()

## READ IN DATA ##

for(file in savFiles){
  
  print("-----------------------------------------------", quote = FALSE)
  print(paste0("File: ", file), quote = FALSE)
  
  spssData <- spss.system.file(file)
  data <- as.data.frame(as.data.set(spssData))
  
  print(paste0("Samples: ", nrow(data), " Columns: ", ncol(data)), quote = FALSE)
  
 
  if(nrow(data) > 13395){
    warning("Found ", nrow(data), " samples in ", file)
    next
  }
  if(is.null(data$pseudoidext)){
    stop("No pseudoidext in ", file)
  }
  if(nrow(data) != length(unique(data$pseudoidext))){
    warning("Duplicate samples in ", file)
    next
  }
  
  data$pseudoidext <- as.character(data$pseudoidext)
  
  drops <- c("project", "subjectid","encountercode")
#  print(colnames(data))
#  print(!(colnames(data) %in% drops))
  data <- data[,!(colnames(data) %in% drops)]
  datas[[file]] <- data
  
  description <- data.frame(column = names(description(spssData)), description = as.vector(unlist(description(spssData))), stringsAsFactors = FALSE)
  description <- description[!(description$column %in% drops),]
  descriptions[[file]] <- description
  
  write.table(data, file = gsub(".sav", ".txt", gsub('./', '../parsed3/', file)), quote = FALSE, sep = "\t", col.names = NA)
}

hist(sapply(datas, nrow))
hist(sapply(datas, ncol))

# str(datas[[1]])
# str(descriptions)

mergeData <- function(x,y){
  return(merge(x,y,by="pseudoidext", all=TRUE))
}

#merge datasets
hist(sapply(datas, function(x){sum(colnames(x) %in% "pseudoidext")}))

sapply(datas, function(x){sum(colnames(x) %in% "pseudoidext")})

dataMerged <- Reduce(mergeData, datas)
# str(dataMerged)
# dim(dataMerged)

allIds <- Reduce(c, sapply(datas, function(x){return(x["pseudoidext"])}))
#length(unique(allIds))

#specifically read in medication file
spssData <- spss.system.file("./Questionnaire - Health - Medication - Prescription (Que_Hea_Med_Pre1).sav")
medicationData <- as.data.frame(as.data.set(spssData))
medicationData$pseudoidext <- as.character(medicationData$pseudoidext)


## ADJUST DATA ##

dataMerged$Biobank <- "LL"
#dataMerged$pheno_id <- ""
dataMerged$Ascertainment_criterion <- "Unselected population based sample"
dataMerged$GWAS_Chip <- "Immunnochip, Cytochip"
#dataMerged$GWAS_DataGeneration_Date <- dataMerged$bl1year
#dataMerged$DNA_BloodSampling_Age <- dataMerged$age_bl1
#dataMerged$DNA_BloodSampling_Date <- dataMerged$bl1year
dataMerged$DNA_BloodSampling_Time <- "Morning" #information from Ettje and Jackie 
dataMerged$DNA_Source <- "Whole peripheral blood"
dataMerged$DNA_Extraction_Method <- NA
#dataMerged$DNA_Extraction_Date <- dataMerged$bl1year
dataMerged$DNA_QuantificationMethod <- NA
dataMerged$DNA_A260A280ratio <- NA
dataMerged$RNA_BloodSampling_Age <- 2013-dataMerged$birthyear
dataMerged$RNA_Sampling_Date <- "dec-2013" #information from Ettje and Jackie
dataMerged$RNA_Sampling_Time <- NA
dataMerged$RNA_Source <- "Whole peripheral blood"
dataMerged$RNA_Extraction_Date <- "dec-2013"
dataMerged$RNA_Extraction_Method <- "PAXGene"
dataMerged$RNA_RIN <- NA
dataMerged$RNA_A260280ratio <- NA
#dataMerged$BirthYear <- dataMerged$birthyear
dataMerged$Sex <- recode(dataMerged$geslacht, 'Male' = '0', 'Female' = '1')
dataMerged$Smoking <- NA
dataMerged$Smoking[ (!is.na(dataMerged$smk1) & dataMerged$smk1 == "No") ] <- 0 #never smoker
dataMerged$Smoking[ (!is.na(dataMerged$smk3) & dataMerged$smk3 == "Yes") ] <- 2 #current smoker
#dataMerged$Smoking[ (!is.na(dataMerged$smk1) & dataMerged$smk1 == "Yes") & (!is.na(dataMerged$smk3) & dataMerged$smk3 == "No") ] <- 1 #former smoker ('Have you ever smoked for a full year? = yes', 'Do you smoke now, or have you smoked in the past month? = no')
dataMerged$Smoking[ (!is.na(dataMerged$smk1) & dataMerged$smk1 == "Yes") & (!is.na(dataMerged$smk5) & dataMerged$smk5 == "Yes") ] <- 1 #former smoker ('Have you ever smoked for a full year? = yes', 'Have you stopped smoking? = yes')
dataMerged$Smoking_Age <- ifelse(dataMerged$Smoking == 2, dataMerged$age_bl1, dataMerged$smk6) #if current smoker, enter current age, otherwise enter answer to question 'How old were you when you stopped smoking?'
#dataMerged$Lipids_BloodSampling_Age <- dataMerged$age_bl1
#dataMerged$Lipids_BloodSampling_Date <- dataMerged$bl1year
dataMerged$Lipids_BloodSampling_Time <- "Morning" #information from Ettje and Jackie
dataMerged$Lipids_BloodSampling_Fasting <- recode(dataMerged$nuchter, 'Yes' = '1', 'No' = '0')
#dataMerged$TotChol <- dataMerged$cho
#dataMerged$HDLchol <- dataMerged$hdc
#dataMerged$Triglycerides <- dataMerged$tgl
#dataMerged$LDLchol <- dataMerged$ldc
dataMerged$LDLcholMethod <- "2"
#dataMerged$LipidsMed_Age <- dataMerged$health39
dataMerged$LipidMed <- 0
dataMerged$LipidMed[dataMerged$pseudoidext %in% medicationData$pseudoidext[grep("^C10", medicationData$atccode)]] <- 2
dataMerged$LipidMed[dataMerged$pseudoidext %in% medicationData$pseudoidext[grep("^C10AA", medicationData$atccode)]] <- 1
dataMerged$LipidMed[dataMerged$pseudoidext %in% medicationData$pseudoidext[grep("^C10B", medicationData$atccode)]] <- 1
#dataMerged$Anthropometry_Age <- dataMerged$age_bl1
dataMerged$Height <- NA
dataMerged$Height[ dataMerged$meetstand == "Standing" ] <- dataMerged$lengte #NA if dataMerged$meetstand = "Sitting"
#dataMerged$Weight <- dataMerged$gewicht
#dataMerged$CRP_BloodSampling_Age <- dataMerged$age_bl1
#dataMerged$CRP_BloodSampling_Date <- dataMerged$bl1year
dataMerged$CRP_BloodSampling_Time <- NA
#dataMerged$hsCRP <- dataMerged$lcrp
#dataMerged$CellCount_BloodSampling_Age <- dataMerged$age_bl1
#dataMerged$CellCount_BloodSampling_Date <- dataMerged$bl1year
dataMerged$CellCount_BloodSampling_Time <- NA
#dataMerged$WBC <- dataMerged$leu
#dataMerged$RBC <- dataMerged$er
#dataMerged$HGB <- dataMerged$hb
#dataMerged$HCT <- dataMerged$ht
dataMerged$MCV <- NA
dataMerged$MCH <- NA
dataMerged$MCHC <- NA
dataMerged$CHCM <- NA
dataMerged$CH <- NA
dataMerged$RDW <- NA
dataMerged$HDW <- NA
#dataMerged$PLT <- dataMerged$tr
dataMerged$MPV <- NA
#dataMerged$Neut <- dataMerged$gr
#dataMerged$Lymph <- dataMerged$ly
#dataMerged$Mono <- dataMerged$mo
#dataMerged$Eos <- dataMerged$eo
#dataMerged$Baso <- dataMerged$ba
dataMerged$LUC <- NA
#dataMerged$Neut_Perc <- dataMerged$grp
#dataMerged$Lymph_Perc <- dataMerged$lyp
#dataMerged$Mono_Perc <- dataMerged$mop
#dataMerged$Eos_Perc <- dataMerged$eop
#dataMerged$Baso_Perc <- dataMerged$bap
dataMerged$LUC_Perc <- NA


## MAP TO LLDEEP IDs ##

mapping1 <- fread("/groups/umcg-lld/tmp03/umcg-aclaringbould/phenotypes/LLDeep_all_phenos/OV11_0098_DAG1_Linking_pseudoidext_external_20150518.txt")
mapping2 <- fread("/groups/umcg-lld/tmp03/umcg-aclaringbould/phenotypes/LLDeep_all_phenos/LLDeep_samples_V09.txt")
mapping <- merge(mapping1, mapping2, by="PSEUDOID")


pheno <- merge(dataMerged, mapping, by.x = "pseudoidext", by.y = "PSEUDOIDEXT")
#dim(dataMerged)
#dim(pheno)

## SELECT COLUMNS OF INTEREST ##

pheno <- select(pheno,
                pheno_id = LLDeep_SampleID, 
                waist_circumference = taille,
                Smoking)
                

## WRITE DATA ##

write.table(pheno, "/groups/umcg-lld/tmp03/umcg-ndeklein/EWAS/EWAS_scripts/LLD_waistCircumference.txt", sep = "\t", quote = FALSE, row.names = FALSE)
