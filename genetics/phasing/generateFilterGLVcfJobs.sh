module load tabix
glDir=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/results/prePhasing/PLtoGL
INITIALVCFDIR=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged/results/filterGQ20_callRate50_noRnaEditSites/
while read line
do

  echo $line
  if [ "$line" = "CHR,chromosomeChunk" ];
  then
    continue
  fi
  CHR=`echo $line | awk '{print $1}' FS=","`
  START=`echo $line | awk '{print $2}' FS=":" | awk '{print $1}' FS="-"`
  END=`echo $line | awk '{print $2}' FS=":" | awk '{print $2}' FS="-"`
  numberOfSnps=$(tabix $INITIALVCFDIR/genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr$CHR.filter.GQ20_callRate50.BiallelicSNVsOnly.noRnaEditSites.gg.vcf.gz $CHR:${START}-${END} | wc -l)
  if [ $numberOfSnps -eq 0 ];
  then
    echo "Chunk $CHR:${START}-${END} did not have any SNPs, skipping"
    continue
  fi
  jobsDir=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/jobs/prePhasing/filterGL/chr${CHR}
  mkdir -p $jobsDir
  INPUTDIR=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/results/beagleChunksNoRNAedit/chr${CHR}
  INTERVALDIR=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/results/prePhasing/intervals/chr${CHR}
  mkdir -p $INTERVALDIR
  RESULTSDIR=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/results/prePhasing/genotypeVcfGL_filtered/chr${CHR}
  mkdir -p ${RESULTSDIR}
  vcfPrefix=genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr
  vcfPostfix=.filter.GQ20_callRate50.BiallelicSNVsOnly.noRnaEditSites.probs.gg.vcf.gz

echo chr$CHR
echo "#!/bin/bash
#SBATCH --job-name=chr$CHR.$START.${END}_BIOS_freeze2.1_LLDeep_noRnaEditSites_FilterGenotypeVcfGL
#SBATCH --output=chr$CHR.$START.${END}_BIOS_freeze2.1_LLDeep_noRnaEditSites_FilterGenotypeVcfGL.out
#SBATCH --error=chr$CHR.$START.${END}_BIOS_freeze2.1_LLDeep_noRnaEditSites_FilterGenotypeVcfGL.err
#SBATCH --time=5:59:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 10gb
#SBATCH --export=NONE
#SBATCH --get-user-env=30L
#SBATCH --qos=regular

set -e
set -u

ENVIRONMENT_DIR='.'

module load Python/3.5.1-foss-2015b
module load tabix
module list

echo \"## \"\$(date)\" Start \$0\"

if [ ! -f $INTERVALDIR/all_positions_chr$CHR.$START.${END}.intervals ];
then
  if [ ! -f $INPUTDIR/${vcfPrefix}$CHR.$START.${END}${vcfPostfix} ];
  then
    echo "$INPUTDIR/${vcfPrefix}$CHR.$START.${END}${vcfPostfix} does not exist"
    exit 1;
  fi
  zcat $INPUTDIR/${vcfPrefix}$CHR.$START.${END}${vcfPostfix} | \\
            grep -v '^#' | awk '{print \$1 \":\" \$2 \"-\" \$2 }' > $INTERVALDIR/all_positions_chr$CHR.$START.${END}.intervals
else
  echo "$INTERVALDIR/all_positions_chr$CHR.$START.${END}.intervals already exists, using that"
fi

if [ ! -f $RESULTSDIR/genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr$CHR.$START.${END}.genotypeGVCF.gg.vcf ];
then
  python /groups/umcg-bios/tmp03/projects/phasing/selectVariants.py \\
      $glDir//genotypes_BIOSfreeze2.1_LLDeep_noRnaEditSites_chr${CHR}.vcf.gz \\
      $INTERVALDIR/all_positions_chr$CHR.$START.${END}.intervals \\
      $RESULTSDIR/genotypes_BIOSfreeze2.1_LLDeep_noRnaEditSites_merged.chr$CHR.$START.${END}.genotypeGVCF.gg.vcf
else
  echo "$RESULTSDIR/genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr$CHR.$START.${END}.genotypeGVCF.gg.vcf already exists, using that"
fi

cd $RESULTSDIR
echo "bgzipping..."
bgzip -f $RESULTSDIR/genotypes_BIOSfreeze2.1_LLDeep_noRnaEditSites_merged.chr$CHR.$START.${END}.genotypeGVCF.gg.vcf
echo "tabix..."
tabix -f $RESULTSDIR/genotypes_BIOSfreeze2.1_LLDeep_noRnaEditSites_merged.chr$CHR.$START.${END}.genotypeGVCF.gg.vcf.gz
echo "md5sum..."
md5sum $RESULTSDIR/genotypes_BIOSfreeze2.1_LLDeep_noRnaEditSites_merged.chr$CHR.$START.${END}.genotypeGVCF.gg.vcf.gz > \
      $RESULTSDIR/genotypes_BIOSfreeze2.1_LLDeep_noRnaEditSites_merged.chr$CHR.$START.${END}.genotypeGVCF.gg.vcf.gz.md5
md5sum $RESULTSDIR/genotypes_BIOSfreeze2.1_LLDeep_noRnaEditSites_merged.chr$CHR.$START.${END}.genotypeGVCF.gg.vcf.gz.tbi > \
      $RESULTSDIR/genotypes_BIOSfreeze2.1_LLDeep_noRnaEditSites_merged.chr$CHR.$START.${END}.genotypeGVCF.gg.vcf.gz.tbi.md5

echo \"## \"\$(date)\" Done \$0\"

">$jobsDir/chr$CHR.$START.${END}_BIOS_freeze2.1_LLDeep_noRnaEditSites_FilterGenotypeVcfGL.sh
done<../geneChunks.20170519.csv
