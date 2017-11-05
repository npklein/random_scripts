INITIALVCFDIR=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged/results/filterGQ20_callRate50_noRnaEditSites/
module load tabix
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
  jobDir=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/jobs/prePhasing/convertBeagleToShapeit/chr$CHR/
  mkdir -p $jobDir

  INPUTDIR=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/results/beagleChunksNoRNAedit/chr$CHR
  GLDIR=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/results/prePhasing/genotypeVcfGL_filtered/chr$CHR/
  OUTPUTDIR=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/results/prePhasing/beagleToShapeit/chr$CHR/
  mkdir -p $OUTPUTDIR
  LIKELIHOODVCF=genotypes_BIOSfreeze2.1_LLDeep_noRnaEditSites_merged.chr${CHR}.${START}.${END}.genotypeGVCF.gg.vcf.gz
  BEAGLEVCF=genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr${CHR}.${START}.${END}.filter.GQ20_callRate50.BiallelicSNVsOnly.noRnaEditSites.probs.gg.vcf.gz
echo "#!/bin/bash
#SBATCH --job-name=chr$CHR.$START.${END}_BIOS_LLDeep_noRNAedit_freeze2_ConvertBeagle
#SBATCH --output=chr$CHR.$START.${END}_BIOS_LLDeep_noRNAedit_freeze2_ConvertBeagle.out
#SBATCH --error=chr$CHR.$START.${END}_BIOS_LLDeep_noRNAedit_freeze2_ConvertBeagle.err
#SBATCH --time=5:59:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 10gb
#SBATCH --export=NONE
#SBATCH --get-user-env=30L
#SBATCH --qos=regular
mkdir -p ${OUTPUTDIR}

set -e
set -u

ENVIRONMENT_DIR='.'


module load GATK
module list

echo \"## \"\$(date)\" Start \$0\"


module load tabix
module load prepareGenFromBeagle4
# Glib is also set as dependency of prepareGenFromBeagle4 but still needs to be loaded after
module load GLib
module list

# had to be gzipped again otherwise it will not read the last line
zcat $INPUTDIR/$BEAGLEVCF \\
    > \$TMPDIR/${BEAGLEVCF%.gz}
gzip \$TMPDIR/${BEAGLEVCF%.gz}
zcat $GLDIR/$LIKELIHOODVCF \\
   > \$TMPDIR/${LIKELIHOODVCF%.gz}
gzip \$TMPDIR/${LIKELIHOODVCF%.gz}

#Run conversion script beagle vcf to .hap.gz format
prepareGenFromBeagle4 \
 --likelihoods \$TMPDIR/${LIKELIHOODVCF} \\
 --posteriors \$TMPDIR/${BEAGLEVCF} \\
 --output $OUTPUTDIR/genotypes_BIOSfreeze2.1_LLDeep_noRNAedit_merged.chr$CHR.$START.$END.beagle.Filteredgenotype

echo "returncode: $?";
# these output files are NOT gzipped, so rename them to filename without gz

echo "succes moving files";


">$jobDir/chr$CHR.$START.$END\_BIOS_LLDeep_noRNAedit_freeze2_ConvertBeagle.sh
done<../geneChunks.20170519.csv
