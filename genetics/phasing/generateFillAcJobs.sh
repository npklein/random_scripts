module load tabix
phasedGeneChunks=../geneChunks.20170519.csv
INITIALVCFDIR=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged/results/filterGQ20_callRate50_noRnaEditSites/
while read line
do
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
    echo $line
    jobDir=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/jobs/phasing/fillAC/chr${CHR}/
    mkdir -p $jobDir
    INPUTDIR=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/results/phasing/shapeitVCF/chr${CHR}/
    RESULTSDIR=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/results/phasing/shapeitVCFwithAC/chr${CHR}/
    CHUNK=$CHR.$START.$END
      echo "#!/bin/bash
#SBATCH --job-name=chr$CHUNK.fillAC
#SBATCH --output=fillAC.chr$CHUNK.out
#SBATCH --error=fillAC.chr$CHUNK.err
#SBATCH --time=5:59:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 8gb
#SBATCH --nodes 1
#SBATCH --export=NONE
#SBATCH --get-user-env=30L
#SBATCH --qos=regular

set -e
set -u

ENVIRONMENT_DIR='.'


module load tabix/0.2.6-foss-2015b
module load BCFtools/1.5-foss-2015b
module list

echo \"## \"\$(date)\" Start \$0\"

mkdir -p ${RESULTSDIR}

bcftools +fill-tags ${INPUTDIR}/genotypes_BIOS_LLDeep_noRNAedit_merged.chr$CHUNK.shapeit.phased.vcf.gz \\
          -o $RESULTSDIR/genotypes_BIOS_LLDeep_noRNAedit_merged.chr$CHUNK.shapeit.phased.withAC.vcf.gz \\
          -O z

tabix -f $RESULTSDIR/genotypes_BIOS_LLDeep_noRNAedit_merged.chr$CHUNK.shapeit.phased.withAC.vcf.gz

echo \"returncode: \$?\";
cd ${RESULTSDIR}
name=\$(basename ${RESULTSDIR}/genotypes_BIOS_LLDeep_noRNAedit_merged.chr$CHUNK.shapeit.phased.withAC.vcf.gz)
md5sum $name > $name.md5

cd -
echo \"succes moving files\";


echo \"## \"\$(date)\" ##  \$0 Done \"

">$jobDir/fillAC.chr$CHUNK.sh

done<$phasedGeneChunks
