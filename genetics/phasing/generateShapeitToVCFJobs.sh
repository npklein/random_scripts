module load tabix
phasedGeneChunks=../geneChunks.20170519.csv
INITIALVCFDIR=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged/results/filterGQ20_callRate50_noRnaEditSites/
while read line
do
    #echo $line
    CHR=`echo $line | awk '{print $1}' FS=","`
    START=`echo $line | awk '{print $2}' FS=":" | awk '{print $1}' FS="-"`
    END=`echo $line | awk '{print $2}' FS=":" | awk '{print $2}' FS="-"`
    numberOfSnps=$(tabix $INITIALVCFDIR/genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr$CHR.filter.GQ20_callRate50.BiallelicSNVsOnly.noRnaEditSites.gg.vcf.gz $CHR:${START}-${END} | wc -l)
    if [ $numberOfSnps -eq 0 ];
    then
      echo "Chunk $CHR:${START}-${END} did not have any SNPs, skipping"
      continue
    fi
    jobDir=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/jobs/phasing/shapeitToVCF/chr${CHR}/
    mkdir -p $jobDir
    INPUTDIR=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/results/phasing/shapeit/chr${CHR}/
    RESULTSDIR=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/results/phasing/shapeitVCF/chr${CHR}/
    phasedGeneChunks=../beagledGeneChunks.13092017.csv
    CHUNK=$CHR.$START.$END
      echo "#!/bin/bash
#SBATCH --job-name=chr$CHUNK.ConvertShapeitToVCF
#SBATCH --output=ConvertShapeitToVCF.chr$CHUNK.out
#SBATCH --error=ConvertShapeitToVCF.chr$CHUNK.err
#SBATCH --time=0:30:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 8gb
#SBATCH --nodes 1
#SBATCH --export=NONE
#SBATCH --get-user-env=30L
#SBATCH --qos=dev

set -e
set -u

ENVIRONMENT_DIR='.'


module load tabix/0.2.6-foss-2015b
module load shapeit/v2.r837-static
module list

echo \"## \"\$(date)\" Start \$0\"

mkdir -p ${RESULTSDIR}

zcat ${INPUTDIR}/genotypes_BIOS_LLDeep_Diagnostics_merged.chr$CHUNK.shapeit.phased.hap.gz \\
        > ${INPUTDIR}/genotypes_BIOS_LLDeep_Diagnostics_merged.chr$CHUNK.shapeit.phased.haps

if [ -f  ${INPUTDIR}/genotypes_BIOS_LLDeep_Diagnostics_merged.chr$CHUNK.shapeit.phased.hap.gz.sample ];
then
    mv ${INPUTDIR}/genotypes_BIOS_LLDeep_Diagnostics_merged.chr$CHUNK.shapeit.phased.hap.gz.sample \\
        ${INPUTDIR}/genotypes_BIOS_LLDeep_Diagnostics_merged.chr$CHUNK.shapeit.phased.sample
fi

shapeit \\
 -convert \\
 --input-haps ${INPUTDIR}/genotypes_BIOS_LLDeep_Diagnostics_merged.chr$CHUNK.shapeit.phased \\
 --output-vcf ${RESULTSDIR}/genotypes_BIOS_LLDeep_Diagnostics_merged.chr$CHUNK.shapeit.phased.vcf.gz

 echo \"returncode: \$?\";
 cd ${RESULTSDIR}
 bname=\$(basename ${RESULTSDIR}/genotypes_BIOS_LLDeep_Diagnostics_merged.chr$CHUNK.shapeit.phased.vcf.gz)

 # has to be bgzipped
 gunzip \${bname}
 bgzip \${bname%.gz}
 tabix \${bname}
 echo \"making md5sums\"
 md5sum \${bname} > \${bname}.md5

 cd -
 echo \"succes moving files\";


echo \"## \"\$(date)\" ##  \$0 Done \"

">$jobDir/ConvertShapeitToVCF.chr$CHUNK.sh

done<$phasedGeneChunks
