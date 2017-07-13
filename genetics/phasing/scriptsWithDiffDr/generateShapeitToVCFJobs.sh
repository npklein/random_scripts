for dr in $(seq 0 0.1 1);
do
  echo ${dr}
  for i in {1..22}
  do

    jobDir=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/jobs/phasingDifferentFilters/shapeitToVCF/DR${dr}/chr${i}
    mkdir -p $jobDir
    INPUTDIR=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/results/beagleDifferentFilters/shapeit_DR${dr}/chr${i}/
    OUTPUTDIRGENECHUNKS=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/results/beagleDifferentFilters/shapeitVCF_DR${dr}/
    mkdir -p $OUTPUTDIRGENECHUNKS
    RESULTSDIR=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/results/beagleDifferentFilters/shapeitVCF_DR${dr}/chr${i}/
    phasedGeneChunks=${OUTPUTDIRGENECHUNKS}/phaseGeneChunks_chr${i}.DR${dr}.csv
    echo "CHR,chromosomeChunk" > $phasedGeneChunks
    for chunk in $(find $INPUTDIR/*.hap.gz -not -empty -ls | awk '{print $11}');
    do
      CHR=$(echo $chunk | awk -F"." '{print $3}' | awk -F"chr" '{print $2}')
      START=$(echo $chunk | awk -F"." '{print $4}')
      END=$(echo $chunk | awk -F"." '{print $5}')
      echo "$CHR,$CHR:$START-$END" >>  $phasedGeneChunks
    done
    echo "wrote to $phasedGeneChunks"
    while read line
    do
      CHR=`echo $line | awk '{print $1}' FS=","`
      START=`echo $line | awk '{print $2}' FS=":" | awk '{print $1}' FS="-"`
      END=`echo $line | awk '{print $2}' FS=":" | awk '{print $2}' FS="-"`
      CHUNK=$CHR.$START.$END.DR${dr}
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
done
done
