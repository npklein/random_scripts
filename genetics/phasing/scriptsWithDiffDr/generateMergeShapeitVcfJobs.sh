VCFPREFIX=genotypes_BIOS_LLDeep_Diagnostics_merged.chr
VCFPOSTFIX=.shapeit.phased.vcf.gz
for dr in $(seq 0 0.1 1);
do
  echo ${dr}
  for CHR in {1..22}
  do
    INPUTGENECHUNKS=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/results/beagleDifferentFilters/shapeitVCF_DR${dr}/phaseGeneChunks_chr${CHR}.DR${dr}.csv
    jobDir=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/jobs/phasingDifferentFilters/mergeShapeitVCF/DR${dr}/
    mkdir -p $jobDir
    INPUTDIR=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/results/beagleDifferentFilters/shapeitVCF_DR${dr}/chr${CHR}/
    OUTPUTDIR=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/results/beagleDifferentFilters/shapeitVCF_DR${dr}/
    OUTPUTVCF=${OUTPUTDIR}/${VCFPREFIX}${CHR}.DR${dr}.concat${VCFPOSTFIX%.gz}
    echo "#!/bin/bash
#SBATCH --job-name=MergeVCF.chr${CHR}
#SBATCH --output=MergeVCF.chr${CHR}.out
#SBATCH --error=MergeVCF.chr${CHR}.err
#SBATCH --time=5:59:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 8gb
#SBATCH --nodes 1
#SBATCH --export=NONE
#SBATCH --get-user-env=30L
#SBATCH --qos=dev

set -e
set -u

ENVIRONMENT_DIR='.'


module load BCFtools
set -e

# skip header and then parse the chunks
MERGEVCFINPUT=()
echo \"Looping over the chunks...\"

while read line
do
  echo \$line

  if [ "\$line" = \"CHR,chromosomeChunk\" ];
  then
    echo "skipping header"
    continue
  fi

  CHR=\$(echo \$line | awk '{print \$1}' FS=\",\")
  START=\$(echo \$line | awk '{print \$2}' FS=\":\" | awk '{print \$1}' FS=\"-\")
  END=\$(echo \$line | awk '{print \$2}' FS=\":\" | awk '{print \$2}' FS=\"-\")
  CHUNK=\$CHR.\$START.\$END.DR${dr}

  if [ -f ${INPUTDIR}/${VCFPREFIX}\$CHUNK${VCFPOSTFIX} ];
  then
    MERGEVCFINPUT+=(\"${INPUTDIR}/${VCFPREFIX}\$CHUNK${VCFPOSTFIX}\")
  else
    echo ${INPUTDIR}/${VCFPREFIX}\$CHUNK${VCFPOSTFIX} \" does not exist, exitting\"
    exit 1;
  fi
done < ${INPUTGENECHUNKS}

bcftools concat \${MERGEVCFINPUT[@]} > ${OUTPUTVCF}

grep '^#' ${OUTPUTVCF} > ${OUTPUTVCF%.vcf}.sorted.vcf \\
         && grep -v '^#' ${OUTPUTVCF} \\
                        | LC_ALL=C sort -t\$'\t' -k1,1 -k2,2n >> ${OUTPUTVCF%.vcf}.sorted.vcf
mv ${OUTPUTVCF%.vcf}.sorted.vcf ${OUTPUTVCF}

bgzip ${OUTPUTVCF}
tabix ${OUTPUTVCF}.gz


 echo \"returncode: \$?\";
 echo \"succes moving files\";


echo \"## \"\$(date)\" ##  \$0 Done \"

">$jobDir/MergeVCF.chr${CHR}.sh
done
done
