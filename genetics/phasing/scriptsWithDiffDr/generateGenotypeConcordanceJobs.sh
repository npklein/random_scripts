concordanceFolder=/groups/umcg-bios/tmp03/projects/genotypeConcordance/
wgs_folder=${concordanceFolder}/GoNL_WGS_calls/
GQ=20
callrate=0

for dr in $(seq 0 0.1 0.9);
do
  echo ${dr}
  for CHR in {1..22}
  do

    jobDir=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/jobs/phasingDifferentFilters/compareGenotypes/DR${dr}/
    mkdir -p $jobDir
    INPUTDIR=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/results/beagleDifferentFilters/shapeitVCF_DR${dr}/
    shapeitVCF=${INPUTDIR}/genotypes_BIOS_LLDeep_Diagnostics_merged.chr${CHR}.DR${dr}.concat.shapeit.phased.vcf.gz
    wgsVCF=${wgs_folder}/gonl-abc_samples.chr${CHR}.release5.NoMAFSelection.vcf.gz

    RESULTDIR=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/results/beagleDifferentFilters/genotypeComparison/DR${dr}/
    echo "#!/bin/bash
#SBATCH --job-name=chr$CHR.DR${dr}.CompareGenotypes
#SBATCH --output=chr$CHR.DR${dr}.CompareGenotypes.out
#SBATCH --error=chr$CHR.DR${dr}.CompareGenotypes.err
#SBATCH --time=5:59:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 20gb
#SBATCH --nodes 1
#SBATCH --export=NONE
#SBATCH --get-user-env=30L
#SBATCH --qos=leftover

set -e
set -u

ENVIRONMENT_DIR='.'


module load tabix/0.2.6-foss-2015b
module load shapeit/v2.r837-static
module list

set -e
module load tabix/0.2.6-foss-2015b
module load VCFtools/0.1.14-foss-2015b-Perl-5.22.0-bare 


# INCASE IT IS RERUN, REMOVE PREVIOUS RESULT AS NEW DATA GETS APPENDED TO IT
#	merged calling
rm -f ${RESULTDIR}/RNA_GQ${GQ}_call${callrate}_count.txt
rm -f ${RESULTDIR}/skipped*
rm -f ${RESULTDIR}/WGS_filtered_count.txt

mkdir -p ${RESULTDIR}

#	Step 2: Compare
java -jar ${concordanceFolder}/CompareGenotypeCalls-1.5-SNAPSHOT/CompareGenotypeCalls.jar \\
	-d1 ${wgsVCF} \\
	-D1 VCF \\
	-d2 ${shapeitVCF} \\
	-D2 VCF \\
	-o ${RESULTDIR}/intersected_filter_GQ${GQ}_call${callrate}_chr${CHR}_DR${dr} \\
	-s ${concordanceFolder}/${rnaseq_rare_variants}linking_file.txt 

#awk -v x=\${CHR} '/Number of SNPs compared:/{print x, "compared", \$NF} /Skipped vars due to incosistant alleles:/{print x, "skipped", \$NF}' ${projectDir}${output_folder}intersected_GQ${GQ}_call${callrate}_chr${CHR}.log  >> ${projectDir}${output_folder}skipped_merged_no_filter.txt
#" > ${jobDir}/chr${CHR}.DR${dr}.CompareGenotypes.sh

done
done
