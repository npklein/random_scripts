jobDir=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/jobs/phasingDifferentFilters/selectGonlFilteredGL/
mkdir -p $jobDir

for i in {1..22}
do

for dr in $(seq 0 0.1 1); do
INPUTDIR=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/results/beagleDifferentFilters/genotypeVcfGL_filtered_DR${dr}
OUTPUTDIR=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/results/beagleDifferentFilters/genotypeVcfGL_filtered_DR${dr}_GoNL
mkdir -p $OUTPUTDIR
INPUTVCF=genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr${i}.genotypeGVCF.DR${dr}.gg.vcf.gz
OUTPUTVCF=genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr${i}.genotypeGVCF.GoNL.DR${dr}.gg.vcf.gz
echo "#!/bin/bash
#SBATCH --job-name=chr${i}_BIOS_LLDeep_Diagnostics_freeze2_SelectGonl_DR${dr}
#SBATCH --output=chr${i}_BIOS_LLDeep_Diagnostics_freeze2_SelectGonl_DR${dr}.out
#SBATCH --error=chr${i}_BIOS_LLDeep_Diagnostics_freeze2_SelectGonl_DR${dr}.err
#SBATCH --time=5:59:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 10gb
#SBATCH --export=NONE
#SBATCH --get-user-env=30L
#SBATCH --qos=leftover

set -e
set -u

ENVIRONMENT_DIR='.'


module load BCFtools
module list

echo \"## \"\$(date)\" Start \$0\"

bcftools view -S /groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/GoNL_variants/gonlSamplesRNAids.txt \\
             -o ${OUTPUTDIR}/${OUTPUTVCF} \\
             -O z \\
             --force-samples \\
              ${INPUTDIR}/${INPUTVCF} 

">$jobDir/chr${i}_BIOS_LLDeep_Diagnostics_freeze2_SelectGonl_DR${dr}.sh
done
done
