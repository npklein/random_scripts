resultsdir=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/results/
RESULTSDIR=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/results/genotypeVcfGL_filtered_tmp

mkdir -p ${RESULTSDIR}

for i in {1..22}
do

echo "#!/bin/bash
#SBATCH --job-name=chr${i}_BIOS_freeze2.1_LLDeep_Diagnostics_FilterGenotypeVcfGL
#SBATCH --output=chr${i}_BIOS_freeze2.1_LLDeep_Diagnostics_FilterGenotypeVcfGL.out
#SBATCH --error=chr${i}_BIOS_freeze2.1_LLDeep_Diagnostics_FilterGenotypeVcfGL.err
#SBATCH --time=5:59:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 10gb
#SBATCH --export=NONE
#SBATCH --get-user-env=30L
#SBATCH --qos=dev

set -e
set -u

ENVIRONMENT_DIR='.'


module load Python
module load tabix
module list

echo \"## \"\$(date)\" Start \$0\"

zcat $resultsdir/beagle_vcf_filtered/genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr$i.beagle.genotype.probs.gg.vcf.gz | \\
            awk '{print \$1 \":\" \$2 \"-\" \$2 }' > $resultsdir/beagle_vcf_filtered/all_positions_chr${i}.intervals

python /groups/umcg-bios/tmp03/projects/phasing/selectVariants.py \\
    $resultsdir/genotypeVcfGL//genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr${i}.genotypeGVCF.gg.vcf.gz \\
    $resultsdir/beagle_vcf_filtered/all_positions_chr${i}.intervals \\
    $RESULTSDIR/genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr${i}.genotypeGVCF.DR2Filtered.gg.vcf

cd $RESULTSDIR
bgzip ${RESULTSDIR}/genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr${i}.genotypeGVCF.DR2Filtered.gg.vcf
tabix ${RESULTSDIR}/genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr${i}.genotypeGVCF.DR2Filtered.gg.vcf.gz
md5sum ${RESULTSDIR}/genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr${i}.genotypeGVCF.DR2Filtered.gg.vcf.gz > \\
        ${RESULTSDIR}/genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr${i}.genotypeGVCF.DR2Filtered.gg.vcf.gz.md5
md5sum ${RESULTSDIR}/genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr${i}.genotypeGVCF.DR2Filtered.gg.vcf.gz.tbi > \\
        ${RESULTSDIR}/genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr${i}.genotypeGVCF.DR2Filtered.gg.vcf.gz.tbi.md5
        
echo \"## \"\$(date)\" Done \$0\"

">./jobs/filter_GLVCF_tmp/chr${i}_BIOS_freeze2.1_LLDeep_Diagnostics_FilterGenotypeVcfGL.sh
done
