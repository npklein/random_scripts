jobsDir=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/jobs/phasingDifferentFilters/filterGL
mkdir -p $jobsDir

for dr in $(seq 0 0.1 1); 
do
resultsdir=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/results/beagleDifferentFilters/beagle_vcf_filtered_DR${dr}
glDir=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/results/genotypeVcfGL/
RESULTSDIR=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/results/beagleDifferentFilters/genotypeVcfGL_filtered_DR${dr}
mkdir -p ${RESULTSDIR}
vcfPrefix=genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr
vcfPostfix=.beagle.DR${dr}Filteredgenotype.probs.gg.vcf.gz
for i in {1..22}
do
echo chr$i
echo "#!/bin/bash
#SBATCH --job-name=chr${i}_BIOS_freeze2.1_LLDeep_Diagnostics_FilterGenotypeVcfGL_DR${dr}
#SBATCH --output=chr${i}_BIOS_freeze2.1_LLDeep_Diagnostics_FilterGenotypeVcfGL_DR${dr}.out
#SBATCH --error=chr${i}_BIOS_freeze2.1_LLDeep_Diagnostics_FilterGenotypeVcfGL_DR${dr}.err
#SBATCH --time=5:59:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 10gb
#SBATCH --export=NONE
#SBATCH --get-user-env=30L
#SBATCH --qos=leftover

set -e
set -u

ENVIRONMENT_DIR='.'


module load Python
module load tabix
module list

echo \"## \"\$(date)\" Start \$0\"

if [ ! -f $resultsdir/all_positions_chr${i}.intervals ];
then
  zcat $resultsdir/${vcfPrefix}${i}${vcfPostfix} | \\
            grep -v '^#' | awk '{print \$1 \":\" \$2 \"-\" \$2 }' > $resultsdir/all_positions_chr${i}.intervals
else
  echo "resultsdir/all_positions_chr${i}.intervals already exists, using that"
fi

if [ ! -f $RESULTSDIR/genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr${i}.genotypeGVCF.DR${dr}.gg.vcf ];
then
  python /groups/umcg-bios/tmp03/projects/phasing/selectVariants.py \\
      $glDir//genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr${i}.genotypeGVCF.gg.vcf.gz \\
      $resultsdir/all_positions_chr${i}.intervals \\
      $RESULTSDIR/genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr${i}.genotypeGVCF.DR${dr}.gg.vcf
else 
  echo "$RESULTSDIR/genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr${i}.genotypeGVCF.DR${dr}.gg.vcf already exists, using that"
fi

cd $RESULTSDIR
echo "bgzipping..."
bgzip $RESULTSDIR/genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr${i}.genotypeGVCF.DR${dr}.gg.vcf
echo "tabix..."
tabix $RESULTSDIR/genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr${i}.genotypeGVCF.DR${dr}.gg.vcf.gz
echo "md5sum..."
md5sum $RESULTSDIR/genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr${i}.genotypeGVCF.DR${dr}.gg.vcf.gz > $RESULTSDIR/genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr${i}.genotypeGVCF.DR${dr}.gg.vcf.gz.md5

echo \"## \"\$(date)\" Done \$0\"

">$jobsDir/chr${i}_BIOS_freeze2.1_LLDeep_Diagnostics_FilterGenotypeVcfGL_DR${dr}.sh
done
done
