jobsDir=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/jobs/filterGLchr21
mkdir -p $jobsDir

INPUTDIR=RESULTSDIR="/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/results/PLtoGLchr21/"
glDir=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/results/genotypeVcfGL/
RESULTSDIR=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/results/beagleDifferentFilters/genotypeVcfGL_filtered
mkdir -p ${RESULTSDIR}
vcfPrefix=genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr
vcfPostfix=.beagle.Filteredgenotype.probs.gg.vcf.gz
for i in {21..21}
do
echo chr$i
echo "#!/bin/bash
#SBATCH --job-name=chr${i}_BIOS_freeze2.1_LLDeep_Diagnostics_FilterGenotypeVcfGL
#SBATCH --output=chr${i}_BIOS_freeze2.1_LLDeep_Diagnostics_FilterGenotypeVcfGL.out
#SBATCH --error=chr${i}_BIOS_freeze2.1_LLDeep_Diagnostics_FilterGenotypeVcfGL.err
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

if [ ! -f $INPUTDIR/all_positions_chr${i}.intervals ];
then
  zcat $INPUTDIR/${vcfPrefix}${i}${vcfPostfix} | \\
            grep -v '^#' | awk '{print \$1 \":\" \$2 \"-\" \$2 }' > $INPUTDIR/all_positions_chr${i}.intervals
else
  echo "$INPUTDIR/all_positions_chr${i}.intervals already exists, using that"
fi

if [ ! -f $RESULTSDIR/genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr${i}.genotypeGVCF.gg.vcf ];
then
  python /groups/umcg-bios/tmp03/projects/phasing/selectVariants.py \\
      $glDir//genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr${i}.genotypeGVCF.gg.vcf.gz \\
      $INPUTDIR/all_positions_chr${i}.intervals \\
      $RESULTSDIR/genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr${i}.genotypeGVCF.gg.vcf
else 
  echo "$RESULTSDIR/genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr${i}.genotypeGVCF.gg.vcf already exists, using that"
fi

cd $RESULTSDIR
echo "bgzipping..."
bgzip $RESULTSDIR/genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr${i}.genotypeGVCF.gg.vcf
echo "tabix..."
tabix $RESULTSDIR/genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr${i}.genotypeGVCF.gg.vcf.gz
echo "md5sum..."
md5sum $RESULTSDIR/genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr${i}.genotypeGVCF.gg.vcf.gz > $RESULTSDIR/genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr${i}.genotypeGVCF.gg.vcf.gz.md5

echo \"## \"\$(date)\" Done \$0\"

">$jobsDir/chr${i}_BIOS_freeze2.1_LLDeep_Diagnostics_FilterGenotypeVcfGL.sh
done
