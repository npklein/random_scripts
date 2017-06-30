RESULTSDIR=/groups/umcg-bios/tmp03/projects/phasing/results_GQ20/beagle_vcf_filtered
mkdir -p /groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/results/beagle_vcf_filtered/
for i in {1..22}
do


echo "#!/bin/bash
#SBATCH --job-name=chr${i}_BIOS_LLDeep_Diagnostics_freeze2_FilterBeagle
#SBATCH --output=chr${i}_BIOS_LLDeep_Diagnostics_freeze2_FilterBeagle.out
#SBATCH --error=chr${i}_BIOS_LLDeep_Diagnostics_freeze2_FilterBeagle.err
#SBATCH --time=5:59:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 10gb
#SBATCH --export=NONE
#SBATCH --get-user-env=30L
#SBATCH --qos=dev

set -e
set -u

ENVIRONMENT_DIR='.'


module load GATK
module list

echo \"## \"\$(date)\" Start \$0\"


java -Xmx8g -jar \${EBROOTGATK}/GenomeAnalysisTK.jar \\
  -T SelectVariants \\
  -R /apps/data/ftp.broadinstitute.org/bundle/2.8/b37/human_g1k_v37.fasta \\
  -V /groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/results/beagle/genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr${i}.beagle.genotype.probs.gg.vcf.gz \\
  -o /groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/results/beagle_vcf_filtered/genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr${i}.beagle.DR2Filteredgenotype.probs.gg.vcf.gz \\
  -select \"DR2 > 0.989999999999999999\"

">./jobs/beagle_vcf_filtered/chr${i}_BIOS_LLDeep_Diagnostics_freeze2_FilterBeagle_FilterBeagle.sh

done
