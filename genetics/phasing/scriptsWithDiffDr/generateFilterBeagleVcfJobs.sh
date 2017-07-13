for i in {1..22}
do

for dr in $(seq 0 0.1 1); do 
mkdir -p jobs/phasingDifferentFilters/beagleDifferentFilters/beagle_vcf_filtered_${dr}/
echo "#!/bin/bash
#SBATCH --job-name=chr${i}_BIOS_LLDeep_Diagnostics_freeze2_FilterBeagle_DR${dr}
#SBATCH --output=chr${i}_BIOS_LLDeep_Diagnostics_freeze2_FilterBeagle_DR${dr}.out
#SBATCH --error=chr${i}_BIOS_LLDeep_Diagnostics_freeze2_FilterBeagle_DR${dr}.err
#SBATCH --time=5:59:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 10gb
#SBATCH --export=NONE
#SBATCH --get-user-env=30L
#SBATCH --qos=leftover
mkdir -p /groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/results/beagleDifferentFilters/beagle_vcf_filtered_DR${dr}/

set -e
set -u

ENVIRONMENT_DIR='.'


module load GATK
module list

echo \"## \"\$(date)\" Start \$0\"


java -Xmx8g -jar \${EBROOTGATK}/GenomeAnalysisTK.jar \\
  -T SelectVariants \\
  -R /apps/data/ftp.broadinstitute.org/bundle/2.8/b37/human_g1k_v37.fasta \\
  -V /groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/results/beagle//genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr${i}.beagle.genotype.probs.gg.vcf.gz \\
  -o /groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/results/beagleDifferentFilters/beagle_vcf_filtered_DR${dr}/genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr${i}.beagle.DR${dr}Filteredgenotype.probs.gg.vcf.gz \\
  -select \"DR2 > $dr\"

">./jobs/phasingDifferentFilters/beagleDifferentFilters/beagle_vcf_filtered_${dr}/chr${i}_BIOS_LLDeep_Diagnostics_freeze2_FilterBeagle_DR${dr}.sh
done
done
