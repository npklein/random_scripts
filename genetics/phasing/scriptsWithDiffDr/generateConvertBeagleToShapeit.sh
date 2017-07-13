jobDir=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/jobs/phasingDifferentFilters/convertBeagleToShapeit/
mkdir -p $jobDir
for i in {1..22}
do

for dr in $(seq 0 0.1 1); do
INPUTDIR=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/results/beagleDifferentFilters/beagle_vcf_filtered_DR${dr}_gonlOnly/
GLDIR=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/results/beagleDifferentFilters/genotypeVcfGL_filtered_DR${dr}_GoNL//
OUTPUTDIR=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/results/beagleDifferentFilters/beagleToShapeit_DR${dr}/
LIKELIHOODVCF=genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr${i}.genotypeGVCF.GoNL.DR${dr}.gg.ORDERED.vcf.gz
BEAGLEVCF=genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr${i}.beagle.DR${dr}Filteredgenotype.GoNL.probs.gg.vcf.gz
echo "#!/bin/bash
#SBATCH --job-name=chr${i}_BIOS_LLDeep_Diagnostics_freeze2_ConvertBeagle_DR${dr}
#SBATCH --output=chr${i}_BIOS_LLDeep_Diagnostics_freeze2_ConvertBeagle_DR${dr}.out
#SBATCH --error=chr${i}_BIOS_LLDeep_Diagnostics_freeze2_ConvertBeagle_DR${dr}.err
#SBATCH --time=5:59:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 10gb
#SBATCH --export=NONE
#SBATCH --get-user-env=30L
#SBATCH --qos=regular
mkdir -p ${OUTPUTDIR}

set -e
set -u

ENVIRONMENT_DIR='.'


module load GATK
module list

echo \"## \"\$(date)\" Start \$0\"


module load tabix
module load prepareGenFromBeagle4
# Glib is also set as dependency of prepareGenFromBeagle4 but still needs to be loaded after
module load GLib
module list

# had to be gzipped again otherwise it will not read the last line
zcat $INPUTDIR/$BEAGLEVCF \\
    > \$TMPDIR/${BEAGLEVCF%.gz}
gzip \$TMPDIR/${BEAGLEVCF%.gz}
zcat $GLDIR/$LIKELIHOODVCF \\
   > \$TMPDIR/${LIKELIHOODVCF%.gz}
gzip \$TMPDIR/${LIKELIHOODVCF%.gz}

#Run conversion script beagle vcf to .hap.gzeit format
prepareGenFromBeagle4 \
 --likelihoods \$TMPDIR/${LIKELIHOODVCF} \\
 --posteriors \$TMPDIR/${BEAGLEVCF} \\
 --output $OUTPUTDIR/genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr${i}.beagle.DR${dr}.Filteredgenotype

echo "returncode: $?";
# these output files are NOT gzipped, so rename them to filename without gz

echo "succes moving files";




">$jobDir//chr${i}_BIOS_LLDeep_Diagnostics_freeze2_ConvertBeagle_DR${dr}.sh
done
done
