
jobsDir=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/jobs/convertPLtoGLchr21/
mkdir -p $jobsDir

echo "#!/bin/bash
#SBATCH --job-name=ConvertPLtoGLchr21
#SBATCH --output=ConvertPLtoGLchr21.out
#SBATCH --error=ConvertPLtoGLchr21.err
#SBATCH --time=23:59:00
#SBATCH --cpus-per-task 2
#SBATCH --mem 12gb
#SBATCH --nodes 1
#SBATCH --export=NONE
#SBATCH --get-user-env=30L
#SBATCH --qos=dev

set -e
set -u

ENVIRONMENT_DIR='.'

module load Python
module load ngs-utils/16.12.1
module list

echo \"## \"\$(date)\" Start \$0\"


INPUTDIR=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged/results/filterGQ20_callRate50/
RESULTSDIR="/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/results/PLtoGLchr21/"
mkdir -p \$RESULTSDIR
#Run conversion script beagle vcf to .hap.gzeit format
python \${EBROOTNGSMINUTILS}/PL_to_GL_reorder.py \
    --vcf \$INPUTDIR/genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr21.filter.GQ20_callRate50.PASSonly.BiallelicSNVsOnly.noDiagnostics.gg.vcf.gz \
    --out \$RESULTSDIR/genotypes_BIOSfreeze2.1_LLDeep_noDiagnostics_chr21.vcf.gz

 echo \"returncode: \$?\";
 echo \"succes moving files\";


echo \"## \"\$(date)\" ##  \$0 Done \"

">$jobsDir/ConvertPLtoGLchr21.sh



