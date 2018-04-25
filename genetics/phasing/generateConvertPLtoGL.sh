
jobsDir=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/jobs/prePhasing/convertPLtoGL/
mkdir -p $jobsDir

for i in {1..22}
do
echo "#!/bin/bash
#SBATCH --job-name=ConvertPLtoGLchr$i
#SBATCH --output=ConvertPLtoGLchr$i.out
#SBATCH --error=ConvertPLtoGLchr$i.err
#SBATCH --time=23:59:00
#SBATCH --cpus-per-task 2
#SBATCH --mem 12gb
#SBATCH --nodes 1
#SBATCH --export=NONE
#SBATCH --get-user-env=30L
#SBATCH --qos=leftover

set -e
set -u

ENVIRONMENT_DIR='.'

module load Biopython/1.65-foss-2015b-Python-3.4.1
module load ngs-utils/17.08.3
module list

echo \"## \"\$(date)\" Start \$0\"


INPUTDIR=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged/results/filterGQ20_callRate50_noRnaEditSites/
RESULTSDIR="/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/results/prePhasing/PLtoGL/"
mkdir -p \$RESULTSDIR
#Run conversion script beagle vcf to .hap.gzeit format
python \${EBROOTNGSMINUTILS}/PL_to_GL_reorder.py \
    --vcf \$INPUTDIR/genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr$i.filter.GQ20_callRate50.BiallelicSNVsOnly.noRnaEditSites.gg.vcf.gz \
    --out \$RESULTSDIR/genotypes_BIOSfreeze2.1_LLDeep_noRnaEditSites_chr$i.vcf.gz

 echo \"returncode: \$?\";
 echo \"succes moving files\";


echo \"## \"\$(date)\" ##  \$0 Done \"

">$jobsDir/ConvertPLtoGLchr$i.sh
done


