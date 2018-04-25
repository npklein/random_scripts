#!/bin/bash
for CHR in {1..22}
do  # to skip header
    INPUTDIR="/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/results/phasing/gene_matrix/chr$CHR"
    jobsDir=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/jobs/phasing/merge_gene_matrix/
    OUTPUTDIR=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/results/phasing/gene_matrix/
    mkdir -p $OUTPUTDIR
    mkdir -p $jobsDir
    cd $inputDir
        echo "#!/bin/bash
#SBATCH --job-name=merge_gene_matrix.chr$CHR
#SBATCH --output=merge_gene_matrix.chr$CHR.out
#SBATCH --error=merge_gene_matrix.chr$CHR.err
#SBATCH --time=5:59:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 15gb
#SBATCH --nodes 1
#SBATCH --export=NONE
#SBATCH --get-user-env=30L
#SBATCH --qos=dev
set -e
set -u

ENVIRONMENT_DIR='.'

ml purge
module load Python
module list

echo \"## \"\$(date)\" Start \$0\"

python ~/random_scripts/genetics/ASE/merge_geneAE_matrices.py \\
    $INPUTDIR/ \\
    $OUTPUTDIR/BIOS.RNAgenos.chr$CHR

echo \"## \"\$(date)\" ##  \$0 Done \"

">$jobsDir/merge_gene_matrix.chr$CHR.sh
done
