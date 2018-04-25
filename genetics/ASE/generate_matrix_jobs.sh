#!/bin/bash
OUTNAME="BIOS.RNAgenos"
for CHR in {1..22}
do  # to skip header
    RESULTSDIR="/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/results/phasing/gene_matrix/chr$CHR/"
    jobsDir=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/jobs/phasing/gene_matrix/chr$CHR/
    inputDir=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/results/phasing/geneAE_metaGenes/chr$CHR/
    mkdir -p $RESULTSDIR
    mkdir -p $jobsDir
    cd $inputDir
    chunks=`ls */`
    cd -
    echo $chunks
    for chunk in $chunks;
    do
        echo "#!/bin/bash
#SBATCH --job-name=create_gene_matrix.chr$CHR.chunk$chunk
#SBATCH --output=create_gene_matrix.chr$CHR.chunk$chunk.out
#SBATCH --error=create_gene_matrix.chr$CHR.chunk$chunk.err
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

python ~/random_scripts/genetics/ASE/create_geneAE_matrix.py \\
    $inputDir/$chunk/ \\
    $RESULTSDIR/$OUTNAME.chr$CHR.chunk$chunk \\
    --removeNoCounts \\
    --chr $CHR

echo \"## \"\$(date)\" ##  \$0 Done \"

">$jobsDir/create_gene_matrix.chr$CHR.chunk$chunk.sh
echo $jobsDir/create_gene_matrix.chr$CHR.chunk$chunk.sh
    done
done
