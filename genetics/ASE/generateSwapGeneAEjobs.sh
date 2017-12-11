for CHR in {1..22}
do
    for CHUNK in $(ls /groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/results/phasing/readbackedPhasing/chr$CHR/)
    do
      jobsDir="/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/jobs/phasing/swapGeneAeCounts/chr$CHR/"
      RESULTSDIR="/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/results/phasing/geneAE_switched//chr$CHR/chunk$CHUNK/"
      GENEAEDIR="/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/results/phasing/geneAE_metaGenes/chr$CHR/$CHUNK/"
      SWITCHDIR="/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/results/phasing/switchAndErrors/chr$CHR/"
      mkdir -p $jobsDir
    echo "#!/bin/bash
#SBATCH --job-name=swapGeneAE.chr$CHR.$CHUNK
#SBATCH --output=swapGeneAE.chr$CHR.$CHUNK.out
#SBATCH --error=swapGeneAE.chr$CHR.$CHUNK.err
#SBATCH --time=5:59:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 5gb
#SBATCH --nodes 1
#SBATCH --export=NONE
#SBATCH --get-user-env=30L
#SBATCH --qos=regular
set -e
set -u

module load PyVCF/0.6.8-foss-2015b-Python-3.5.1

  mkdir -p $RESULTSDIR
  python ~/haploSwapper/swapGeneAE.py \\
    --switchErrorDir $SWITCHDIR \\
    --geneAeCountsDir $GENEAEDIR \\
    --out_dir $RESULTSDIR \\
    --chunk $CHUNK
" > $jobsDir/swapGeneAE.chr$CHR.$CHUNK.sh
    echo $jobsDir/swapGeneAE.chr$CHR.$CHUNK.sh
done
done
