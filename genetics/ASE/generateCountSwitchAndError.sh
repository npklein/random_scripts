for chr in {1..22}
do
  echo chr$chr
  createJobsDir=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/jobs/phasing/countSwitchAndErrors_noRnaEditing/chr$chr
  mkdir -p $createJobsDir
  echo "making jobs in $createJobsDir"
  for CHUNK in /groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing//results/phasing/readbackedPhasing/chr$chr/*
  do
    CHUNKNAME=$(echo $CHUNK | awk -F"/" '{print $12}')
    referenceVCF="/groups/umcg-gonl/tmp03/projects/passBiallelicPhased/mergedVCFs/annotated.with.GoNL.AFandAC/06_IL_haplotype_panel/gonl.chr$chr.snps_indels.r5.3.vcf.gz"
    RESULTSDIR="/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/results/phasing/switchAndErrors/chr$chr/"
    INPUTDIR="/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing//results/phasing/readbackedPhasing/chr$chr/$CHUNKNAME/"
    LINKINGFILE="/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/freeze2_GoNL_related_GTE_30092016_QCpassed.csv"
    SAMPLEFILE="/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/GoNL_sample_IDs.txt"
  echo "#!/bin/bash
#SBATCH --job-name=chr$chr.chunk$CHUNKNAME.countSwapsAndErrors
#SBATCH --output=chr$chr.chunk$CHUNKNAME.countSwapsAndErrors.out
#SBATCH --error=chr$chr.chunk$CHUNKNAME.countSwapsAndErrors.err
#SBATCH --time=5:59:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 4gb
#SBATCH --nodes 1
#SBATCH --export=NONE
#SBATCH --get-user-env=30L
#SBATCH --qos=regular

module load PyVCF/0.6.8-foss-2015b-Python-3.5.1

mkdir -p $RESULTSDIR
python ~/haploSwapper/switchAndErrorCounter.py \\
    --chunkDir $INPUTDIR \\
    --vcfReference $referenceVCF \\
    --linking_file $LINKINGFILE \\
    --samples_file $SAMPLEFILE \\
    --out_file $RESULTSDIR/chr$chr.chunk$CHUNKNAME.switchAndError.txt \\
    --chr $chr
" > $createJobsDir/chr$chr.chunk$CHUNKNAME.countSwapsAndErrors.sh
done
done
