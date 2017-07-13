for dr in $(seq 0 0.1 1); 
do
  echo ${dr}
  jobDir=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/jobs/phasingDifferentFilters/shapeit/DR${dr}/
  mkdir -p $jobDir
  for i in {1..22}
  do
    mkdir -p /groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/jobs/phasingDifferentFilters/shapeit/DR${dr}/chr${i}/
  done
  
 while read line
 do

    #echo $line
    CHR=`echo $line | awk '{print $1}' FS=","`
    START=`echo $line | awk '{print $2}' FS=":" | awk '{print $1}' FS="-"`
    END=`echo $line | awk '{print $2}' FS=":" | awk '{print $2}' FS="-"`
    RESULTSDIR="/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/results/beagleDifferentFilters/shapeit_DR${dr}/chr$CHR/"
    INPUTDIR=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/results/beagleDifferentFilters/beagleToShapeit_DR${dr}/
    INPUTPREFIX=$INPUTDIR/genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr${CHR}.beagle.DR${dr}.Filteredgenotype


    echo "#!/bin/bash
#SBATCH --job-name=chr$CHR.$START.$END.DR${dr}.genotypes_BIOS_LLDeep_Diagnostics_merged_ShapeitPhasing
#SBATCH --output=SBATCH --output=ShapeitPhasing.chr$CHR.$START.$END.DR${dr}.out
#SBATCH --error=ShapeitPhasing.chr$CHR.$START.$END.DR${dr}.err
#SBATCH --time=5:59:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 8gb
#SBATCH --nodes 1
#SBATCH --export=NONE
#SBATCH --get-user-env=30L
#SBATCH --qos=leftover

set -e
set -u

ENVIRONMENT_DIR='.'


module load tabix/0.2.6-foss-2015b
module load shapeit/v2.r837-static
module list

echo \"## \"\$(date)\" Start \$0\"

mkdir -p $RESULTSDIR
#Run shapeit
# The shaping is scaffolded using the chip-based or wgs phased genotypes (--input-init). For data without this information (like
# vcfs from public rnaseq) this pipeline needs to be different OR it needs to be phased together with BIOS samples (using BIOS
# samples as scaffolding, but could give population problems)
shapeit \\
 -call \\
 --input-gen $INPUTPREFIX.gen.gz \\
             $INPUTPREFIX.gen.sample \\
 --input-init $INPUTPREFIX.hap.gz \\
              $INPUTPREFIX.hap.sample \\
 --input-map /apps/data/www.shapeit.fr/genetic_map_b37//genetic_map_chr$CHR\_combined_b37.txt \\
 --input-scaffold /groups/umcg-lld/tmp03/projects/genotypingRelease3/selectionLldeep/lldeepPhased_RNA_IDs//chr_$CHR.haps \\
                  /groups/umcg-lld/tmp03/projects/genotypingRelease3/selectionLldeep/lldeepPhased_RNA_IDs//chr_$CHR.sample \\
 --input-thr 1.0 \\
 --thread 1 \\
 --window 0.1 \\
 --states 400 \\
 --states-random 200 \\
 --burn 0 \\
 --run 12 \\
 --prune 4 \\
 --main 20 \\
 --output-max $RESULTSDIR/genotypes_BIOS_LLDeep_Diagnostics_merged.chr$CHR.$START.$END.DR${dr}.shapeit.phased.hap.gz \\
             $RESULTSDIR/genotypes_BIOS_LLDeep_Diagnostics_merged.chr$CHR.$START.$END.DR${dr}.shapeit.phased.hap.gz.sample \\
 --output-log $RESULTSDIR/genotypes_BIOS_LLDeep_Diagnostics_merged.chr$CHR.$START.$END.DR${dr}.shapeit.phased.log \\
 --input-from $START \\
 --input-to $END

 echo \"returncode: \$?\";
 cd $RESULTSDIR
 bname=\$(basename $RESULTSDIR/genotypes_BIOS_LLDeep_Diagnostics_merged.chr$CHR.$START.$END.DR${dr}.shapeit.phased.hap.gz)
 md5sum \${bname} > \${bname}.md5
 bname=\$(basename $RESULTSDIR/genotypes_BIOS_LLDeep_Diagnostics_merged.chr$CHR.shapeit.phased.hap.gz.sample)
 md5sum \${bname} > \${bname}.md5
 bname=\$(basename $RESULTSDIR/genotypes_BIOS_LLDeep_Diagnostics_merged.chr$CHR.$START.$END.DR${dr}.shapeit.phased.log)
 bname=\$(basename $RESULTSDIR/genotypes_BIOS_LLDeep_Diagnostics_merged.chr$CHR.shapeit.phased.log)
 md5sum \${bname} > \${bname}.md5
 cd -
 echo \"succes moving files\";

touch ShapeitPhasing.chr$CHR.$START.$END.DR${dr}.sh.finished

echo \"## \"\$(date)\" ##  \$0 Done \"

">$jobDir/chr$CHR/ShapeitPhasing.chr$CHR.$START.$END.DR${dr}.sh

done</groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/geneChunks.20170519.csv
done
