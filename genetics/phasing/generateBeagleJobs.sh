jobsDir=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/jobs/beagleChunksNoRNAedit/
INPUTDIR=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged/results/filterGQ20_callRate50_noRnaEditSites/
module load tabix
for i in {1..22}
do
  mkdir -p $jobsDir/chr$i/
done

while read line
do

  echo $line
  if [ "$line" = "CHR,chromosomeChunk" ];
  then
    continue
  fi
  CHR=`echo $line | awk '{print $1}' FS=","`
  START=`echo $line | awk '{print $2}' FS=":" | awk '{print $1}' FS="-"`
  END=`echo $line | awk '{print $2}' FS=":" | awk '{print $2}' FS="-"`

  numberOfSnps=$(tabix $INPUTDIR/genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr$CHR.filter.GQ20_callRate50.BiallelicSNVsOnly.noRnaEditSites.gg.vcf.gz $CHR:${START}-${END} | wc -l)
  if [ $numberOfSnps -eq 0 ];
  then
    echo "Chunk $CHR:${START}-${END} did not have any SNPs, skipping"
    continue
  fi
  RESULTSDIR="/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/results/beagleChunksNoRNAedit/chr$CHR/"
  randomNumber=$(( $RANDOM % 2 ))
#  if [ $randomNumber -eq 1 ];
#  then
#     qos="regular"
#  else
     qos="leftover"
#  fi
#  echo "qos: $qos"
echo "#!/bin/bash
#SBATCH --job-name=chr$CHR.$START.$END.genotypes_BIOS_LLDeep_Diagnostics_merged_BeagleGenotyping
#SBATCH --output=BeagleGenotyping.chr$CHR.$START.$END.out
#SBATCH --error=BeagleGenotyping.chr$CHR.$START.$END.err
#SBATCH --time=23:59:00
#SBATCH --cpus-per-task 8
#SBATCH --mem 16gb
#SBATCH --nodes 1
#SBATCH --export=NONE
#SBATCH --get-user-env=30L
#SBATCH --qos=$qos

set -e
set -u

ENVIRONMENT_DIR='.'


module load tabix/0.2.6-foss-2015b
module load beagle/27Jul16.86a-Java-1.8.0_45
module list

echo \"## \"\$(date)\" Start \$0\"

mkdir -p ${RESULTSDIR}


java -Xmx15g -Djava.io.tmpdir=\$TMPDIR -XX:ParallelGCThreads=2 -jar \$EBROOTBEAGLE/beagle.27Jul16.86a.jar \\
 gl=$INPUTDIR/genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr$CHR.filter.GQ20_callRate50.BiallelicSNVsOnly.noRnaEditSites.gg.vcf.gz \\
 out=${RESULTSDIR}/genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr$CHR.filter.GQ20_callRate50.BiallelicSNVsOnly.noRnaEditSites.probs.gg \\
 chrom=${CHR}\:${START}\-${END} \

 echo \"returncode: \$?\";
 echo \"succes moving files\";

touch BeagleGenotyping.chr${CHR}.${START}.${END}.sh.finished

echo \"## \"\$(date)\" ##  \$0 Done \"

">$jobsDir/chr$CHR/BeagleGenotyping.chr$CHR.$START.$END.sh

done</groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/geneChunks.20170519.csv


