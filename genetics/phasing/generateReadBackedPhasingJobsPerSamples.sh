module load tabix
mapq=255
baseq=10
INITIALVCFDIR=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged/results/filterGQ20_callRate50_noRnaEditSites/
INDIVIDUALBAMLINK=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/individual_bam_link.txt
CHUNKFILE=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/geneChunks.20170519.csv
jobsDir="/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/jobs/phasing/readbackedPhasingPerSample/"
RESULTSDIR="/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/results/phasing/readbackedPhasing/"
INPUTDIR="/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/results/phasing/shapeitVCFwithAC/"
mkdir -p $jobsDir

# make directories here because it takes too long in the jobs
while read line;
do
  if [ "$line" == "CHR,chromosomeChunk" ];
  then
      continue
  fi
  CHR=`echo $line | awk '{print $1}' FS=","`
  START=`echo $line | awk '{print $2}' FS=":" | awk '{print $1}' FS="-"`
  END=`echo $line | awk '{print $2}' FS=":" | awk '{print $2}' FS="-"`
  numberOfSnps=$(tabix $INITIALVCFDIR/genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr$CHR.filter.GQ20_callRate50.BiallelicSNVsOnly.noRnaEditSites.gg.vcf.gz $CHR:${START}-${END} | wc -l)
  if [ $numberOfSnps -eq 0 ];
  then
    echo "Chunk $CHR:${START}-${END} did not have any SNPs, skipping"
    continue
  fi
  echo "Making dirs for chr$CHR-$START:$END"
  RESULTSDIRCHUNK=$RESULTSDIR/chr$CHR/$START.$END
  mkdir -p $RESULTSDIRCHUNK/variant_connections/
  mkdir -p $RESULTSDIRCHUNK/allelic_counts/
  mkdir -p $RESULTSDIRCHUNK/haplotypes/
  mkdir -p $RESULTSDIRCHUNK/haplotypic_counts/
  mkdir -p $RESULTSDIRCHUNK/allele_config/
  mkdir -p $RESULTSDIRCHUNK/vcf_per_sample/
  mkdir -p $RESULTSDIRCHUNK/vcf_per_sample_no_counts/
done<$CHUNKFILE

while read individual_bam_link
do
  if [ "$individual_bam_link" == "individualID,sampleName,bam" ];
  then
    continue
  fi

  bam=$(echo $individual_bam_link | awk -F"," '{ print $3 }')
  bamBasename=$(basename $bam)
  SAMPLENAME=$(echo $individual_bam_link | awk -F"," '{ print $2}')

  if ! grep -Fxq "$SAMPLENAME" /groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/samples_in_vcf.txt;
  then
      echo "Sample '$SAMPLENAME' not in VCF"
      continue
  fi
  echo "#!/bin/bash
#SBATCH --job-name=$SAMPLENAME.ReadbackPhasing
#SBATCH --output=ReadbackedPhasing.$SAMPLENAME.out
#SBATCH --error=ReadbackedPhasing.$SAMPLENAME.err
#SBATCH --time=23:59:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 4gb
#SBATCH --nodes 1
#SBATCH --export=NONE
#SBATCH --get-user-env=30L
#SBATCH --qos=regular

set -e
set -u
ENVIRONMENT_DIR='.'

ml purge
module load tabix/0.2.6-foss-2015b
module load phASER/20171017-7f5377f
module load BCFtools/1.3-foss-2015b
module load BEDTools/2.25.0-foss-2015b
module list

echo \"## \"\$(date)\" Start \$0\"
TMPBAM=\$TMPDIR/$bamBasename
TMPBLACKLIST=\$TMPDIR/hg19_hla.bed.gz
rsync -vP ${bam} \$TMPBAM
rsync -vP ${bam%bam}bai \${TMPBAM}.bai
rsync -vP /apps/data/ftp.nygenome.org/sec/phaser/hg19_hla.bed.gz \$TMPBLACKLIST

while read line;
do
  if [ \"\$line\" == \"CHR,chromosomeChunk\" ];
  then
      continue
  fi
  CHR=\`echo \$line | awk '{print \$1}' FS=\",\"\`
  START=\`echo \$line | awk '{print \$2}' FS=\":\" | awk '{print \$1}' FS=\"-\"\`
  END=\`echo \$line | awk '{print \$2}' FS=\":\" | awk '{print \$2}' FS=\"-\"\`
  numberOfSnps=\$(tabix $INITIALVCFDIR/genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr\$CHR.filter.GQ20_callRate50.BiallelicSNVsOnly.noRnaEditSites.gg.vcf.gz \$CHR:\${START}-\${END} | wc -l)
  if [ \$numberOfSnps -eq 0 ];
  then
    echo \"Chunk \$CHR:\${START}-\${END} did not have any SNPs, skipping\"
    continue
  fi
  echo \"Phasing chr\$CHR-\$START:\$END\"
  RESULTSDIRCHUNK=$RESULTSDIR/chr\$CHR/\$START.\$END
  VCFDIR=\"/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/results/phasing/shapeitVCFwithAC//chr\$CHR\"
  VCFNAME=\"genotypes_BIOS_LLDeep_noRNAedit_merged.chr\$CHR.\$START.\$END.shapeit.phased.withAC.vcf.gz\"
  VCF=\"\$VCFDIR/\$VCFNAME\"
  TMPVCF=\"\$TMPDIR/\$VCFNAME\"
  rsync -vP \$VCF \$TMPVCF
  rsync -vP \$VCF.tbi \$TMPVCF.tbi
  phaserOutPrefix=\$RESULTSDIRCHUNK/BIOS_LLDeep_noRNAeditSites_phASER.$SAMPLENAME.chr\$CHR.\$START.\$END
  VCFOUT=\$(basename \$phaserOutPrefix)
  echo "Running \$VCF..."
#  if [ -f \$RESULTSDIRCHUNK/vcf_per_sample/\${VCFOUT}.vcf.gz ] || [ -f \$RESULTSDIRCHUNK/vcf_per_sample_no_counts/\${VCFOUT}.vcf.gz ];
#  then
#    echo \"\$RESULTSDIRCHUNK/vcf_per_sample/\${VCFOUT}.vcf.gz or  \$RESULTSDIRCHUNK/vcf_per_sample_no_counts/\${VCFOUT}.vcf.gz already exists, skipping\"
#    continue
#  fi
  #Set output prefix per sample for statistics etc.
   python \$EBROOTPHASER/phaser/phaser.py \\
    --paired_end 1 \\
    --bam \$TMPBAM \\
    --vcf \$TMPVCF \\
    --mapq $mapq \\
    --sample $SAMPLENAME \\
    --baseq $baseq \\
    --o \$phaserOutPrefix \\
    --temp_dir \$TMPDIR \\
    --threads 1 \\
    --gw_phase_method 1 \\
    --chr \$CHR \\
    --gw_phase_vcf 1 \\
    --show_warning 1 \\
    --debug 1 \\
    --blacklist \$TMPBLACKLIST \\
    --haplo_count_blacklist /apps/data/ftp.nygenome.org/sec/phaser/hg19_haplo_count_blacklist.chr\$CHR.bed.gz

    echo \"moving log files from \$phaserOutPrefix to corresponding directories\"
    # if there are no heteroyzgous sites there is no output to move, so use || true that it doesnt error on move
    mv \$phaserOutPrefix.variant_connections.txt \$RESULTSDIRCHUNK/variant_connections/ || true
    mv \$phaserOutPrefix.allelic_counts.txt \$RESULTSDIRCHUNK/allelic_counts/ || true
    mv \$phaserOutPrefix.haplotypes.txt \$RESULTSDIRCHUNK/haplotypes/ || true
    mv \$phaserOutPrefix.haplotypic_counts.txt \$RESULTSDIRCHUNK/haplotypic_counts/ || true
    mv \$phaserOutPrefix.allele_config.txt \$RESULTSDIRCHUNK/allele_config/ || true
    mv \$phaserOutPrefix.vcf.gz \$RESULTSDIRCHUNK/vcf_per_sample/ || \
            (bcftools view --samples \$SAMPLENAME -Oz -o \$RESULTSDIRCHUNK/vcf_per_sample_no_counts/\${VCFOUT}.vcf \$VCF \
            && bgzip \$RESULTSDIRCHUNK/vcf_per_sample_no_counts/\${VCFOUT}.vcf && tabix \$RESULTSDIRCHUNK/vcf_per_sample_no_counts/\${VCFOUT}.vcf.gz)
    mv \$phaserOutPrefix.vcf.gz.tbi \$RESULTSDIRCHUNK/vcf_per_sample/ || true

done<$CHUNKFILE

echo \"phaser done\"
echo \"returncode: \$?\"
echo \"succes moving files\";
echo \"## \"\$(date)\" ##  \$0 Done \"

">$jobsDir/ReadbackedPhasing.$SAMPLENAME.sh
echo $jobsDir/ReadbackedPhasing.$SAMPLENAME.sh
done<$INDIVIDUALBAMLINK



