mapq=255
baseq=10
for chr in {1..22}
do
  echo chr$chr
  createJobsDir=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/jobs/phasing/generateReadbackedPhasing/
  mkdir -p $createJobsDir
  echo "#!/bin/bash
#SBATCH --job-name=chr$chr.createReadbackJobs
#SBATCH --output=chr$chr.createReadbackJobs.out
#SBATCH --error=chr$chr.createReadbackJobs.err
#SBATCH --time=0:59:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 1gb
#SBATCH --nodes 1
#SBATCH --export=NONE
#SBATCH --get-user-env=30L
#SBATCH --qos=dev

echo \"## \"\$(date)\" Start \$0\"

INDIVIDUALBAMLINK=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/individual_bam_link.txt
NUMBEROFBAMS=\$(wc -l \$INDIVIDUALBAMLINK | awk '{print \$1}')
while read line
do
  CHR=\`echo \$line | awk '{print \$1}' FS=\",\"\`
  if [ ! \"\$CHR\" = \"$chr\" ];
  then
      continue
  fi
  echo \$CHR
  START=\`echo \$line | awk '{print \$2}' FS=\":\" | awk '{print \$1}' FS=\"-\"\`
  END=\`echo \$line | awk '{print \$2}' FS=\":\" | awk '{print \$2}' FS=\"-\"\`
  echo \"chr\$CHR-\$START:\$END\"

  RESULTSDIR=\"/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/results/phasing/readbackedPhasing/chr\$CHR/\$START.\$END/\"
  INPUTDIR=\"/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/results/phasing/shapeitVCFwithAC/chr\$CHR\"
  mkdir -p /groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/jobs/phasing/readbackedPhasing/chr\$CHR/
  BAMSTART=1
  STEPSIZE=200
  for index in \$(seq \$BAMSTART \$STEPSIZE \$NUMBEROFBAMS)
  do
    BAMEND=\$((BAMSTART+STEPSIZE))
    echo \"#!/bin/bash
#SBATCH --job-name=chr\$CHR.\$START.\$END.BAM\$BAMSTART.\$BAMEND.ReadbackPhasing
#SBATCH --output=ReadbackedPhasing.chr\$CHR.\$START.\$END.BAM\$BAMSTART.\$BAMEND.out
#SBATCH --error=ReadbackedPhasing.chr\$CHR.\$START.\$END.BAM\$BAMSTART.\$BAMEND.err
#SBATCH --time=23:59:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 1gb
#SBATCH --nodes 1
#SBATCH --export=NONE
#SBATCH --get-user-env=30L
#SBATCH --qos=dev

set -e
set -u
ENVIRONMENT_DIR='.'

ml purge
module load tabix/0.2.6-foss-2015b
module load phASER/20170714-cd7daba
module load BCFtools/1.3-foss-2015b
module load BEDTools/2.25.0-foss-2015b
module list

echo \\\"## \\\"\\\$(date)\\\" Start \\\$0\\\"
mkdir -p \$RESULTSDIR

mkdir -p \$RESULTSDIR/variant_connections/
mkdir -p \$RESULTSDIR/allelic_counts/
mkdir -p \$RESULTSDIR/haplotypes/
mkdir -p \$RESULTSDIR/haplotypic_counts/
mkdir -p \$RESULTSDIR/allele_config/
mkdir -p \$RESULTSDIR/vcf_per_sample/
mkdir -p \$RESULTSDIR/vcf_per_sample_no_counts/

VCFDIR=\\\"/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/results/phasing/shapeitVCFwithAC//chr\$CHR\\\"
VCFNAME=\\\"genotypes_BIOS_LLDeep_Diagnostics_merged.chr\$CHR.\$START.\$END.shapeit.phased.withAC.vcf.gz\\\"
VCF=\\\"\\\$VCFDIR/\\\$VCFNAME\\\"
echo "BAMSTART: \$BAMSTART"
echo "BAMEND: \$BAMEND"
INDEX=0
while read individual_bam_link; do

  if [ \\\"\$individual_bam_link\\\" == \\\"individualID,sampleName,bam\\\" ];
  then
    continue
  fi

  INDEX=\\\$((INDEX+1))
  echo "INDEX: \\\$INDEX"
  if [ "\\\$INDEX" -lt "\$BAMSTART" ];
  then
    continue
  fi
  if [ "\\\$INDEX" -gt "\$BAMEND" ];
  then
    echo "included all BAMS for this chunk"
    exit;
  fi

  bam=\\\$(echo \\\$individual_bam_link | awk -F\\\",\\\" '{ print \\\$3 }')
  echo \\\$bam
  bname=\\\$(basename \\\$bam)
  SAMPLENAME=\\\$(echo \\\$individual_bam_link | awk -F\\\",\\\" '{ print \\\$2}')

  if ! grep -Fxq \\\"\\\$SAMPLENAME\\\" /groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/samples_in_vcf.txt;
  then
      echo \\\"Sample '\\\$SAMPLENAME' not in VCF\\\"
      continue
  fi
  echo \\\"Phasing sample \\\$SAMPLENAME\\\"
  phaserOutPrefix=\$RESULTSDIR/BIOS_LLDeep_Diagnostics_phASER.\\\$SAMPLENAME.chr\$CHR.\$START.\$END
  #Set output prefix per sample for statistics etc.
  python \\\$EBROOTPHASER/phaser/phaser.py \\\\
    --paired_end 1 \\\\
    --bam \\\$bam \\\\
    --vcf \\\$VCF \\\\
    --mapq $mapq \\\\
    --sample \\\$SAMPLENAME \\\\
    --baseq $baseq \\\\
    --o \\\$phaserOutPrefix \\\\
    --temp_dir \\\$TMPDIR \\\\
    --threads 1 \\\\
    --gw_phase_method 1 \\\\
    --chr \$CHR \\\\
    --gw_phase_vcf 1 \\\\
    --show_warning 1 \\\\
    --debug 1 \\\\
    --blacklist /apps/data/ftp.nygenome.org/sec/phaser/hg19_hla.bed.gz \\\\
    --haplo_count_blacklist /apps/data/ftp.nygenome.org/sec/phaser/hg19_haplo_count_blacklist.chr\$CHR.bed.gz

  # phaser does't send appropriate exit signal so try like this
  if tail -n5 ReadbackedPhasing.chr\$CHR.\$START.\$END.BAM\$BAMSTART.\$BAMEND.out | grep -q ERROR;
  then
      if tail -n5 ReadbackedPhasing.chr\$CHR.\$START.\$END.BAM\$BAMSTART.\$BAMEND.out | grep -q \\\"No heterozygous sites that passed all filters were included in the analysis\\\";
      then
        # SHOULD THIS BE CHANGED TO REMOVE THE EMPTY FILES?
        echo \\\"No heterozygous sites to phase\\\"
      elif tail -n5 ReadbackedPhasing.chr\$CHR.\$START.\$END.BAM\$BAMSTART.\$BAMEND.out | grep -q \\\"No reads could be matched to variants\\\";
      then
        echo \\\"No reads could be matched to variants, probably because none passed baseq andor mapq threshold.\\\"
      else
          echo \\\"returncode: 1\\\";
          echo \\\"fail\\\";
          >&2 echo \\\"Tried to use bam '\\\$bam'\\\"
          >&2 echo \\\"Tried to use VCF '\\\$VCF'\\\"
          >&2 echo \\\"exit, phASER error\\\"
          exit 1;
      fi
      echo \\\"Just writing VCF away without readbackedphasing\\\"
      VCFOUT=\\\$(basename \\\$phaserOutPrefix)
      bcftools view --samples \\\$SAMPLENAME \\\$VCF > \$RESULTSDIR/\\\${VCFOUT}.vcf
      bgzip -f \$RESULTSDIR/\\\${VCFOUT}.vcf
      tabix -f \$RESULTSDIR/\\\${VCFOUT}.vcf.gz
      mv \$RESULTSDIR/\\\${VCFOUT}.vcf.gz* \$RESULTSDIR/vcf_per_sample_no_counts/
  else
    echo \\\"Success: moving log files from \\\$phaserOutPrefix to corresponding directories\\\"
    mv \\\$phaserOutPrefix.variant_connections.txt \$RESULTSDIR/variant_connections/
    mv \\\$phaserOutPrefix.allelic_counts.txt \$RESULTSDIR/allelic_counts/
    mv \\\$phaserOutPrefix.haplotypes.txt \$RESULTSDIR/haplotypes/
    mv \\\$phaserOutPrefix.haplotypic_counts.txt \$RESULTSDIR/haplotypic_counts/
    mv \\\$phaserOutPrefix.allele_config.txt \$RESULTSDIR/allele_config/
    mv \\\$phaserOutPrefix.vcf.gz \$RESULTSDIR/vcf_per_sample/
    mv \\\$phaserOutPrefix.vcf.gz.tbi \$RESULTSDIR/vcf_per_sample/
  fi
done<\$INDIVIDUALBAMLINK

echo \\\"phaser done\\\"
echo \\\"returncode: \\\$?\\\"
echo \\\"succes moving files\\\";
echo \\\"## \\\"\\\$(date)\\\" ##  \\\$0 Done \\\"

\">/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/jobs/phasing/readbackedPhasing/chr\$CHR/ReadbackedPhasing.chr\$CHR.\$START.\$END.BAM\$BAMSTART.\$BAMEND.sh
BAMSTART=\$((BAMEND+1))
done
done</groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/beagledGeneChunks.18072017.csv
" > $createJobsDir/chr$chr.createReadbackJobs.sh
done

