mapq=10
baseq=10
OneKgPhase3VCF="/apps/data/1000G/phase3/20130502//ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz"


mkdir -p /groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/jobs/readbackedPhasingGoNLWGS/chr$1/

for CHR in {1..22}
do  # to skip header
  RESULTSDIR="/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/results/readbackedPhasingGoNLWGS/chr$CHR/"
  mkdir -p /groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/jobs/readbackedPhasingGoNLWGS/chr$CHR/
  SKIPPED=0

  while read individual_bam_link; do
    if [ "$individual_bam_link" == "individualID,sampleName,bam" ];
    then
        continue
    fi
    bam=$(echo $individual_bam_link | awk -F"," '{ print $3 }')
    bname=$(basename $bam)
    SAMPLENAME=$(echo $individual_bam_link | awk -F"," '{ print $2}')

    if ! grep -Fxq "$SAMPLENAME" GoNL_variants/gonlSamplesRNAids.txt
    then
#        echo "$SAMPLENAME not in GoNL"
        continue
    fi
    if ! grep -Fxq "$SAMPLENAME" /groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/samples_in_vcf.txt;
    then
        echo "$SAMPLENAME not in VCF"
        continue
    fi
    echo "$SAMPLENAME in GoNL, make job"

    WGSID=$(grep -w $SAMPLENAME GoNL_variants/linking_file.txt | awk '{ print $1 }')

    echo "#!/bin/bash
#SBATCH --job-name=$SAMPLENAME.chr$CHR.ReadbackPhasing
#SBATCH --output=ReadbackedPhasing.$SAMPLENAME.chr$CHR.out
#SBATCH --error=ReadbackedPhasing.$SAMPLENAME.chr$CHR.err
#SBATCH --time=05:59:00
#SBATCH --cpus-per-task 4
#SBATCH --mem 1gb
#SBATCH --nodes 1
#SBATCH --export=NONE
#SBATCH --get-user-env=30L
#SBATCH --qos=regular

set -e
set -u

ENVIRONMENT_DIR='.'

ml purge
module load tabix/0.2.6-foss-2015b
module load phASER/20170606-e75699d
module load BCFtools/1.3-foss-2015b
module load BEDTools/2.25.0-foss-2015b
module list

echo \"## \"\$(date)\" Start \$0\"
mkdir -p $RESULTSDIR

mkdir -p $RESULTSDIR/variant_connections/
mkdir -p $RESULTSDIR/allelic_counts/
mkdir -p $RESULTSDIR/haplotypes/
mkdir -p $RESULTSDIR/haplotypic_counts/
mkdir -p $RESULTSDIR/allele_config/
mkdir -p $RESULTSDIR/vcf_per_sample/

VCFDIR=\"/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/GoNL_variants/\"
VCFNAME=\"gonl-abc_samples.chr$CHR.release5.NoMAFSelection.vcf.gz\"
VCF=\"\$VCFDIR/\$VCFNAME\"

phaserOutPrefix=$RESULTSDIR/BIOS_LLDeep_Diagnostics_phASER.$SAMPLENAME.chr$CHR
#Set output prefix per sample for statistics etc.
output=\$(python \$EBROOTPHASER/phaser/phaser.py \\
    --paired_end 1 \\
    --bam $bam \\
    --vcf \$VCF \\
    --mapq $mapq \\
    --sample $WGSID \\
    --baseq $baseq \\
    --o \$phaserOutPrefix \\
    --temp_dir \$TMPDIR \\
    --threads 1 \\
    --gw_phase_method 1 \\
    --chr $CHR \\
    --gw_af_vcf $OneKgPhase3VCF \\
    --gw_phase_vcf 1 \\
    --show_warning 1 \\
    --debug 1)

# blacklist takes too long
#    --blacklist /apps/data/ftp.nygenome.org/sec/phaser/hg19_haplo_count_blacklist.bed.gz \\

# phaser does't send appropriate exit signal so try like this
if echo \$output | grep -q ERROR;
then
    >&2 echo \"\$output\"
    if echo \$output | grep -q \"No heterozygous sites that passed all filters were included in the analysis\";
    then
        echo "No heterozygous sites to phase, so just writing sample away without readbackedphasing"
        VCFOUT=\$(basename \$phaserOutPrefix)
        bcftools view --samples $SAMPLENAME \$VCF > $RESULTSDIR/vcf_per_sample/\${VCFOUT}.vcf
        bgzip $RESULTSDIR/vcf_per_sample/\${VCFOUT}.vcf
        tabix $RESULTSDIR/vcf_per_sample/\${VCFOUT}.vcf.gz
    elif echo \$output | grep -q \"No reads could be matched to variants\";
    then
        echo "No reads could be matched to variants, probably because none passed baseq andor mapq threshold. Just writing sample away without readbackedphasing"
        VCFOUT=\$(basename \$phaserOutPrefix)
        bcftools view --samples $SAMPLENAME \$VCF > $RESULTSDIR/vcf_per_sample/\${VCFOUT}.vcf
        bgzip $RESULTSDIR/vcf_per_sample/\${VCFOUT}.vcf
        tabix $RESULTSDIR/vcf_per_sample/\${VCFOUT}.vcf.gz
    else
        echo \"returncode: 1\";
        echo \"fail\";
        >&2 echo \"exit, phASER error\"
        exit 1;
    fi
else
    echo \"\$output\"
    echo \" moving log files from \$phaserOutPrefix to corresponding directories\"
    mv \$phaserOutPrefix.variant_connections.txt $RESULTSDIR/variant_connections/
    mv \$phaserOutPrefix.allelic_counts.txt $RESULTSDIR/allelic_counts/
    mv \$phaserOutPrefix.haplotypes.txt $RESULTSDIR/haplotypes/
    mv \$phaserOutPrefix.haplotypic_counts.txt $RESULTSDIR/haplotypic_counts/
    mv \$phaserOutPrefix.allele_config.txt $RESULTSDIR/allele_config/
    mv \$phaserOutPrefix.vcf.gz $RESULTSDIR/vcf_per_sample/
    mv \$phaserOutPrefix.vcf.gz.tbi $RESULTSDIR/vcf_per_sample/
fi


echo \"phaser done\"
echo \"returncode: $?\"
echo \"succes moving files\";
echo \"## \"\$(date)\" ##  \$0 Done \"

">/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/jobs/readbackedPhasingGoNLWGS/chr$CHR/ReadbackedPhasing.$SAMPLENAME.chr$CHR.sh
   done</groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/individual_bam_link.txt

done</groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/phasedGeneChunks.21062017.csv


