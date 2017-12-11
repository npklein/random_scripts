HAPLOTYPEDIR="haplotypic_counts"
module load tabix
INITIALVCFDIR=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged/results/filterGQ20_callRate50_noRnaEditSites/
while read line
do
    if [ "$line" = "CHR,chromosomeChunk" ];
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
    CHUNK=$START.$END
    echo $CHUNK
    features=/apps/data/ftp.nygenome.org/sec/phaser/Homo_sapiens.GRCh37.75.metaGenes.chr$CHR.bed
#    features=/apps/data/ftp.nygenome.org/sec/phaser/hg19_ensembl.chr$CHR.bed
    RESULTSDIR="/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/results/phasing/geneAE_metaGenes/chr$CHR/$CHUNK/"
    INPUTDIR="/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/results/phasing/readbackedPhasing/chr$CHR/$CHUNK/"
    jobsDir=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/jobs/phasing/geneAE_metaGenes/chr$CHR/
    mkdir -p $jobsDir

    echo "#!/bin/bash
#SBATCH --job-name=chr$CHR.$CHUNK.geneAE.metaGenes
#SBATCH --output=geneAE.chr$CHR.$CHUNK.metaGenes.out
#SBATCH --error=geneAE.chr$CHR.$CHUNK.metaGenes.err
#SBATCH --time=05:59:00
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
module load phASER/20170714-cd7daba
module list

echo \"## \"\$(date)\" Start \$0\"
mkdir -p $RESULTSDIR

TMPFEATURES=\$TMPDIR/featureSubset.bed
python ~/random_scripts/genetics/phasing/subset_bedfile.py $features \$TMPFEATURES $CHR $START $END

for haplotype_count in \`ls ${INPUTDIR}/${HAPLOTYPEDIR}/*.txt\`
do
    SAMPLENAME=\$(echo \$haplotype_count | awk -F"." '{print \$3}')
    mkdir -p $RESULTSDIR/
    echo \"haplotype count file: \$haplotype_count\"
    echo \"samplename: \$SAMPLENAME\"
    echo \"start\"
      python \$EBROOTPHASER/phaser_gene_ae/phaser_gene_ae.py \\
                        --haplotypic_counts \$haplotype_count \\
                        --features \$TMPFEATURES \\
                        --o $RESULTSDIR/chr$CHR.$CHUNK.\$SAMPLENAME.geneAE.txt
    echo \"finish\"
done


echo \"returncode: $?\"
echo \"succes moving files\";
echo \"## \"\$(date)\" ##  \$0 Done \"

">$jobsDir/geneAE.chr$CHR.$CHUNK.metaGenes.sh
echo $jobsDir/geneAE.chr$CHR.$CHUNK.metaGenes.sh
done<../geneChunks.20170519.csv
