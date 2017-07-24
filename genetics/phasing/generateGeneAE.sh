set -e
set -u

ENVIRONMENT_DIR='.'

ml purge
module load phASER/20170714-cd7daba
module list


#usage: phaser_gene_ae.py [-h] --haplotypic_counts HAPLOTYPIC_COUNTS --features
#                         FEATURES --o O [--id_separator ID_SEPARATOR]
#                         [--gw_cutoff GW_CUTOFF] [--min_cov MIN_COV]
#                         [--min_haplo_maf MIN_HAPLO_MAF]

#optional arguments:
#  -h, --help            show this help message and exit
#  --haplotypic_counts HAPLOTYPIC_COUNTS
#                        Output file from phASER containing read counts for
#                        haplotype blocks. NOTE: unphased_vars must have been
#                        enabled when phASER was run.
#  --features FEATURES   File in BED format (0 BASED COORDINATES -
#                        chr,start,stop,name) containing the features to
#                        produce counts for.
#  --o O                 Output file
#  --id_separator ID_SEPARATOR
#                        Separator used for generating unique variant IDs when
#                        phASER was run.
#  --gw_cutoff GW_CUTOFF
#                        Minimum genome wide phase confidence for phASER
#                        haplotype blocks.
#  --min_cov MIN_COV     Minimum total coverage for a feature to be outputted.
#  --min_haplo_maf MIN_HAPLO_MAF
#                        The minimum MAF used to phase a haplotype for it to be
#                        considered genome wide phased when generating gene
#                        level counts. Setting this number higher will result
#                        in more confident phasing if genotypes were population
#                        prephased. Value must be between 0 and 0.5.
mkdir -p results/gene_ae/
for i in {1..22}
do
  for haplotype_count in chr$i/haplotypic_counts/*.txt
  do
    SAMPLENAME=$(echo $haplotype_count | awk -F"." '{print $2}')
    echo $haplotype_count
    echo $SAMPLENAME
    if [ ! -f results/gene_ae/chr$i.$SAMPLENAME.geneAE.txt ];
    then
      python $EBROOTPHASER/phaser_gene_ae/phaser_gene_ae.py --haplotypic_counts $haplotype_count --features hg19_ensembl.bed \
                           --o results/gene_ae/chr$i.$SAMPLENAME.geneAE.txt
    else
        echo "exists"
    fi
  done
done




