#!/bin/bashq

# Merge genes into meta-genes in a GTF such that if
# Gene A      ----------------------
# Gene B                     -----------
# It becomes  AAAAAAAAAAAAAAAXXXXXXXBBBB
# where X is gene A + B

# NOTE: based on Ensemble GTF. This assumes columns split
#       on tabs + spaces where column 1 = CHR, $4 = start,
#       $5 = end, $10 = gene ID (see example.gtf in this dir)


module load Python/3.5.1-foss-2015b
module load BEDTools/2.25.0-foss-2015b
module load BEDOPS/2.4.20
# select gene features only, sort by chr then start position, change to bed file (chr, start, stop first columns),
# partition using bedops, then add extra column info back again
# has to write bed.tmp because bedmap and bedops both have to read the file, don't know how to pipe the input to both
python cutStrangeChr.py /apps/data/ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf \
    | awk '{if ($3 == "gene") print $0;}' \
    | bedtools sort -i stdin \
    | awk  '{print $1 "\t" $4 "\t" $5 "\t" $10}'  \
    | sed 's/"//g' - | sed 's/;//g' \
    > /tmp/tmp.bed
bedops --partition /tmp/tmp.bed \
    | bedmap \
        --echo \
        --echo-map-id-uniq \
        --delim '\t' \
        - /tmp/tmp.bed \
        > /apps/data/ftp.nygenome.org/sec/phaser/Homo_sapiens.GRCh37.75.metaGenes.bed
