#!/bin/bash

# Merge genes into meta-genes in a GTF such that if
# Gene A      ----------------------------
# Gene B             ------        -----------
# It becomes  AAAAAAAXXXXXXAAAAAAAAXXXXXXXBBBB
# where X is gene A + B
module load Python/3.5.1-foss-2015b
module load BEDTools/2.25.0-foss-2015b
module load BEDOPS/2.4.20
# select gene features only, sort by chr then start position, change to bed file (chr, start, stop first columns),
# partition using bedops, then add extra column info back again
# has to write bed.tmp because bedmap and bedops both have to read the file, don't know how to pipe the input to both
python cutStrangeChr.py /apps/data/ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf \
    | awk '{if ($3 == "gene") print $0;}' \
    | bedtools sort -i stdin \
    | awk  -F $'\t' 'BEGIN {OFS = FS} {print $1, $4, $5, $2, $3, $7, $9}'  \
    > /apps/data/ftp.nygenome.org/sec/phaser/tmp.bed
#bedops --partition tmp.bed \
#    | bedmap --echo --echo-map-id-uniq --delim '\t' - tmp.bed \
#    > /apps/data/ftp.nygenome.org/sec/phaser/Homo_sapiens.GRCh37.75.metaGenes.bed

#rm tmp.bed
