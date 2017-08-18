#!/bin/bashq

# Merge genes into meta-genes in a GTF such that if
# Gene A      ----------------------
# Gene B                     -----------
# It becomes  AAAAAAAAAAAAAAAXXXXXXXBBBB
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
    | awk  '{print $1 "\t" $4 "\t" $5 "\t" $10}'  \
    | sed 's/"//g' - | sed 's/;//g' \
    > /apps/data/ftp.nygenome.org/sec/phaser/tmp.bed
bedops --partition /apps/data/ftp.nygenome.org/sec/phaser/tmp.bed \
    | bedmap \
        --echo \
        --echo-map-id-uniq \
        --delim '\t' \
        - /apps/data/ftp.nygenome.org/sec/phaser/tmp.bed \
        > /apps/data/ftp.nygenome.org/sec/phaser/Homo_sapiens.GRCh37.75.metaGenes.bed

rm tmp.bed
