echo -e "CHR\tPOS\tchunkSize"
for i in chr21/*.vcf.gz;
do
  numberOfSnps=$(zcat $i | grep -v '^#' | wc -l)
  zcat $i | grep -v '^#' | awk -v numberOfSnps=$numberOfSnps '{print $1 "\t" $2 "\t" numberOfSnps}'
done
