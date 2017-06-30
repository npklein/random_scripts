name=phasedGeneChunks.21062017.csv
echo "CHR,chromosomeChunk" > $name
ls results/shapeit/chr*/*sample | awk 'FS="." {print $2 "," $2 ":" $3 "-" $4}' >> $name
sed -i 's;chr;;g' $name
