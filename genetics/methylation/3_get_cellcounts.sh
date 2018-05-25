# niekdeklein@gmail.com
OUTDIR=/groups/umcg-lld/tmp03/umcg-ndeklein/EWAS/output/CPACOR
mkdir -p $OUTDIR

module load R/3.3.3-foss-2015b
module list

R --vanilla << "EOF"

setwd('/groups/umcg-lld/tmp03/umcg-ndeklein/EWAS/output/CPACOR')

require(minfi)
require(IlluminaHumanMethylation450kmanifest)

filenames <- Sys.glob("/groups/umcg-lld/tmp03/rawdata/IlluminaHumanMethylation450k/*/*_Grn.idat")
fileBaseNames <- unlist(strsplit(filenames,'_Grn.idat'))
print('read data')
RGset <- read.metharray(fileBaseNames, verbose=TRUE)
counts <- estimateCellCounts (RGset)
write.table(counts,file="cellcounts.txt")
print(paste('Saved to',getwd()))
EOF

