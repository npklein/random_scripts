import glob
import os

directory = 'merged'
if not os.path.exists(directory):
    os.makedirs(directory)

print('reading')
gene_info = {}
geneCounts = glob.glob('chr*/*/*.txt')
genecountsLength = str(len(geneCounts))
x = 0
header = ''
for f in geneCounts:
    if x % 1000 == 0:
        print(str(x)+'/'+genecountsLength)
    x+=1
    sample = f.split('.')[4]
    if sample not in gene_info:
        gene_info[sample] = {}
    with open(f) as input_file:
        header = input_file.readline() #header
        for line in input_file:
            line = line.strip().split('\t')
            # line[6] == totalCount
            if int(line[6]) != 0:
                gene = line[3]
                if gene in gene_info[sample]:
                    raise RuntimeError('2x same gene has counts for this sample.\nsample: '+sample+'\ngene: '+gene)
                gene_info[sample][gene] = line
print('finished reading')

print('writing')
for sample in gene_info:
    print(sample)
    with open(directory+'/'+sample+'.geneAE.txt','w') as out:
        out.write(header)
        for gene in sorted(gene_info[sample]):
            out.write('\t'.join(gene_info[sample][gene])+'\n')
