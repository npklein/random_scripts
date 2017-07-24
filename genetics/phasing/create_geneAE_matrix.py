import glob
import argparse
import os

parser = argparse.ArgumentParser(description='Make matrix (genes x samples) of the gene counts and total depth for all results in given dir + subdirs. IMPORTANT: sample names are taken by splitting bam name on . and taking first [0] element.')
parser.add_argument('geneAE_rootdir', help='Rootdir where all gene counts can be found. Will walk through subdirs and use all *.txt files')
parser.add_argument('out_prefix', help='Write output to <out_prefix>.geneCounts.txt and <out_prefix>.totalDepth.txt')
parser.add_argument('--verbose', action='store_true',help='If set, print out which files are being used')

args = parser.parse_args()
print('Combining gene count files from '+args.geneAE_rootdir)
print('Writing results to '+args.out_prefix)

print('Start combining...')

gene_data = {}
gene_names = []
genes = []
print('reading data')
#contig    start    stop    name    aCount    bCount    totalCount    log2_aFC    n_variants    variants    gw_phased    bam
for subdir, dirs, files in os.walk(args.geneAE_rootdir):
    for file in files:
        if not file.endswith('.txt') or not 'chr22' in file:
            continue
        geneAE_file = subdir+'/'+file
        if args.verbose:
            print(geneAE_file)
        with open(geneAE_file) as input_file:
            input_file.readline()
            input_file.readline()
            for line in input_file:
                line = line.strip().split('\t')
                name = line[-1].split('.')[0]
                if not name in gene_data:
                    gene_data[name] = {}
                    gene_names.append(name)
                gene = line[3]
                log2_aFC = line[7]
                totalCount = line[6]
                genes.append(gene)

                gene_data[name][gene] = {'log2_aFC':log2_aFC, 'totalCount':totalCount}

print('writing data')
with open(args.out_prefix+'.geneCounts.txt','w') as outLog2_aFC, open(args.out_prefix+'.totalDepth.txt','w') as outTotalCount:
    for name in gene_names:
        outLog2_aFC.write('\t'+name)
        outTotalCount.write('\t'+name)
    outLog2_aFC.write('\n')
    outTotalCount.write('\n')
    for gene in genes:
        outLog2_aFC.write(gene)
        outTotalCount.write(gene)
        for name in gene_names:
            outLog2_aFC.write('\t'+str(gene_data[name][gene]['log2_aFC']))
            outTotalCount.write('\t'+str(gene_data[name][gene]['totalCount']))
        outLog2_aFC.write('\n')
        outTotalCount.write('\n')
