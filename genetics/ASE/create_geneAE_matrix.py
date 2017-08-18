import glob
import argparse
import os
import sys
parser = argparse.ArgumentParser(description='Make matrix (genes x samples) of the gene counts and total depth for all results in given dir + subdirs. IMPORTANT: sample names are taken by splitting bam name on . and taking first [0] element.')
parser.add_argument('geneAE_rootdir', help='Rootdir where all gene counts can be found. Will walk through subdirs and use all *.txt files')
parser.add_argument('out_prefix', help='Write output to <out_prefix>.geneCounts.txt, <out_prefix>.totalDepth.txt, and <out_prefix>.alleleCounts.txt')
parser.add_argument('--verbose', action='store_true',help='If set, flush_print out which files are being used')
parser.add_argument('--removeNoCounts', action='store_true', help='Remove genes for which totalDepth is 0 for all samples')
parser.add_argument('--chr', help='Which chromosome to make table for')

args = parser.parse_args()

def flush_print(message):
    print(message)
    sys.stdout.flush()


flush_print('Combining gene count files from '+args.geneAE_rootdir)
if args.chr:
    flush_print('for choromosome '+args.chr)
flush_print('Writing results to '+args.out_prefix)

flush_print('Start combining...')

gene_data = {}
sample_names = []
genes = []
flush_print('reading data')
#contig    start    stop    name    aCount    bCount    totalCount    log2_aFC    n_variants    variants    gw_phased    bam
x = 0
for subdir, dirs, files in os.walk(args.geneAE_rootdir):
    flush_print(subdir)
    flush_print(str(len(files))+' files in this subdir')
    for file in files:
        if not file.endswith('.txt'):
            continue
        x += 1
        if x % 1000 == 0:
            flush_print(str(x)+' files read')
        geneAE_file = subdir+'/'+file
        if args.verbose:
            flush_print(geneAE_file)
        with open(geneAE_file) as input_file:
            input_file.readline()
            for line in input_file:                
                line = line.strip().split('\t')
                
                if args.chr and line[0] != args.chr:
                    continue
                name = line[-1].split('.')[0]
                if not name in gene_data:
                    gene_data[name] = {}
                    sample_names.append(name)
                gene = line[3]
                log2_aFC = line[7]
                totalCount = line[6]
                aCount = line[4]
                bCount = line[5]
                if not gene in genes:
                    genes.append(gene)
                gene_data[name][gene] = {'log2_aFC':log2_aFC, 'totalCount':totalCount,
                                         'aCount':aCount,'bCount':bCount}

flush_print('Filling in values for samples that do not include all genes')
if args.removeNoCounts:
    flush_print('removing genes where all samples have 0 depth')
for gene in list(genes):
    noZeros = False
    for name in sample_names:
        if gene not in gene_data[name]:
            gene_data[name][gene] = {'log2_aFC':'inf', 'totalCount':0,
                                         'aCount':0,'bCount':0}
        if args.removeNoCounts:
            if int(gene_data[name][gene]['totalCount']) > 0:
                noZeros = True
    if not noZeros:
        flush_print('Gene has 0 totalCount for all samples, removing '+gene)
        for name in sample_names:
            del(gene_data[name][gene])
        genes.remove(gene)

flush_print('writing data')
with open(args.out_prefix+'.geneCounts.txt','w') as outLog2_aFC, open(args.out_prefix+'.totalDepth.txt','w') as outTotalCount:
    with open(args.out_prefix+'.alleleCounts.txt','w') as outAlleleCounts:
        for name in sample_names:
            outLog2_aFC.write('\t'+name)
            outTotalCount.write('\t'+name)
            outAlleleCounts.write('\t'+name+'_aCount')
            outAlleleCounts.write('\t'+name+'_bCount')
        outLog2_aFC.write('\n')
        outTotalCount.write('\n')
        outAlleleCounts.write('\n')
        for gene in sorted(genes):
            outLog2_aFC.write(gene)
            outTotalCount.write(gene)
            outAlleleCounts.write(gene)
            for name in sample_names:
                outLog2_aFC.write('\t'+str(gene_data[name][gene]['log2_aFC']))
                outTotalCount.write('\t'+str(gene_data[name][gene]['totalCount']))
                outAlleleCounts.write('\t'+str(gene_data[name][gene]['aCount']))
                outAlleleCounts.write('\t'+str(gene_data[name][gene]['bCount']))
            outLog2_aFC.write('\n')
            outTotalCount.write('\n')
            outAlleleCounts.write('\n')
