import glob
import argparse
import os
import sys
parser = argparse.ArgumentParser(description='Merge matrices from phASER geneAE together')
parser.add_argument('matrix_rootdir', help='Rootdir where all matrices can be found. Will walk through subdirs and use all *.geneCounts.txt, *.totalDepth.txt, and *.alleleCounts.txt files')
parser.add_argument('out_prefix', help='Write output to <out_prefix>.geneCounts.txt, <out_prefix>.totalDepth.txt, and <out_prefix>.alleleCounts.txt')

args = parser.parse_args()

def flush_print(message):
    print(message)
    sys.stdout.flush()



flush_print('Combining matrix files from '+args.matrix_rootdir)
flush_print('Writing results to '+args.out_prefix)

flush_print('Start combining...')

gene_counts = {}
total_depth = {}
allele_counts = {}
sample_names = set([])
sample_names_alleleCounts = set([])
genes = set([])
flush_print('reading data')
x = 0
for subdir, dirs, files in os.walk(args.matrix_rootdir):
    flush_print(subdir)
    for file in files:
        if not file.endswith('.txt'):
            continue
        matrix_file = subdir+'/'+file
        with open(matrix_file) as input_file:
            header = input_file.readline().strip()
            if len(header) == 0:
                continue
            samples = header.split('\t')
            sample_index = {}
            for index,sample in enumerate(samples):
                if file.endswith('alleleCounts.txt'):
                    if sample not in allele_counts:
                        allele_counts[sample] = {}
                        sample_names_alleleCounts.add(sample)
                elif sample not in gene_counts:
                    gene_counts[sample] = {}
                    total_depth[sample] = {}
                    sample_names.add(sample)
                sample_index[index] = sample

            for line in input_file:
                line = line.strip().split('\t')
                gene = line[0]
                genes.add(gene)
                for index, element in enumerate(line[1:]):
                    if file.endswith('alleleCounts.txt'):
                        allele_counts[sample_index[index]][gene] =  element                        
                    elif file.endswith('geneCounts.txt'):
#                        if gene in gene_counts[sample_index[index]]:
#                            raise RuntimeError('Each gene should only occur once per sample, '+gene+' happened twice for '+sample)
                        print(gene)
                        gene_counts[sample_index[index]][gene] =  element
                    elif file.endswith('totalDepth.txt'):
                        total_depth[sample_index[index]][gene] = element
                    elif file.endswith('snps.txt'):
                        #not implemented yet
                        continue
                    else:
                        raise RuntimeError('Wrong ending of input file for '+file)


flush_print('Filling in values for samples that do not include all genes')
for gene in list(genes):
    for name in sample_names:
        if gene not in gene_counts[name]:
            gene_counts[name][gene] = 'inf'
            total_depth[name][gene] = 0
    for name in sample_names_alleleCounts:
        if gene not in allele_counts[name]:
            allele_counts[name][gene] = 0



flush_print('writing data')
with open(args.out_prefix+'.geneCounts.txt','w') as outLog2_aFC, open(args.out_prefix+'.totalDepth.txt','w') as outTotalCount:
    with open(args.out_prefix+'.alleleCounts.txt','w') as outAlleleCounts:
        for name in sorted(sample_names):
            outLog2_aFC.write('\t'+name)
            outTotalCount.write('\t'+name)
        for name in sorted(sample_names_alleleCounts):
            outAlleleCounts.write('\t'+name)
        outLog2_aFC.write('\n')
        outTotalCount.write('\n')
        outAlleleCounts.write('\n')
        for gene in sorted(genes):
            outLog2_aFC.write(gene)
            outTotalCount.write(gene)
            outAlleleCounts.write(gene)
            for name in sorted(sample_names):
                outLog2_aFC.write('\t'+str(gene_counts[name][gene]))
                outTotalCount.write('\t'+str(total_depth[name][gene]))
            for name in sorted(sample_names_alleleCounts):
                outAlleleCounts.write('\t'+str(allele_counts[name][gene]))
            outLog2_aFC.write('\n')
            outTotalCount.write('\n')
            outAlleleCounts.write('\n')
