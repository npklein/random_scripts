import glob
from pyensembl import EnsemblRelease

data = EnsemblRelease(75) # will take a while to get the genome data installed in your system

def list_index(l):
    d = {}
    for index, element in enumerate(l):
        d[element] = index
    return d

snps = set([])
for chr in glob.glob('chr*'):
    print(chr)
    for f in glob.glob(chr+'/*.switchAndError.txt'):
        with open(f) as input_file:
            header_index = list_index(input_file.readline().strip().split('\t'))
            for line in input_file:
                line = line.strip().split('\t')
                for snp in line[header_index['overlapping_snps']].split(','):
                    if len(snp.strip()) == 0:
                        continue
                    snps.add(snp)
                for snp in line[header_index['snps_only_chunkVCF']].split(','):
                    if len(snp.strip()) == 0:
                        continue
                    snps.add(snp)
                for snp in line[header_index['snps_only_refVCF']].split(','):
                    if len(snp.strip()) == 0:
                        continue
                    snps.add(snp)

with open('snp_gene_locations.txt','w') as out:
    number_of_snps = len(snps)
    for index, snp in enumerate(snps):
        if index % 5000 == 0:
            print(str(index)+'/'+str(number_of_snps))
        gene_names = data.gene_names_at_locus(contig=chr.replace("chr",""), position=int(snp.split('_')[1]))
        gene_ids = []
        for gene_name in gene_names:
            for gene_id in data.gene_ids_of_gene_name(gene_name):
                gene_ids.append(gene_id)
        out.write(snp+'\t'+','.join(gene_names)+'\t'+','.join(gene_ids)+'\n')
