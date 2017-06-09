chen_data = {}
chen_genes = {}
#1:169581130_G_A rs3917729   ENSG00000000460.11  1.748e-06   0.7945  0.00937304607425    1.290e-05   0.1075  0.1662
with open('/groups/umcg-bios/tmp03/projects/deconvolution/chenReplication/mono/monoSummaryTopGenes.txt') as input_file:
    for line in input_file:
        line = line.strip().split('\t')
        rs_id = line[1]
        gene = line[2].split('.')[0]
        qtl = gene+'_'+rs_id
        chen_data[qtl] = line
        if gene not in chen_genes:
            chen_genes[gene] = []
        chen_genes[gene].append(line)

with open('/groups/umcg-bios/tmp03/projects/deconvolution/chenReplication/mono/replication/intersectGenes.txt') as input_file:
    intersect_genes = input_file.read().split('\n')



gene_in_mono = set([])
gene_in_blood = set([])
gene_in_decon = set([])
with open('BaseCellPerc_merged_filtered.csv') as input_file:
    input_file.readline()
    for line in input_file:
        line = line.strip().split('\t')
        gene = line[0].split('_')[0]
        if gene in chen_genes:
            if gene not in intersect_genes:
                print(gene)
            gene_in_decon.add(gene)
            if float(line[3]) < 0.05:
                gene_in_mono.add(gene)
            if float(line[-1]) < 0.05:
                gene_in_blood.add(gene)

print(len(gene_in_decon))
#print('total genes:',len(chen_genes.keys()))
#print('blood + mono sign:',len(gene_in_mono & gene_in_blood))
#print('blood_only:',len(gene_in_blood -  gene_in_mono))
#print('mono_only:',len(gene_in_mono - gene_in_blood))

