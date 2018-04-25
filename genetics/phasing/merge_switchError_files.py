import glob


def header_index(header):
    d = {}
    for index, h in enumerate(header):
        d[h] = index
    return d

headerToMake = ['chunk','sample','overlapping_snps','snps_only_chunkVCF','snps_only_refVCF',
                'switch','no switch','switchSnps','noSwitchSnps','totalHetsChunk',
                'totalHetsRef','chunkVCF hom alt TO refVCF hom ref','chunkVCF hom ref - refVCF hom alt',
                'chunkVCF hom - refVCF het','chunkVCF het - refVCF hom','chunkVCF het - refVCF het',
                'chunkVCF hom alt TO refVCF hom ref SNPs',
                'chunkVCF hom ref - refVCF hom alt SNPs','chunkVCF hom - refVCF het SNPs',
                'chunkVCF het - refVCF hom SNPs','chunkVCF het - refVCF het SNPs',
                'totalError','haploSize']

# chunk    sample    switch    no switch    switchSnps    noSwitchSnps    totalHetsChunk    
# totalHetsRef    chunkVCF hom alt TO refVCF hom ref    chunkVCF hom ref - refVCF hom alt    
# chunkVCF hom - refVCF het    chunkVCF het - refVCF hom    chunkVCF hom alt TO refVCF hom ref SNPs    
# chunkVCF hom ref - refVCF hom alt SNPs    chunkVCF hom - refVCF het SNPs    chunkVCF het - refVCF hom SNPs    
# totalError    haploSize
with open('switchError.txt','w') as out:
    out.write('\t'.join(headerToMake)+'\tchr\n')
    skip_totalHetsRef = False
    for f in glob.glob('chr*/*switchAndError.txt'):
        chr = f.split('/')[0]
        with open(f) as input_file:
            headerStr = input_file.readline().strip().split('\t')
            headerIndex = header_index(headerStr)
            for line in input_file:
                line = line.strip().split('\t')
                line[0] = line[1].split('.')[0]+'.'+line[0]
                line[1] = line[1].split('.')[1]
                for c in headerToMake:
#                    print(c)
#                    print(c in headerIndex)
                    if c in headerIndex:
#                        if len(line[headerIndex[c]]) > 10:
#                            line[headerIndex[c]] = line[headerIndex[c]][:10]
                        out.write(line[headerIndex[c]])
                    else:
                        out.write('')
                    out.write('\t')
                out.write(chr+'\n')



