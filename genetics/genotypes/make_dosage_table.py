import argparse
import gzip

parser = argparse.ArgumentParser(description='Extract Chr, Pos and dosage from vcf.gz file.')
parser.add_argument('in_vcf',
                    help='vcf to make dosage table from (has to be gzipped)')
parser.add_argument('out_table', 
                    help='output file name')

args = parser.parse_args()
with gzip.open(args.in_vcf,'rb') as input_file, open(args.out_table,'w') as out:
    vcf_length = None
    for line in input_file:
        line = line.decode('utf-8')
        if line.startswith('#CHR'):
            _line = line.strip().lstrip('#').split('\t')
            # vcf length is the number of column headers in the vcf file. Appended are bedfile columns, but these do not have headers
            vcf_length = len(_line)
            out.write(_line[0]+'\t'+_line[1]+'\t'+_line[3]+'\t'+_line[4]+'\tgene\t'+'\t'.join(_line[9:])+'\n')
        if line.startswith('#'):
            continue
        line = line.strip().split('\t')
        out.write(line[0]+'\t'+line[1]+'\t'+line[3]+'\t'+line[4])
        # the VCF file is merged with a bed file, of which the 4th column is the ensemble ID. This means that
        # vcf length (determined by the columns with a header)+3 gives the ensembl ID (+3 cause index starts at 0)
        out.write('\t'+line[vcf_length+3])
        for index, element in enumerate(line[9:]):
            if index + 9 >= vcf_length:
                break
            out.write('\t')
            element = element.split(':')[0]
            element = element.replace('|','/')
            if element == '0/1' or element == '1/0':
                out.write('1')
            elif element == '0/0':
                out.write('0')
            elif element == '1/1':
                out.write('2')
            elif element == './.':
                out.write('.')
            else:
                raise RuntimeError('should have covered all posibilities, check code. element was: '+element)
        out.write('\n')
