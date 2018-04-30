import argparse
import gzip

parser = argparse.ArgumentParser(description='Extract Chr, Pos and dosage from vcf.gz file.')
parser.add_argument('in_vcf',
                    help='vcf to make dosage table from (has to be gzipped)')
parser.add_argument('out_table', 
                    help='output file name')

args = parser.parse_args()
with gzip.open(args.in_vcf,'rb') as input_file, open(args.out_table,'w') as out:
    for line in input_file:
        line = line.decode('utf-8')
        if line.startswith('#CHR'):
            _line = line.strip().lstrip('#').split('\t')
            out.write(_line[0]+'\t'+_line[1]+'\t'+_line[3]+'\t'+_line[4]+'\t'+'\t'.join(_line[9:])+'\n')
        if line.startswith('#'):
            continue
        line = line.strip().split('\t')
        out.write(line[0]+'\t'+line[1]+'\t'+line[3]+'\t'+line[4]+'\t')
        for element in line[9:]:
            out.write('\t')
            element = element.split(':')[0]
            if element == '0/1':
                out.write('1')
            elif element == '0/0':
                out.write('0')
            elif element == '1/1':
                out.write('2')
            elif element == './.':
                out.write('.')
            else:
                raise RuntimeError('should have covered all posibilities, check code')
        out.write('\n')
