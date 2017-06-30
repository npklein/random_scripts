import argparse
import gzip
import sys

parser = argparse.ArgumentParser(description='Select variants based on position.')
parser.add_argument('V',help='input vcf')
parser.add_argument('L',help='interval list')
parser.add_argument('o',help='output vcf')

args = parser.parse_args()

intervals = []
with open(args.L) as input_file:
    for interval in input_file:
        intervals.append(interval.strip())
with gzip.open(args.V,'rb') as input_file, open(args.o,'w') as output_file:
    x = 0
    for line in input_file:
        line = line.decode('utf-8')
        if line.startswith('#'):
            output_file.write(line)
        else:
            line_spl = line.split('\t')
            if line_spl[0]+':'+line_spl[1]+'-'+line_spl[1] in intervals:
                output_file.write(line)
                x += 1
                if x % 5000 == 0:
                    print(str(x)+'/'+str(len(intervals)))
                    sys.stdout.flush()
print(str(x)+'/'+str(len(intervals)))
sys.stdout.flush()

