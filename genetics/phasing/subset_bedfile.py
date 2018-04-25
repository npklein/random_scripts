from natsort import natsorted, ns
import argparse

parser = argparse.ArgumentParser('Subset a bedfile to only contain smaller set of region for a chromosome')
parser.add_argument("bed_in", help="bedfile to subset")
parser.add_argument("bed_out",help="bedfile to write")
parser.add_argument("chr", help="chromosome to subset on")
parser.add_argument("start",help="start of chunk", type=int)
parser.add_argument("end",help="end of chunk", type=int)

args = parser.parse_args()

bed_info = {}
with open(args.bed_in) as input_file:
    for line in input_file:
        split_line = line.split('\t')
        bed_info[split_line[0]+':'+split_line[1]+'-'+split_line[2]] = line


with open(args.bed_out,'w') as out:
    for chunk in natsorted(list(bed_info.keys())):
        chr = chunk.split(':')[0]
        start = int(chunk.split(':')[1].split('-')[0])
        end = int(chunk.split('-')[1])
        if chr == args.chr and end > args.start and start < args.end:
            out.write(bed_info[chunk])

