#!/usr/bin/python

import sys
import gzip
import argparse

parser = argparse.ArgumentParser(description='Remove additional contigs (GL*, LRG* etc). Output is streamed, not written to file')
parser.add_argument('in_gtf', help='Input gtf to remove contigs from')

args = parser.parse_args()

fname = args.in_gtf
if fname.endswith(".gz"):
    f = gzip.open(fname)
else:
    f = open(fname)

for l in f:
    l = l.strip()
    spl = l.split("\t")
    if (spl[0].isdigit()) or (spl[0] in ["X", "Y", "MT", "M"]):
        print (l)
f.close()
