#!/usr/bin/python3
import argparse
import re
import sys
import numpy as np

parser = argparse.ArgumentParser(description = 'Filters a fasta file removing sequences that do not fullfil a length criterion')

parser.add_argument('-i','--input', help='The input file')
parser.add_argument('-U', '--unique_longer', help="Get only a single pep per gene, the longer one", action="store_true", default=True)
parser.add_argument('-header', '--header', help="Print the header of percentiles", action='store_true', default=False)

args = parser.parse_args()
seqre = re.compile(r'^>([^\[]+)')
seqreGene = re.compile(r'^>.*(?:chromosome|primary_assembly):[^:]+:([^:]+):([^:]+):([^:]+).*gene:([^ ]+)')

name = ""
skipname = ""
seq = ""
data = {}
datalength = {}
dataseqname = {}
dataposition = {}
g = False
skippingSeqs = 0
totalSeqs = 0
skipped = False
with open(args.input, 'r') as f:
    while True:
        ln = f.readline() # why not using readlines?
        ln = ln.strip()
        m = seqre.search(ln)
        if m:
            totalSeqs += 1
            skipped = False
            if args.unique_longer and g:
                if geneName not in data:
                    data[geneName] = seq
                    datalength[geneName] = len(seq)
                    dataseqname[geneName] = name
                    dataposition[geneName] = geneInfo
                elif datalength[geneName] < len(seq):
                    data[geneName] = seq
                    datalength[geneName] = len(seq)
                    dataseqname[geneName] = name
                    dataposition[geneName] = geneInfo
                                    
            name = m.group(1)
            if args.unique_longer:
                g = seqreGene.search(ln)
                if g:
                    geneName = g.group(4)
                    geneInfo = (g.group(1), int(g.group(2)), int(g.group(3)))
                else:
                    # skip
                    skipname = name
                    name = ""
                
            seq = ""
        else: 
            if (seq == "") and (name != ""):
                seq = ln
            elif(seq != "") and (name != ""):
                seq += ln
            elif name == "" and not skipped:
                skippingSeqs += 1
                skipped = True
                #sys.stderr.write("name: "+name+"\nseq: "+seq+"\n")
                #sys.stderr.write(args.input+": Skipping sequence: "+
                #skipname + "\n")
                        
        if not ln:
            break
lengths = [datalength[k] for k in datalength.keys()]
perc = [0, 5, 10, 50, 90, 95, 100]
percentiles = {perc[i]:np.percentile(lengths, perc)[i] for i in range(len(perc))}
if args.header:
    for p in perc:
        print(str(p)+"\t")
print(args.input+"\t", end="")
for v in percentiles.values():
    print(int(v), "\t", end="")
print()
    

