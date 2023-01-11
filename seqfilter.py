#!/usr/bin/python3
import argparse
import re
import sys

parser = argparse.ArgumentParser(description = 'Filters a fasta file removing sequences that do not fullfil a length criterion')

parser.add_argument('-i','--input', help='The input file')
parser.add_argument('-l', '--length', type=int, help='The minimum length')
parser.add_argument('-m', '--methionine', help="The first aminoacid should be methionine", action="store_true")
parser.add_argument('-U', '--unique_longer', help="Get only a single pep per gene, the longer one", action="store_true")
parser.add_argument('-r', '--unique_random', help="Get onle a single pep per gene, a random one", action="store_true")

args = parser.parse_args()

if args.unique_longer and args.unique_random:
    sys.exit("--unique_longer and --unique_random cannot be used simultaneously")

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
        ## let's suppose we are reading seq name i. Then match is
        ## True, 
        
        if m:
            totalSeqs += 1
            skipped = False
            if len(seq) >= args.length  and ( (args.methionine and seq[0] == "M") or (args.methionine == False) ):
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
               # elif args.unique_longer and not g:
               #     seq = ""
                    #sys.stderr.write("WARNING:\n"+ln+" NO CHROMOSOME INFO\n")
                elif args.unique_longer == False and args.unique_random == False:
                    print(">"+name+"\n"+seq+'\n', end='')
        
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

if seq != "":
    if len(seq) >= args.length  and ( (args.methionine and seq[0] == "M") or (args.methionine == False) ):
        if args.unique_longer:
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
        elif args.unique_longer == False and args.unique_random == False:
            print(">"+name+"\n"+seq+'\n', end='')


sortedkeys = sorted(dataposition.keys(), key = lambda x: dataposition[x], reverse=False)

sys.stderr.write("\n"+args.input+": Total sequences: "+str(totalSeqs)+"\n"+str(skippingSeqs)+" were skipped. Probably because of scaffold\n")

if args.unique_longer:
    #print("TRUE")
    for k in sortedkeys:
        #print(k)
        #print(k+":"+str(dataposition[k])+str(dataseqname[k]))
        print(">"+str(dataseqname[k]))
        print(str(data[k]))
