#!/usr/bin/python3
import argparse
import re
import sys

parser = argparse.ArgumentParser()
parser.add_argument('-gl', '--genome_list', help="Provide a list of genomes that you want to analyze. This should be a file with a single filename in each line that represents a single proteome (or genome). A PAIR of proteomes should be provided here")
parser.add_argument('-og', '--orthology_results', help="Provide the results of Orthofinder which refer to the orthogroups")

# Parse and print the results
args = parser.parse_args()

orthfile = args.orthology_results

delimiter = re.compile("[^\t]+")
indelimiter = re.compile("[^,]+")

data = {}

with open(orthfile, 'r') as orth:
    while True:
        ln = orth.readline()
        ln = ln.strip()
        tokens = delimiter.findall(ln)
        index = 0
        for t in tokens:
            if index == 0:
                key = t
                index = 1
                continue
            g = indelimiter.findall(t)
            if index == 1:
                sp1 = g
            elif index == 2:
                sp2 = g
            index+=1
        if key in data:
            sys.stderr.write("Error key is already present")
        data[key] = (len(sp1), len(sp2))

        if not ln:
            break

cnt = 0


for k in data.keys():
    if data[k][0] == 1 and data[k][1] == 1:
        print(str(k)+":"+str(data[k]))
        cnt+=1

print(cnt)
print(len(data))
