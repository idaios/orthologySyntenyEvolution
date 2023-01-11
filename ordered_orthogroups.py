#!/usr/bin/python3
import argparse
import re
import sys
import glob

# change this if needed
fadir="fa_files" 

parser = argparse.ArgumentParser()
parser.add_argument('-gl', '--genome_list', help="Provide a list of genomes that you want to analyze. This should be a file with a single filename in each line that represents a single proteome (or genome). A PAIR of proteomes should be provided here", required=True)
parser.add_argument('-og', '--orthology_results', help="Provide the results of Orthofinder which refer to the orthogroups", required=True)
parser.add_argument('-n', '--neighbors', help="The number of neighboring genes to assess the synteny proportion")

# Parse and print the results
args = parser.parse_args()

orthfile = args.orthology_results

# set as a delimiter either the tab or the comma. Here, I don't care
# currently to know at which species each gene belongs to. Thus, I
# will just keep a track of all genes that belong to an orthogroup
betweenSpeciesDelimiter = re.compile("[^\t]+")
inSpeciesDelimiter = re.compile("[^,]+")

data = {}
orthogroupID = ""
orthoIndex=0
orthoData=[]
orthoCount=[]
oneToManyIndexes=[]
tmpGenes=[None]*2
with open(orthfile, 'r') as orth:
    while True:
        ln = orth.readline()
        if not ln:
            break
        ln = ln.strip()
        tokens = betweenSpeciesDelimiter.findall(ln)
        index = 0
        
        for t in tokens:
            if index == 0:
                orthogroupID = t
                index=1
                continue

            genes = inSpeciesDelimiter.findall(t)
            tmpGenes[index-1] = genes
            for g in genes:
                data[g] = orthogroupID
                
            index+=1
            
        ## append as a list
        orthoCount.append( [len(tmpGenes[i]) for i in range(index-1) ] )
        ## append as a list. Not sure if this is necessary, but let's
        ## try it
        orthoData.append( [ tmpGenes[i] for i in range(index-1) ] )

        if max(orthoCount[-1]) > 1 and min(orthoCount[-1]) == 1:
            ## one to many relation
            oneToManyIndexes.append(len(orthoCount)-1)


## read the genomes
with open(args.genome_list, 'r') as inplist:
    while True:
        f = inplist.readline()
        f = f.strip()
        if not f:
            break
        
        path=fadir+"/"+f+"*"
        fastaFile = glob.glob(path)[0]

        with open(fastaFile, 'r') as fa:
            for ln in fa:
                ln.strip()
                

# for i in oneToManyIndexes:
#     print(orthoData[i])

    

