#!/usr/bin/python3
import argparse
import re
import sys
import glob
import numpy as np

# change this if needed
fadir="fa_files" 

n = 10
parser = argparse.ArgumentParser()
parser.add_argument('-gl', '--genome_list', help="Provide a list of genomes that you want to analyze. This should be a file with a single filename in each line that represents a single proteome (or genome). A PAIR of proteomes should be provided here", required=True)
parser.add_argument('-og', '--orthology_results', help="Provide the results of Orthofinder which refer to the orthogroups", required=True)
parser.add_argument('-n', '--neighbors', help="The number of neighboring genes to assess the synteny proportion")


### nw function #################################
def nw(x, y, match = 1, mismatch = 1, gap = 2):
    score = 0
    nx = len(x)
    ny = len(y)
    # Optimal score at each possible pair of characters.
    F = np.zeros((nx + 1, ny + 1))
    F[:,0] = np.linspace(0, -nx * gap, nx + 1)
    F[0,:] = np.linspace(0, -ny * gap, ny + 1)
    # Pointers to trace through an optimal aligment.
    P = np.zeros((nx + 1, ny + 1))
    P[:,0] = 3
    P[0,:] = 4
    # Temporary scores.
    t = np.zeros(3)
    for i in range(nx):
        for j in range(ny):
            if x[i] == y[j]:
                t[0] = F[i,j] + match
            else:
                t[0] = F[i,j] - mismatch
            t[1] = F[i,j+1] - gap
            t[2] = F[i+1,j] - gap
            tmax = np.max(t)
            F[i+1,j+1] = tmax
            if t[0] == tmax:
                P[i+1,j+1] += 2
            if t[1] == tmax:
                P[i+1,j+1] += 3
            if t[2] == tmax:
                P[i+1,j+1] += 4
    score = tmax
    # Trace through an optimal alignment.
    i = nx
    j = ny
    rx = []
    ry = []
    while i > 0 or j > 0:
        if P[i,j] in [2, 5, 6, 9]:
            rx.append(x[i-1])
            ry.append(y[j-1])
            i -= 1
            j -= 1
        elif P[i,j] in [3, 5, 7, 9]:
            rx.append(x[i-1])
            ry.append('-')
            i -= 1
        elif P[i,j] in [4, 6, 7, 9]:
            rx.append('-')
            ry.append(y[j-1])
            j -= 1
    # Reverse the strings.
    rx = ''.join(rx)[::-1]
    ry = ''.join(ry)[::-1]
    return score #'\n'.join([rx, ry])


#################################################



# Parse and print the results
args = parser.parse_args()

orthfile = args.orthology_results

# set as a delimiter either the tab or the comma. Here, I don't care
# currently to know at which species each gene belongs to. Thus, I
# will just keep a track of all genes that belong to an orthogroup
betweenSpeciesDelimiter = re.compile("[^\t]+")
inSpeciesDelimiter = re.compile("[^,]+")
## the specification of a fasta name in order to parse the names of sequences
fastaSpec = re.compile(r"^>([^ ]+)")

data = {}
orthogroupID = ""
orthoIndex=0
orthoData=[]
orthoCount=[]
oneToManyIndexes=[]
tmpGenes=[None]*2
geneOrder=[]
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
                g = g.strip()
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



dictGenes = {}
vGenes = []
orgDictGenes = []
orgVGenes = []
seq={}
## read the genomes

with open(args.genome_list, 'r') as inplist:
    while True:
        f = inplist.readline()
        f = f.strip()
        if not f:
            break
        
        path=fadir+"/"+f+"*"
        fastaFile = glob.glob(path)[0]
        curIndex = 0
        vGenes = []
        aGene = ""
        with open(fastaFile, 'r') as fa:
            for ln in fa:
                ln = ln.strip()
                m = fastaSpec.match(ln)
                if m:
                    if aGene != "":
                        sys.stderr.write("Error parsing sequences\n")
                        sys.exit()
                    aGene = m.group(1)
                    dictGenes[aGene] = curIndex
                    vGenes.append(aGene)
                    curIndex+=1
                else:
                    if aGene == "":
                        sys.stderr.write("Error 2 parsing sequences\n")
                        sys.exit()
                    seq[aGene] = ln
                    aGene = ""
                        
        orgDictGenes.append(dictGenes)
        orgVGenes.append(vGenes)

##print(orgDictGenes[0]['ENSPTRP00000071475.1'])

## Now, process the results 1-many initially since, I guess they are
## the most interesting
for i in range(len(oneToManyIndexes)):
    ind = oneToManyIndexes[i]
    oneIndex = 0
    
    for j in range(len(orthoCount[ind])):
        if orthoCount[ind][j] == 1:
            oneIndex = j
            break
    ## at this point we know the reference gene
    refGene = orthoData[ind][oneIndex]
    #print(refGene)
    indexRefGene = orgDictGenes[oneIndex][refGene[0]]
    minIndexRefGene = max(0, indexRefGene - n)
    maxIndexRefGene = min(indexRefGene + n, len(orgVGenes[oneIndex])-1)
    #print(str(refGene) + " " + str(indexRefGene)+" "+str(minIndexRefGene)+" "+str(maxIndexRefGene))

    neighborOrthRefGene = []
    for j in range(minIndexRefGene, maxIndexRefGene+1, 1):
        #if j != indexRefGene:
        key = orgVGenes[oneIndex][j]
        if key in data:
            orthogroup = data[ key ]
            if j == indexRefGene:
                orthogroup = "*"+orthogroup+"*"
            neighborOrthRefGene.append(orthogroup)

    ## for all other organisms but the reference (the reference is the organism that has one copy)
    for j in range(len(orthoCount[ind])):
        ## skip the organism that is used as reference
        if j == oneIndex:
            continue
        ## now, j is the index for another organism
        ## This is the list of orthogroup partners of the reference gene for
        ## the organism j
        listOfOrtho = orthoData[ind][j]
        
        #print("The list of orthologues for the "+refGene[0]+"("+data[refGene[0]]+") is: ", end=" ")
        #print(listOfOrtho)
        for orthogene in listOfOrtho:
            ## HERE
            orthogene = orthogene.strip()
            indexOrthGene = orgDictGenes[j][orthogene]
            minIndexOrthGene = max(0, indexOrthGene - n)
            maxIndexOrthGene = min(indexOrthGene + n, len(orgVGenes[j])-1)
            #print(str(orthogene) + " " + str(indexOrthGene)+" "+str(minIndexOrthGene)+" "+str(maxIndexOrthGene))

            neighborOrthOrthGene = []
            for jj in range(minIndexOrthGene, maxIndexOrthGene+1, 1):
                #if jj != indexOrthGene:
                key = orgVGenes[j][jj]
                if jj == indexOrthGene and key not in data:
                    sys.stderr.write("Cannot be that the gene --"+key+"-- is not in the data\n")
                    print(data[key])
                    sys.exit()
                                     
                if key in data:
                    orthogroup = data[ key ]
                    if jj == indexOrthGene:
                        orthogroup = "**"+orthogroup+"**"
                    neighborOrthOrthGene.append(orthogroup)
            #print ("indexOrth:"+str(indexOrthGene)+" min:"+str(minIndexOrthGene)+" max:"+str(maxIndexOrthGene))
            common = list(set(neighborOrthRefGene).intersection(neighborOrthOrthGene))
            percentageCommon = len(common)/len(set(neighborOrthOrthGene+neighborOrthRefGene))
            alScore = nw(seq[refGene[0]], seq[orthogene])/(len(seq[refGene[0]])+len(seq[orthogene]))
            print(refGene[0]+" "+orthogene+" "+str(percentageCommon)+" "+str(alScore)+" "+str(len(seq[refGene[0]]))+" "+str(len(seq[orthogene])))






            

            
            
    #print(len(neighborOrthRefGene))
    #print(neighborOrthRefGene)

# for i in oneToManyIndexes:
#     print(orthoData[i])

    

