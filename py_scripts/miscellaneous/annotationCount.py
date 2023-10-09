import sys
import os
import pandas as pd

df = pd.read_csv('/Users/mmeynadier/Documents/PhD/species/Clytia/raw/ChGenes.gtf',sep='\t',header=None)

listDf = df[8].to_list() 
geneList = [] ; transcriptList = []
for i in listDf:
    firstParse = i.split(';')[1]
    geneList.append(i.split(';')[0])
    transcriptList.append(firstParse.split(';')[0])
print(len(geneList)) ; print(len(transcriptList))
geneList = set(geneList)
transcriptList = set(transcriptList)
print(len(geneList)) ; print(len(transcriptList))