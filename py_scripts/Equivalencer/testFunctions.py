import pandas as pd
import statistics
import math
import numpy as np
import pickle
from scipy.stats import shapiro 
import statsmodels.api as sm
from functools import reduce
import scipy.stats as stats
import seaborn as sns
import random
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from collections import defaultdict
import os
import sys
import re
from pheatmap import pheatmap
from functools import reduce
from Equivalencer_basic import *
from Equivalencer_stat import *

"""
def readPfamFile(path,file):
    with open(path+file,'r') as f:
        l = f.readline()
        while l != "":
            if l[0] == "#":
                print("test")
                l = f.readline()
                print(l)
"""
#readPfamFile('/Users/mmeynadier/Documents/PhD/species/crossSpecies/RBH_x_markers_identity/','Clytia_Nematostella_Xenia_automatic_cnidocytes.pfam')

#table = [[2, 1], [184, 316]]

#a=[601, 218, 737, 247, 694, 467, 433, 203, 70, 186, 592, 552, 627, 369, 774, 570, 227, 317, 148, 215, 486, 386, 63, 169, 42, 311, 158, 220, 76, 93, 628, 41]
#b=[399, 394, 263, 300, 306, 532, 452, 196, 65, 317, 408, 448, 373, 230, 226, 409, 309, 254, 249, 262, 513, 325, 186, 185, 84, 437, 157, 150, 162, 141, 372, 93]
#c=[6810, 7193, 6674, 7164, 6717, 6944, 6978, 7208, 7341, 7225, 6819, 6859, 6784, 7042, 6637, 6841, 7184, 7094, 7263, 7196, 6925, 7025, 7348, 7242, 7369, 7100, 7253, 7191, 7335, 7318, 6783, 7370]
#d=[18917, 18922, 19053, 19016, 19010, 18784, 18864, 19120, 19251, 18999, 18908, 18868, 18943, 19086, 19090, 18907, 19007, 19062, 19067, 19054, 18803, 18991, 19130, 19131, 19232, 18879, 19159, 19166, 19154, 19175, 18944, 19223]

pvalDict = {}

"""
for i in range(len(a)):
    table = [[a[i],b[i]],[c[i],d[i]]]
    pval = stats.fisher_exact(table,alternative='greater')
    pvalDict[i] = pval[1]
print(pvalDict)    
""" 

def combine_hex_values(d):
    """
    Description
    -----------
    

    Parameters
    ----------
    d

    Returns
    ----------
    zpad(hex(red)[2:]) + zpad(hex(green)[2:]) + zpad(hex(blue)[2:])
    """
    d_items = sorted(d.items())
    tot_weight = sum(d.values())
    red = int(sum([int(k[:2], 16)*v for k, v in d_items])/tot_weight)
    green = int(sum([int(k[2:4], 16)*v for k, v in d_items])/tot_weight)
    blue = int(sum([int(k[4:6], 16)*v for k, v in d_items])/tot_weight)
    zpad = lambda x: x if len(x)==2 else '0' + x
    return zpad(hex(red)[2:]) + zpad(hex(green)[2:]) + zpad(hex(blue)[2:])

#print(stats.fisher_exact(table))
#print(stats.fisher_exact(table,alternative='greater'))
#print(stats.fisher_exact(table,alternative='less'))

#print(stats.chi2_contingency(table,correction='True'))
#print(stats.chi2_contingency(table,correction='False'))


#print(len(nematostellaContingency[0]))

"""
print("Clytia")
for i in range(len(clytiaContingency)):
    for j in range(len(clytiaContingency[i])):
        for k in range(len(clytiaContingency[i][j])):
            if j == 16:
                if k == 9:
                    print(clytiaContingency[i][j][k]) 

print("Nematostella")
for i in range(len(nematostellaContingency)):
    for j in range(len(nematostellaContingency[i])):
        for k in range(len(nematostellaContingency[i][j])):
            if j == 0:
                if k == 0:
                    print(nematostellaContingency[i][j][k])


print(stats.fisher_exact([[105,431],[206,25705]],alternative='greater'))
print(stats.fisher_exact([[15, 498], [586, 20397]],alternative='greater'))  


with open('results/stat/Clytia_adult_Clytia_larva_markers_matrix_list_1000g_genomeBased_Clytia.csv','rb') as f:
    matrixGenome = pickle.load(f)

with open('results/stat/Clytia_adult_Nematostella_adult_markers_matrix_list_1000g_orthopairsBased.csv','rb') as f:
    matrixOrtho = pickle.load(f)

print(matrixGenome)
print("\n\n\n")
#print(matrixOrtho)

print("Genome")
print(stats.fisher_exact([[68, 932], [667, 25060]],alternative='greater'))  

print("\n")
print("Ortho")
print(stats.fisher_exact([[8, 593], [45, 6765]],alternative='greater'))  

listSp = ["Clytia","Nematostella"]
listLifestage=["adult","adult"]
test="Fisher"
contMethod=["genomeBased","orthopairsBased"]

genomeClytia = pd.read_csv("results/stat/"+listSp[0]+"_"+listLifestage[0]+"_"+listSp[1]+"_"+listLifestage[1]+"_"+test+"_test_1000g_"+contMethod[0]+"_tippett.csv").drop(columns=['Unnamed: 0'])
orthoClytia = pd.read_csv("results/stat/"+listSp[0]+"_"+listLifestage[0]+"_"+listSp[1]+"_"+listLifestage[1]+"_"+test+"_test_1000g_"+contMethod[1]+".csv").drop(columns=['Unnamed: 0'])

print(genomeClytia)
print("ORTHO")
print(orthoClytia)
"""

def unique(list1):
    res = reduce(lambda re, x: re+[x] if x not in re else re, list1, [])
    return res





def manipData(spList,orthoMode,geneSource): 
    try:
        with open('results/basic/'+geneSource+'/'+orthoMode+'/'+spList[0]+'_'+spList[1]+'_orthologs.csv', "rb") as f:
            RBHtupleList = pickle.load(f)
    except FileNotFoundError:
        with open('results/basic/'+geneSource+'/'+orthoMode+'/'+spList[1]+'_'+spList[0]+'_orthologs.csv', "rb") as f:
            RBHtupleList = pickle.load(f)
    return RBHtupleList



def genesCounting(RBHtupleList,markersDictSpList,orthoMode,contMethod,species,spList):
    sp1MarkersDict = markersDictSpList[0]
    sp2MarkersDict = markersDictSpList[1]
    caseA = []
    caseB = []
    caseC = []
    caseD = []
    if RBHtupleList == "Decorator":
        allConcatGenes = []
        for k1,v1 in sp1MarkersDict.items():
            for k2,v2 in sp2MarkersDict.items():
                concatGenes = list(set(v1 + v2))
            allConcatGenes = allConcatGenes + concatGenes
        allConcatGenes = list(set(allConcatGenes))
    if contMethod == "genomeBased":
        proteom = 'input/proteins/longest_isoform_'+species+'Proteins.fasta'
        geneList = longestIsoforms(proteom)
        orthologsList = []
        if species == spList[0]:
            sp1MarkersDict = markersDictSpList[0]
            sp2MarkersDict = markersDictSpList[1]
        elif species == spList[1]:
            sp1MarkersDict = markersDictSpList[1]
            sp2MarkersDict = markersDictSpList[0]
        for i in RBHtupleList:
            orthologsList = np.append(orthologsList,i[0]) ; orthologsList = np.append(orthologsList,i[1])
    for k1, v1 in sp1MarkersDict.items():
        clusterCaseA = []
        clusterCaseB = []
        clusterCaseC = []
        clusterCaseD = []
        for k2, v2 in sp2MarkersDict.items():
            caseAcount = 0
            caseBcount = 0
            caseCcount = 0
            caseDcount = 0
            if RBHtupleList == "Decorator":
                for i in allConcatGenes:
                    if i in v1 and i in v2:
                        caseAcount += 1
                    elif i in v1 and i not in v2:
                        caseBcount += 1
                    elif i not in v1 and i in v2:
                        caseCcount += 1
                    elif i not in v1 and i not in v2:
                        caseDcount += 1 
            else:
                if orthoMode == "one2one" and contMethod == "orthopairsBased":
                    for i in RBHtupleList:
                        if (i[1] in v1 and i[0] in v2):
                            caseAcount += 1 
                        elif (i[1] not in v1 and i[0] in v2):
                            caseBcount += 1
                        elif (i[1] in v1 and i[0] not in v2):
                            caseCcount += 1
                        elif (i[1] not in v1 and i[0] not in v2):
                            caseDcount += 1 
                elif orthoMode == "many2many" and contMethod == "orthopairsBased":
                    for i in RBHtupleList:
                        intersectSp1 = [item for item in i[0] if item in v1]
                        intersectSp2 = [item for item in i[1] if item in v2]
                        ratioSp1 = len(intersectSp1)/len(i[0])
                        ratioSp2 = len(intersectSp2)/len(i[1]) 
                        if ratioSp1 != 0 and ratioSp2 != 0: 
                            caseAcount += ratioSp1 + ratioSp2 
                        elif ratioSp1 != 0 and ratioSp2 == 0:
                            caseBcount += ratioSp1
                        elif ratioSp1 == 0 and ratioSp2 != 0:
                            caseCcount += ratioSp2
                        elif ratioSp1 == 0 and ratioSp2 == 0:
                            caseDcount += 1 
                elif orthoMode == "SCO" and contMethod == "orthopairsBased":
                    pass
                elif orthoMode == "one2many" and contMethod == "orthopairsBased":
                    for i in RBHtupleList:
                        intersect=[item for item in i[1] if item in v2] 
                        if len(intersect) !=0 and i[0] in v1:
                            caseAcount += len(intersect) + 1
                        elif len(intersect) == 0 and i[0] in v1:
                            caseBcount += 1
                        elif len(intersect) != 0 and i[0] not in v1:
                            caseCcount += len(intersect)
                        elif len(intersect) == 0 and i[0] not in v1:
                            caseDcount += 1
                elif orthoMode == "one2one" and contMethod == "genomeBased":
                    for i in geneList:
                        if i in orthologsList:
                            correspondingGene = searchCorrespondingOrth(RBHtupleList,i,orthoMode)
                            if (i in v1 and correspondingGene in v2):
                                caseAcount +=1
                            elif (i in v1 and correspondingGene not in v2):
                                caseBcount +=1
                            elif (i not in v1 and correspondingGene in v2):
                                caseCcount += 1
                            elif (i not in v1 and correspondingGene not in v2):
                                caseDcount += 1
                        else:
                            if i in v1:
                                caseBcount += 1
                            else:
                                caseDcount += 1
                elif orthoMode == "many2many" and contMethod == "genomeBased": 
                    alreadyCountedGenes = []
                    for i in geneList:
                        if i in orthologsList:
                            correspondingManyOrthologsSp1,correspondingManyOrthologsSp2 = searchCorrespondingOrth(RBHtupleList,i,orthoMode)
                            if (correspondingManyOrthologsSp1 not in alreadyCountedGenes) and (correspondingManyOrthologsSp2 not in alreadyCountedGenes): 
                                ratioSp1 = sum(x in correspondingManyOrthologsSp1 for x in v2) / len(correspondingManyOrthologsSp1)                   
                                ratioSp2 = sum(x in correspondingManyOrthologsSp2 for x in v1) / len(correspondingManyOrthologsSp2) 
                                caseAcount += (ratioSp1 + ratioSp2) / 2
                                caseBcount += ratioSp1 
                                caseCcount += ratioSp2 
                                caseDcount += 2 - (ratioSp1 + ratioSp2) 
                                alreadyCountedGenes.append(correspondingManyOrthologsSp1) ; alreadyCountedGenes.append(correspondingManyOrthologsSp2)
            clusterCaseA.append(caseAcount)
            clusterCaseB.append(caseBcount)
            clusterCaseC.append(caseCcount)
            clusterCaseD.append(caseDcount)
        caseA.append(clusterCaseA)
        caseB.append(clusterCaseB)
        caseC.append(clusterCaseC)
        caseD.append(clusterCaseD)  
        print("A :",caseA) ; print("B :",caseB) ; print("C :",caseC) ; print("D :",caseD) 
    matrixList = []
    matrixDfCaseA = pd.DataFrame(caseA)
    matrixDfCaseA = matrixDfCaseA.values.tolist()
    matrixDfCaseB = pd.DataFrame(caseB)
    matrixDfCaseB = matrixDfCaseB.values.tolist()
    matrixDfCaseC = pd.DataFrame(caseC)
    matrixDfCaseC = matrixDfCaseC.values.tolist()
    matrixDfCaseD = pd.DataFrame(caseD)
    matrixDfCaseD = matrixDfCaseD.values.tolist()
    matrixList.append(matrixDfCaseA) ; matrixList.append(matrixDfCaseB)
    matrixList.append(matrixDfCaseC) ; matrixList.append(matrixDfCaseD) 
    return matrixList



"""
def getOrthofinderData(spList):
    many2manyList = []
    pathFile = "input/orthofinder/." 
    fileNames = os.listdir(pathFile)
    for i in fileNames:
        if spList[0]+"Proteins" in i:
            pathFile = pathFile.replace(".",i)
            pathFile += "/."
            fileNames = os.listdir(pathFile) 
    for i in fileNames:
        if spList[1]+"Proteins" in i:
            pathFile = pathFile.replace(".",i) 
            orthofinderOutput = pd.read_csv(pathFile,sep="\t")
            for index, row in orthofinderOutput.iterrows():
                sp1 = row[1] ; sp2 = row[2]
                if ', ' in sp1:
                    sp1 = sp1.split(', ')
                else:
                    sp1 = sp1.split()
                if ', ' in sp2:
                    sp2 = sp2.split(', ')
                else:
                    sp2 = sp2.split()
                many2manyList.append((sp1,sp2))
                with open('results/basic/orthofinder/many2many/'+spList[0]+'_'+spList[1]+'_orthologs.csv', "wb") as f:   
                    pickle.dump(many2manyList,f)
    return 
"""


"""
def basicOrthofinder(spList):
    used_combinations = []
    for i in spList:
        for j in spList:
            if i != j:
                if [i,j] not in used_combinations and [j,i] not in used_combinations:
                    spCombination = [i,j]
                    many2manyList = getOrthofinderData(spCombination)
                    used_combinations.append(spCombination)
"""
                


"""
species = None
spList = ["Clytia","Nematostella"]
lifestageList = ["medusa","polyp"]
contMethod="orthopairsBased"
orthoMode = "many2many"
geneSource = "orthofinder"
test="Fisher"
pValMethod="tippett"

#getOrthofinderData(spList)
RBHtupleList = manipData(spList,orthoMode,geneSource)
markersSubject1 = getAllMarkers(spList[0],lifestageList[0])
markersSubject2 = getAllMarkers(spList[1],lifestageList[1])
markersSubjectList = []
markersSubjectList.append(markersSubject1)
markersSubjectList.append(markersSubject2)
matrixList = genesCounting(RBHtupleList,markersSubjectList,orthoMode,contMethod,species,spList)
#many2manyList = getOrthofinderData(spList)
#print(many2manyList)
"""


listMatrices = []
testMatrixList = []



#print(matrixList)
#print(matrixList2)
"""
listMatrices.append(matrixList) ; listMatrices.append(matrixList2)
for m in range(len(listMatrices)):
    testMatrix = []
    for i in range(len(listMatrices[m][0])):
        listZero = []
        for j in range(len(listMatrices[m][0][i])):
            listZero.append(0)
        testMatrix.append(listZero)
    testMatrixList.append(testMatrix)
outputMatricesList = []

for m in range(len(listMatrices)):
    for i in range(len(testMatrixList[m])):
        for j in range(len(testMatrixList[m][i])):
            testDf = [[listMatrices[m][0][i][j],listMatrices[m][1][i][j]],[listMatrices[m][2][i][j],listMatrices[m][3][i][j]]]
            if (testDf[0][0] > testDf[0][1] and testDf[0][0] > testDf[1][0]):
                testDf[0][0] = max(testDf[0][1],testDf[1][0])
            pStat = ""
            pStat = stats.fisher_exact(testDf,alternative='greater')
            testMatrixList[m][i][j] = pStat[1]
    outputTestMatrix = pd.DataFrame(testMatrixList[m])
    outputMatricesList.append(outputTestMatrix)

resultMatrix = combinePvalues(outputMatricesList[0],outputMatricesList[1],"tippett")
resultMatrix = resultMatrix.values.tolist()
flattenMatrix = reduce(lambda a, b: a + b, resultMatrix) 
pAdjustedMatrix = sm.stats.multipletests(pvals=flattenMatrix,method="bonferroni")
colLen = len(resultMatrix[0]) 
rowLen = len(resultMatrix)
pAdjustedMatrix = pd.DataFrame(np.array(pAdjustedMatrix[1]).reshape(rowLen,colLen))
matrix = np.log10(pAdjustedMatrix) * -1 
matrix = pd.DataFrame.transpose(pAdjustedMatrix)
"""

def comparePPO(input1,input2):
    fasta1PPO = []
    fasta2PPO = []
    with open("results/final/"+input1+"/transcription_factors.fasta","r") as f:
        l = f.readline()
        while l != "":
            if l[0] == ">":
                gene = l.replace("\n","")
                gene = gene.replace(">","") 
                fasta1PPO.append(gene)
            l = f.readline()
    f.close()
    with open("results/final/"+input2+"/transcription_factors.fasta","r") as f:
        l = f.readline()
        while l != "":
            if l[0] == ">":
                gene = l.replace("\n","")
                gene = gene.replace(">","") 
                fasta2PPO.append(gene)
            l = f.readline()
    f.close() 
    intersectPPO = list(set(fasta1PPO) & set(fasta2PPO))
    unsharedPPO1 = list(set(fasta1PPO) - set(fasta2PPO))
    unsharedPPO2 = list(set(fasta2PPO) - set(fasta1PPO))  
    intersectFasta = getProteinSequences(intersectPPO,input1)
    unsharedPPO1Fasta = getProteinSequences(unsharedPPO1,input1)
    unsharedPPO2Fasta = getProteinSequences(unsharedPPO2,input2)
    saveFasta(input1,input2,'intersect',intersectFasta)
    saveFasta(input1,input2,'unshared1',unsharedPPO1Fasta)
    saveFasta(input1,input2,'unshared2',unsharedPPO2Fasta) 
    return 


def getProteinSequences(listPPO,fasta):
    proteomPath = "results/final/"+fasta+"/transcription_factors.fasta"
    with open(proteomPath,'r') as f:
        line = f.readline()
        prot=""
        protList = []
        geneList = []
        while line != "":
            if line[0] == ">":
                gene = line.split('>')[1]
                gene = gene.split('\n')[0]
                if gene in listPPO:
                    line = f.readline()
                    geneList.append(gene)
                    while line != "" and line[0] != ">":
                        prot += line
                        line = f.readline()
                    protList.append(prot)
                    prot = ""
                else:
                    line = f.readline()
            else:
                line = f.readline()
    f.close()
    geneProtDict = {}
    for i in range(len(geneList)):
        geneList[i] = ">" + geneList[i] +"\n"
        geneProtDict[geneList[i]] = protList[i] 
    return geneProtDict


def saveFasta(input1,input2,name,geneProtDict):
    with open("results/final/comparison/"+input1+"_"+input2+"_"+name+".fasta", 'a') as file:
        for k,v in geneProtDict.items():
            file.write(k)
            file.write(v)
    file.close()


comparePPO("RBH_one2one_orthopairsBased_cnidocytes","RBH_one2one_orthopairsBased_neurons")
