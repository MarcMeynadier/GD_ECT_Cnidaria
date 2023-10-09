import pandas as pd
import numpy as np
import random

def getRBHdata(spList):
    """
    Description
    -----------
    Retrieves RBH data between species and merge them in a dataframe.

    Parameters
    ----------
    spList
        list, contains species of subjects for comparison.
    
    Returns
    -------
    RBHdf,
        pandas.core.frame.DataFrame, dataframe of one-to-one orthologs between two species.
    """ 
    sp1 = spList[0] ; sp2 = spList[1]
    forward = pd.read_csv('input/RBH/'+sp1+'_'+sp2+'.txt',sep="\t")
    backward = pd.read_csv('input/RBH/'+sp2+'_'+sp1+'.txt',sep="\t") 
    backward[['qseqid','sseqid']] = backward[['sseqid','qseqid']] 
    RBHdf = pd.merge(forward,backward,on=['qseqid','sseqid'],how='inner')
    return RBHdf

def getAllMarkers(sp,lifestage):
    """
    Description
    -----------
    Retrieves genes markers of all cell clusters from Scanpy output of a subject.

    Parameters
    ----------
    sp
        str, species of subject.
    lifestage
        str, lifestage of subject
    
    Returns
    -------
    markersDictSp
        dict, keys correspond to cell clusters number and values to list of marker genes.
    """ 
    markersDictSp = {}
    pathMarkersSp = 'input/Scanpy/'+sp+'_'+lifestage+'_markers_1000g.csv'
    dfSp = pd.read_csv(pathMarkersSp,sep='\t') 
    markersSp = dfSp.groupby('clus')['markerGene'].apply(list).values
    for j in range(len(markersSp)):
        markersDictSp[j] = markersSp[j]
    return markersDictSp


def createRBHtuple(RBHdf):
    """
    Description
    -----------
    Creates a list of tuple where each tuple corresponds to pairs of one-to-one orthologs.

    Parameters
    ----------
    RBHdf,
        pandas.core.frame.DataFrame, dataframe of one-to-one orthologs between two species.
    
    Returns
    -------
    RBHgenesTupleList
        list, list of tuple where each tuple corresponds to pairs of one-to-one orthologs. 
    """ 
    RBHgenes1 = RBHdf.iloc[:,0].tolist() 
    RBHgenes2 = RBHdf.iloc[:,1].tolist() 
    RBHgenesTupleList = list(tuple(zip(RBHgenes1,RBHgenes2)))
    return RBHgenesTupleList


def getRBHtupleList(spList):
    """
    Description
    -----------
    Retrieves the list of tuple by calling createRBHtuple() from a list of two species.

    Parameters
    ----------
    spList
        list, contains species of subjects for comparison. 
    
    Returns
    -------
    RBHdfTupleList,
        list, list of tuple where each tuple corresponds to pairs of one-to-one orthologs.
    """ 
    if spList[0] != spList[1]:
        RBHdf = getRBHdata(spList) 
        RBHdfTupleList = createRBHtuple(RBHdf)
    else:
        RBHdfTupleList = "Decorator"
    return RBHdfTupleList


def shuffleTuple(RBHgenesTupleList):
    """
    Description
    -----------
    Shuffles the pairs of one-to-one orthologs to calculate a random contingency matrix. 

    Parameters
    ----------
    RBHgenesTupleList
        list, contains tuples of one-to-one orthologs between two species.
    
    Returns
    -------
    shuffleList
        list, list of tuples where the one-to-one orthologs has been shuffled.
    """
    shuffleList = []
    for i in range(len(RBHgenesTupleList)):
        randomTuple1 = random.choice(RBHgenesTupleList)
        randomTuple2 = random.choice(RBHgenesTupleList) 
        randomGene1 = randomTuple1[0]
        randomGene2 = randomTuple2[1]
        shuffleList.append((randomGene1,randomGene2))
    return shuffleList


def countOrthologsPairs(RBHgenesTupleList,markersDictSpList,reciprocality):
    """
    Description
    -----------
    Count the different number of cases for building contingency matrix, where cases depend on the reciprocality parameter. 

    Parameters
    ----------
    RBHgenesTupleList
        list, contains tuples of one-to-one orthologs between two species.
    markersDictSpList
        list, contains dictionnary for each species where keys of the dictionnaries are cluster number and values are list of marker genes.
    reciprocality,
        str, determines which case of the contingency matrix is calculated. 
    Returns
    -------
    matrixDf
        pandas.core.frame.DataFrame, dataframe that contains number of cases between every combination of cell clusters of the two subjects, where case depend on the reciprocality parameter. 
    """ 
    totalCount = []
    sp1MarkersDict = markersDictSpList[0]
    sp2MarkersDict = markersDictSpList[1]
    if RBHgenesTupleList != "Decorator":
        for k1,v1 in sp1MarkersDict.items(): 
            clusterMarkersCount = []
            for k2,v2 in sp2MarkersDict.items():
                count = 0
                countOrthologsPairs = 0
                for i in RBHgenesTupleList: 
                    if reciprocality=="yes":
                        if (i[1] in v1 and i[0] in v2):
                            countOrthologsPairs += 1
                    elif reciprocality=="no":
                        if (i[1] not in v1 and i[0] not in v2):
                            countOrthologsPairs += 1
                    elif reciprocality=="sp1":
                        if (i[1] in v1 and i[0] not in v2):
                            countOrthologsPairs += 1
                    elif reciprocality=="sp2":
                        if (i[1] not in v1 and i[0] in v2):
                            countOrthologsPairs += 1
                    count += 1 
                clusterMarkersCount.append(countOrthologsPairs)
            totalCount.append(clusterMarkersCount)
    else:    
        allConcatGenes = []
        for k1,v1 in sp1MarkersDict.items():
            for k2,v2 in sp2MarkersDict.items():
                concatGenes = list(set(v1 + v2))
                allConcatGenes = allConcatGenes + concatGenes
        allConcatGenes = list(set(allConcatGenes))
        if reciprocality=="yes":
            for k1,v1 in sp1MarkersDict.items(): 
                clusterCountShared = []
                for k2,v2 in sp2MarkersDict.items():
                    countGenesShared = 0 
                    for i in allConcatGenes:
                        if (i in v1 and i in v2):
                            countGenesShared+=1
                    clusterCountShared.append(countGenesShared)
                totalCount.append(clusterCountShared)  
        elif reciprocality=="sp1":
            for k1,v1 in sp1MarkersDict.items():
                clusterCountSp1 = []
                for k2,v2 in sp2MarkersDict.items():
                    countGenesSp1 = 0 
                    for i in allConcatGenes:
                        if (i in v1 and i not in v2):
                            countGenesSp1+=1
                    clusterCountSp1.append(countGenesSp1)
                totalCount.append(clusterCountSp1)   
        elif reciprocality=="sp2": 
            for k1,v1 in sp1MarkersDict.items():
                clusterCountSp2 = []
                for k2,v2 in sp2MarkersDict.items():
                    countGenesSp2 = 0 
                    for i in allConcatGenes:
                        if (i not in v1 and i in v2):
                            countGenesSp2+=1
                    clusterCountSp2.append(countGenesSp2)
                totalCount.append(clusterCountSp2)    
        elif reciprocality=="no": 
            for k1,v1 in sp1MarkersDict.items():
                clusterCountNone = []
                for k2,v2 in sp2MarkersDict.items():
                    countGenesNone = 0 
                    for i in allConcatGenes:
                        if (i not in v1 and i not in v2):
                            countGenesNone+=1
                    clusterCountNone.append(countGenesNone)
                totalCount.append(clusterCountNone)    
    matrixDf = pd.DataFrame(totalCount) 
    return matrixDf


def longestIsoforms(proteom):
    """
    Description
    -----------
    Extracts gene lists that code for longest isoforms by using proteom of species.

    Parameters
    ----------
    proteom
        str, path to proteom file.
    
    Returns
    -------
    geneList,
        list of genes that code for longest isoform, only one gene kept for all isoforms.
    """
    with open(proteom, 'r') as f:
        line = f.readline()
        geneList = np.array([])
        while line != "":
            if line[0] == ">": 
                gene = line.split('\n')[0] 
                gene = gene.replace(">","")
                geneList = np.append(geneList,gene) 
                line = f.readline()
                while line != "" and line[0] != ">":
                    line = f.readline()
        return geneList
    
 
def searchCorrespondingOrth(RBHgenesTupleList,gene):
    """
    Description
    -----------
    Searches for corresponding ortholog from the name of a specific gene.

    Parameters
    ----------
    RBHgenesTupleList
        list, contains tuples of one-to-one orthologs between two species.
    gene
        str, specific gene name.
    
    Returns
    -------
    correspondingOrth
        str, name of the ortholog corresponding to the input gene.
    """
    for i in RBHgenesTupleList: 
        if gene in i:
            if gene == i[0]:
                correspondingOrth = i[1]
                return correspondingOrth
            elif gene == i[1]:
                correspondingOrth = i[0]
                return correspondingOrth
            
             
def countContingencyTable(RBHgenesTupleList,markersDictSpList,species,spList):
    """
    Description
    -----------
    Builds contingency matrix in order to apply statistical test to test equivalence between cell clusters of two subjects.

    Parameters
    ----------
    RBHgenesTupleList
        list, contains tuples of one-to-one orthologs between two species.
    markersDictSpList
        list, contains dictionnary for each species where keys of the dictionnaries are cluster number and values are list of marker genes.
    species
        str, species of subject.
    spList
        list, contains species of subjects for comparison. 
    
    Returns
    -------
    matrixList
        list, contains dataframe for each cases of contingency matrix in order to apply statistical test to test equivalence between cell clusters of two subjects.
    """
    proteom = 'input/proteins/longest_isoform_'+species+'Proteins.fasta'
    geneList = longestIsoforms(proteom)
    if species == spList[0]:
        sp1MarkersDict = markersDictSpList[0]
        sp2MarkersDict = markersDictSpList[1]
    elif species == spList[1]:
        sp1MarkersDict = markersDictSpList[1]
        sp2MarkersDict = markersDictSpList[0]
    orthologsList = []
    if RBHgenesTupleList != "Decorator":
        for i in RBHgenesTupleList:
            orthologsList = np.append(orthologsList,i[0]) ; orthologsList = np.append(orthologsList,i[1])
    caseA = []
    caseB = []
    caseC = []
    caseD = []
    if RBHgenesTupleList != "Decorator":
        for k1,v1 in sp1MarkersDict.items():
            clusterCaseA = []
            clusterCaseB = []
            clusterCaseC = []
            clusterCaseD = []
            for k2,v2 in sp2MarkersDict.items():
                caseAcount = 0
                caseBcount = 0
                caseCcount = 0
                caseDcount = 0
                for i in geneList:
                    if i in orthologsList:
                        correspondingGene = searchCorrespondingOrth(RBHgenesTupleList,i)
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
                clusterCaseA.append(caseAcount)
                clusterCaseB.append(caseBcount) 
                clusterCaseC.append(caseCcount)
                clusterCaseD.append(caseDcount) 
            caseA.append(clusterCaseA)
            caseB.append(clusterCaseB)
            caseC.append(clusterCaseC)
            caseD.append(clusterCaseD)
    else: # Same species different lifestages
        allConcatGenes = []
        for k1,v1 in sp1MarkersDict.items():
            for k2,v2 in sp2MarkersDict.items():
                concatGenes = list(set(v1 + v2))
                allConcatGenes = allConcatGenes + concatGenes
        allConcatGenes = list(set(allConcatGenes))

        for k1,v1 in sp1MarkersDict.items():
            clusterCaseA = []
            clusterCaseB = []
            clusterCaseC = []
            clusterCaseD = []
            for k2,v2 in sp2MarkersDict.items():
                caseAcount = 0
                caseBcount = 0
                caseCcount = 0
                caseDcount = 0
            
                for i in allConcatGenes:
                    if i in v1 and i in v2:
                        caseAcount += 1
                    elif i in v1 and i not in v2:
                        caseBcount += 1
                    elif i not in v1 and i in v2:
                        caseCcount += 1
                    elif i not in v1 and i not in v2:
                        caseDcount += 1
                clusterCaseA.append(caseAcount)
                clusterCaseB.append(caseBcount) 
                clusterCaseC.append(caseCcount)
                clusterCaseD.append(caseDcount)
            caseA.append(clusterCaseA)
            caseB.append(clusterCaseB)
            caseC.append(clusterCaseC)
            caseD.append(clusterCaseD)               
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


                            
                   
        