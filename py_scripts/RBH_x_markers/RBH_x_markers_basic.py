from functools import reduce
import pandas as pd
import numpy as np
import random
import pickle
import os.path


def unique(list1):
    res = reduce(lambda re, x: re+[x] if x not in re else re, list1, [])
    return res


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



def getRBHdata(spList,orthoMode):
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
    outputList = [] 
    sp1 = spList[0] ; sp2 = spList[1]
    forward = pd.read_csv('input/RBH/'+sp1+'_'+sp2+'.txt',sep="\t")
    forward = forward.drop(forward[forward.evalue > 10e-14].index)
    forward = forward.drop(forward[forward.evalue == 0].index)
    backward = pd.read_csv('input/RBH/'+sp2+'_'+sp1+'.txt',sep="\t") 
    backward[['qseqid','sseqid']] = backward[['sseqid','qseqid']]
    backward = backward.drop(backward[backward.evalue > 10e-14].index)
    backward = backward.drop(backward[backward.evalue == 0].index) 
    forwardListTuple = [] ; backwardListTuple = []
    forwardListSimple = [] ; backwardListSimple = [] 
    if orthoMode == "one2one":
        RBHdf = pd.merge(forward,backward,on=['qseqid','sseqid'],how='inner')
        RBHgenes1 = RBHdf.iloc[:,0].tolist() 
        RBHgenes2 = RBHdf.iloc[:,1].tolist() 
        outputList = list(tuple(zip(RBHgenes1,RBHgenes2)))
        with open('results/basic/RBH/one2one/'+spList[0]+'_'+spList[1]+'_orthologs.csv', "wb") as f:   
            pickle.dump(outputList, f)
    elif orthoMode == "one2many":
        for index, row in forward.iterrows():
            tupleOrth = (row['qseqid'],row['sseqid'])
            forwardListTuple.append(tupleOrth) 
            forwardListSimple.append(row['sseqid']) 
        for i in forwardListSimple:
            orthList = []
            for j in forwardListTuple:
                if i in j:
                    if i == j[0]:
                        orthList.append(j[1])
                    elif i == j[1]:
                        orthList.append(j[0])
            if len(orthList) != 0:
                orthTuple = (i,orthList)
                outputList.append(orthTuple)
        with open('results/basic/RBH/one2many/'+spList[0]+'_'+spList[1]+'_orthologs.csv', "wb") as f:   
            pickle.dump(outputList,f)
    elif orthoMode == "many2one":
        for index, row in backward.iterrows():
            tupleOrth = (row['qseqid'],row['sseqid'])
            backwardListTuple.append(tupleOrth) 
            backwardListSimple.append(row['qseqid']) 
        for i in backwardListSimple:
            orthList = []
            for j in backwardListTuple:
                if i in j:
                    if i == j[0]:
                        orthList.append(j[1])
                    elif i == j[1]:
                        orthList.append(j[0])
            if len(orthList) != 0:
                orthTuple = (i,orthList)
                outputList.append(orthTuple)
        with open('results/basic/RBH/many2one/'+spList[0]+'_'+spList[1]+'_orthologs.csv', "wb") as f:   
            pickle.dump(outputList,f)  
    elif orthoMode == "many2many":
        for index, row in forward.iterrows():
            tupleOrth = (row['qseqid'],row['sseqid'])
            forwardListTuple.append(tupleOrth) 
            forwardListSimple.append(row['qseqid'])
        for index, row in backward.iterrows():
            tupleOrth = (row['qseqid'],row['sseqid'])
            backwardListTuple.append(tupleOrth) 
            backwardListSimple.append(row['sseqid']) 
            megaListForward = [] ; megaListBackward = []
        for i in forwardListSimple:
            testList = []
            for j in backwardListTuple:
                if i in j:
                    if j[0] == i:
                        testList.append(j[1])
                    elif j[1] == i:
                        testList.append(j[0])
                    tupleOrth = ([i],testList)
                    megaListForward.append(tupleOrth)
        for i in backwardListSimple:
            testList = []
            for j in forwardListTuple:
                if i in j:
                    if j[0] == i:
                        testList.append(j[1])
                    elif j[1] == i:
                        testList.append(j[0])
                    tupleOrth = (testList,[i])
                    megaListBackward.append(tupleOrth)
        many2manyTuple = ()
        for i in megaListForward:
            for j in megaListBackward:
                addFlag = False
                if len(np.intersect1d(i[0],j[0])) != 0:
                    addFlag = True
                    orthSp1 = i[0] + j[0]
                    orthSp2 = i[1] + j[1]
                if len(np.intersect1d(i[1],j[1])) != 0:
                    addFlag = True
                    orthSp1 = i[0] + j[0]
                    orthSp2 = i[1] + j[1] 
                if addFlag == True:
                    orthSp1 = unique(orthSp1)
                    orthSp2 = unique(orthSp2)
                    many2manyTuple = (orthSp1,orthSp2)
                    outputList.append(many2manyTuple)          
        outputList = unique(outputList)
        with open('results/basic/RBH/many2many/'+spList[0]+'_'+spList[1]+'_orthologs.csv', "wb") as f:   
            pickle.dump(outputList,f)
    return outputList



def getOrthofinderData(spList,orthoMode):
    outputList = []
    pathFile = "input/orthofinder/." 
    fileNames = os.listdir(pathFile)
    for i in fileNames:
        if spList[0]+"Proteins" in i:
            pathFile = pathFile.replace(".",i)
            pathFile += "/."
            fileNames = os.listdir(pathFile) 
    
    if orthoMode == "many2many":
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
                    outputList.append((sp1,sp2))
                    with open('results/basic/orthofinder/many2many/'+spList[0]+'_'+spList[1]+'_orthologs.csv', "wb") as f:   
                        pickle.dump(outputList,f)
    elif orthoMode == "one2one":
        print("\nOrthomode one2one is not supported with Orthofinder data.\n")
        exit()
    return outputList




def getTupleList(spList,orthoMode,geneSource):
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
        pathFile1 = 'results/basic/'+geneSource+'/'+orthoMode+'/'+spList[0]+'_'+spList[1]+'_orthologs.csv' 
        pathFile2 = 'results/basic/'+geneSource+'/'+orthoMode+'/'+spList[1]+'_'+spList[0]+'_orthologs.csv'  
        if os.path.exists(pathFile1) == False:
            if os.path.exists(pathFile2) == False:
                print("\nCurrently processing RBH output in "+orthoMode+" orthologs mode.")
                if geneSource == "RBH":
                    RBHdfTupleList = getRBHdata(spList,orthoMode)
                    return RBHdfTupleList
                elif geneSource == "orthofinder":
                    RBHdfTupleList = getOrthofinderData(spList,orthoMode)
                    return RBHdfTupleList
            else:
                with open(pathFile2, "rb") as f:
                    RBHdfTupleList = pickle.load(f)
                    return RBHdfTupleList 
        else:        
            with open(pathFile1, "rb") as f:
                RBHdfTupleList = pickle.load(f)
                return RBHdfTupleList 
    else:
        RBHdfTupleList = "sameSpecies"
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
    
 


def searchCorrespondingOrth(RBHgenesTupleList,gene,orthoMode):
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
    if orthoMode == "one2one":
        for i in RBHgenesTupleList: 
            if gene in i:
                if gene == i[0]:
                    correspondingOrth = i[1]
                    return correspondingOrth
                elif gene == i[1]:
                    correspondingOrth = i[0]
                    return correspondingOrth
    elif orthoMode == "many2many":
       for i in RBHgenesTupleList: 
            if gene in i[0]:
                correspondingOrth = i[1]
                return correspondingOrth,i[0]
            elif gene in i[1]:
                correspondingOrth = i[0]
                return correspondingOrth,i[1]
            

def genesCounting(RBHtupleList,markersDictSpList,orthoMode,contMethod,species,spList):
    spOrderFlag = blankOrderTest(RBHtupleList,markersDictSpList,orthoMode,contMethod,species,spList)
    sp1MarkersDict = markersDictSpList[0]
    sp2MarkersDict = markersDictSpList[1]
    caseA = []
    caseB = []
    caseC = []
    caseD = []
    if RBHtupleList == "sameSpecies":
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
            if RBHtupleList == "sameSpecies":
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
                elif orthoMode == "many2many":
                    for i in RBHtupleList:
                        if spOrderFlag == True:
                            intersectSp1 = [item for item in i[0] if item in v1]
                            intersectSp2 = [item for item in i[1] if item in v2]
                        else:
                            intersectSp1 = [item for item in i[0] if item in v2]
                            intersectSp2 = [item for item in i[1] if item in v1]
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


def blankOrderTest(RBHtupleList,markersDictSpList,orthoMode,contMethod,species,spList):
    sp1MarkersDict = markersDictSpList[0]
    sp2MarkersDict = markersDictSpList[1]
    caseA = []
    caseB = []
    caseC = []
    caseD = []
    if RBHtupleList == "sameSpecies":
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
            if RBHtupleList == "sameSpecies":
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
                elif orthoMode == "many2many":
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
            clusterCaseA.append(caseAcount)
            clusterCaseB.append(caseBcount)
            clusterCaseC.append(caseCcount)
            clusterCaseD.append(caseDcount)
        caseA.append(clusterCaseA)
        caseB.append(clusterCaseB)
        caseC.append(clusterCaseC)
        caseD.append(clusterCaseD) 
        spOrderFlag = None
        if (caseA == caseB) or (caseA==caseC) or (caseB==caseC):
            spOrderFlag = False
        else:
            spOrderFlag = True
        return spOrderFlag