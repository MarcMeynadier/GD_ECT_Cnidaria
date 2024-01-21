import pandas as pd
import warnings
import numpy as np
import matplotlib.pyplot as plt
import os
import re
import sys
import subprocess
from d3blocks import D3Blocks
from acrylic import Color, RANDOM

# Module importation

import Equivalencer_parameters as param
import Equivalencer_basic as basic


def getOrthologsPairs(RBHgenesTupleList,markersDictSpList,speciesList,lifestageList,geneSource,orthoMode,contMethod):
    """
    Description
    -----------
    Retrieves pairs of orthologs between cell clusters of subject and stores them in a daframe that is saved in CSV.

    Parameters
    ----------
    RBHgenesTupleList
        list, contains tuples of one-to-one orthologs between two species.
    markersDictSpList
        list, contains dictionnary for each species where keys of the dictionnaries are cluster number and values are list of marker genes.
    speciesList
        list, contains species of subjects for comparison.
    lifestageList
        list, contains lifestages of subjects for comparison. 

    Returns
    ----------
    None
    """
    emptyFlag = blankOrthologsPairs(RBHgenesTupleList,markersDictSpList,orthoMode)
    sp1MarkersDict = markersDictSpList[0]
    sp2MarkersDict = markersDictSpList[1]
    clusterOrthologsDf = []
    if RBHgenesTupleList != "sameSpecies":
        if orthoMode == "one2one": 
            for k1,v1 in sp1MarkersDict.items(): 
                clusterOrthologsTuple = []
                for k2,v2 in sp2MarkersDict.items():
                    orthologsPairsList = []
                    for i in RBHgenesTupleList: 
                        if (i[1] in v1 and i[0] in v2):
                            orthologsTuple = (i[1],i[0])
                            orthologsPairsList.append(orthologsTuple)
                    clusterOrthologsTuple.append(orthologsPairsList)
                clusterOrthologsDf.append(clusterOrthologsTuple) 
        elif orthoMode == "many2many":
            for k1,v1 in sp1MarkersDict.items(): 
                clusterOrthologsTuple = []
                for k2,v2 in sp2MarkersDict.items():
                    orthologsPairsList = []
                    for i in RBHgenesTupleList: 
                        if emptyFlag == False:
                            intersectSp1 = [item for item in i[1] if item in v1]
                            intersectSp2 = [item for item in i[0] if item in v2]
                        else:
                            intersectSp1 = [item for item in i[1] if item in v2]
                            intersectSp2 = [item for item in i[0] if item in v1]
                        if len(intersectSp1) != 0 and len(intersectSp2) != 0:
                            orthologsTuple = (intersectSp1,intersectSp2)
                            orthologsPairsList.append(orthologsTuple)
                    clusterOrthologsTuple.append(orthologsPairsList)
                clusterOrthologsDf.append(clusterOrthologsTuple)
    else:
        clusterOrthologsDf = []
        for k1,v1 in sp1MarkersDict.items(): 
            clusterOrthologsShared = []
            for k2,v2 in sp2MarkersDict.items():
                orthologsSharedList = list(set(v1) & set(v2))
                clusterOrthologsShared.append(orthologsSharedList)
            clusterOrthologsDf.append(clusterOrthologsShared)
    matrixOrthologs = pd.DataFrame(clusterOrthologsDf)
    if RBHgenesTupleList != "sameSpecies":            
        matrixOrthologs.to_csv('results/id/'+geneSource+'/'+orthoMode+'/'+contMethod+'/'+speciesList[0]+'_'+lifestageList[0]+'_'+speciesList[1]+'_'+lifestageList[1]+"_orthologsId.csv")
    else:
        matrixOrthologs.to_csv('results/id/sameSpecies/'+speciesList[0]+'_'+lifestageList[0]+'_'+speciesList[1]+'_'+lifestageList[1]+"_orthologsId.csv")
    print("Orthologs shared lists processed for "+speciesList[0]+' '+lifestageList[0]+' & '+speciesList[1]+' '+lifestageList[1]) 
    return


def blankOrthologsPairs(RBHgenesTupleList,markersDictSpList,orthoMode):
    sp1MarkersDict = markersDictSpList[0]
    sp2MarkersDict = markersDictSpList[1]
    clusterOrthologsDf = []
    if RBHgenesTupleList != "sameSpecies":
        if orthoMode == "one2one": 
            for k1,v1 in sp1MarkersDict.items(): 
                clusterOrthologsTuple = []
                for k2,v2 in sp2MarkersDict.items():
                    orthologsPairsList = []
                    for i in RBHgenesTupleList: 
                        if (i[1] in v1 and i[0] in v2):
                            orthologsTuple = (i[1],i[0])
                            orthologsPairsList.append(orthologsTuple)
                    clusterOrthologsTuple.append(orthologsPairsList)
                clusterOrthologsDf.append(clusterOrthologsTuple) 
        elif orthoMode == "many2many":
            for k1,v1 in sp1MarkersDict.items(): 
                clusterOrthologsTuple = []
                for k2,v2 in sp2MarkersDict.items():
                    orthologsPairsList = []
                    for i in RBHgenesTupleList: 
                        intersectSp1 = [item for item in i[1] if item in v1]
                        intersectSp2 = [item for item in i[0] if item in v2]
                        if len(intersectSp1) != 0 and len(intersectSp2) != 0:
                            orthologsTuple = (intersectSp1,intersectSp2)
                            orthologsPairsList.append(orthologsTuple)
                    clusterOrthologsTuple.append(orthologsPairsList)
                clusterOrthologsDf.append(clusterOrthologsTuple)
                emptyFlag = None
                for i in clusterOrthologsDf:
                    empty = all(not sous_liste for sous_liste in i)
                    if empty:
                        emptyFlag = True
                    else:
                        emptyFlag = False
                return emptyFlag
    else:
        clusterOrthologsDf = []
        for k1,v1 in sp1MarkersDict.items(): 
            clusterOrthologsShared = []
            for k2,v2 in sp2MarkersDict.items():
                orthologsSharedList = list(set(v1) & set(v2))
                clusterOrthologsShared.append(orthologsSharedList)
            clusterOrthologsDf.append(clusterOrthologsShared)
            

def parseAnnotFile(annotationFile):
    """
    Description
    -----------
    

    Parameters
    ----------
    

    Returns
    ----------
    
    """
    annotDict = {}
    with open('input/annotation/'+annotationFile,'r') as f:
        l = f.readline()
        subject = ""
        annot = ""
        while l != "":
            if l[0] != "[":
                subject = l.split("\n")[0]
                l = f.readline()
            else:
                annot += l
                l = f.readline()
            annot = annot.split("\n")[0]
            annot = annot.replace("\"","")
            annot = annot.replace("[","")
            annot = annot.replace("]","")
            annot = annot.split(",")
            annotDict[subject] = annot
            annot = ""
    return annotDict


def selectAnnotFile():
    filesList = []
    filesList2 = []
    for i in os.listdir("input/annotation/"):
        if i.endswith(".txt"):
            filesList2.append(i)
            filesList.append(i.split(".")[0])
    print("\nSelect which file you want as annotation file :\n")
    for i in range(len(filesList)):
        filesList[i] = str(i+1) + " - " + filesList[i]
    print(' | '.join(filesList))
    filesList.append("dummy")
    lock = False
    while lock == False:
        choice = int(input())
        if (choice-1) in list(range(len(filesList))) and choice != 0:
            try:
                annotationFile = filesList2[choice-1]
                lock = True
            except IndexError:
                print("You have to choose a number in the list, or select 0 to return to menu")
        elif choice == 0:
            return
        else:
            print("You have to choose a number in the list, or select 0 to return to menu")
    print("\nAnnotation file chosen : "+annotationFile)
    return annotationFile


def applyAnnotation(matrixOrthologs,speciesList,lifestageList):
    """
    Description
    -----------
    

    Parameters
    ----------
    

    Returns
    ----------
    
    """
    annotationFile = selectAnnotFile()
    annotDict = parseAnnotFile(annotationFile)
    subjectAnnot = []
    sp1 = speciesList[0] ; sp2 = speciesList[1]
    lf1 = lifestageList[0] ; lf2 = lifestageList[1]
    annotSp1 = sp1+" "+lf1
    annotSp2 = sp2+" "+lf2 
    try:
        matrixOrthologs.index = annotDict[annotSp1]
        subjectAnnot.append(annotDict[annotSp1])
    except ValueError:
        matrixOrthologs.index = annotDict[annotSp2]
        subjectAnnot.append(annotDict[annotSp2]) 
    try:
        matrixOrthologs.columns = annotDict[annotSp2]
        subjectAnnot.append(annotDict[annotSp2]) 
    except ValueError:
        matrixOrthologs.columns = annotDict[annotSp1]
        subjectAnnot.append(annotDict[annotSp1])
    return matrixOrthologs,subjectAnnot


def chooseClusters(subjectAnnot,speciesList,lifestageList):
    """
    Description
    -----------
    Retrieves pairs of orthologs between cell clusters of subject and stores them in a daframe that is saved in CSV.

    Parameters
    ----------
    RBHgenesTupleList
        list, contains tuples of one-to-one orthologs between two species.
    markersDictSpList
        list, contains dictionnary for each species where keys of the dictionnaries are cluster number and values are list of marker genes.
    speciesList
        list, contains species of subjects for comparison.
    lifestageList
        list, contains lifestages of subjects for comparison. 

    ----------
    None
    """ 
    print("\nClusters list of "+speciesList[0]+" "+lifestageList[0])
    print("Please select clusters you want to extract by their numbers (delimited by space if several)\n")
    print(' | '.join(subjectAnnot[0]))
    extractClusterList1 = input()
    m1 = re.search("\d+ \d+",extractClusterList1)
    if m1:
        extractClusterList1 = extractClusterList1.rstrip()
        extractClusterList1 = extractClusterList1.split(" ")
        extractClusterList1 = [int(x) for x in extractClusterList1]
    else:
        extractClusterList1 = [int(extractClusterList1)]
    print("\nClusters list of "+speciesList[1]+" "+lifestageList[1])
    print("Please select clusters you want to extract by their numbers (delimited by space if several)\n")
    print(' | '.join(subjectAnnot[1]))
    extractClusterList2 = input()
    m2 = re.search("\d+ \d+",extractClusterList2) 
    if m2:
        extractClusterList2 = extractClusterList2.rstrip()
        extractClusterList2 = extractClusterList2.split(" ")
        extractClusterList2 = [int(x) for x in extractClusterList2]
    else:
        extractClusterList2 = [int(extractClusterList2)]
    return extractClusterList1,extractClusterList2


def automaticIdentityStartingPoint(SPS):
    
    """
    Description
    -----------
    Retrieves from the input folder the txt file of Startin Point Subject (SPS).

    Parameters
    ----------
    None

    ----------
    clusterList
        list, list of cell clusters from the SPS for cell type of interest.
    geneList
        list, list of genes from the SPS for cell type of interest. 
    species
        str, species of SPS.
    lifestage
        str, lifestage of SPS.
    """
    with open('input/SPS/'+SPS) as f: 
        clusterList = f.readline()
        clusterList = clusterList.replace('\n','')
        clusterList=clusterList.split(',')
        geneList = f.readline()
        geneList = geneList.replace('\n','')
        geneList = geneList.split(',')
        species = f.readline()
        species = species.replace('\n','')
        lifestage = f.readline()
        lifestage = lifestage.replace('\n','')
    return clusterList,geneList,species,lifestage
    

def trimDuplicateOrthologs(sharedList):
    """
    Description
    -----------
    Takes a list of orthologs and trim every duplicates. 

    Parameters
    ----------
    sharedList
         list, list of orthologs.

    Returns
    ----------
    sharedList
        list, list of orthologs without duplicate.
    """
    sharedList = [x + ',' for x in sharedList]
    sharedList = [sub.replace(' ','') for sub in sharedList]
    sharedList = ''.join(sharedList)
    sharedList = sharedList.split(",")
    sharedList= sharedList[:-1]
    sharedList = list(set(sharedList))
    sharedList.sort()
    return sharedList


def getProteinSequences(species,orthologsList):
    """
    Description
    -----------
    Retrieves protein sequences of specific genes provided by a list, by using proteom of a species that is available in input folder.

    Parameters
    ----------
    species
        str, species name.
    orthologsList
        list, list of orthologs.

    Returns
    ----------
    geneProtDict
        dict, keys are genes names and values are proteins sequences.
    """
    proteomPath = 'input/proteins/longest_isoform_'+species+'Proteins.fasta'
    with open(proteomPath,'r') as f:
        line = f.readline()
        prot=""
        protList = []
        geneList = []
        while line != "":
            if line[0] == ">":
                gene = line.split('>')[1]
                gene = gene.split('\n')[0]
                if gene in orthologsList:
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


def manageProteinSequences(listOfUsedSpecies,orthologsList):
    """
    Description
    -----------
    Retrieves protein sequences of specific genes provided by a list, by using proteom of every species in the input list.

    Parameters
    ----------
    listOfUsedSpecies
        list, list of species already used for comparison. 
    orthologsList
        list, list of orthologs.

    Returns
    ----------
    fastaFinalDict
        dict, keys are genes names and values are proteins sequences.
    """
    geneProtDictList = []
    for i in listOfUsedSpecies:
        geneProtDict = getProteinSequences(i,orthologsList)
        geneProtDictList.append(geneProtDict)
    fastaFinalDict = {}
    for i in geneProtDictList:
        fastaFinalDict.update(i)
    return fastaFinalDict


def saveProteinsFasta(fastaFinalDict,outputName):
    """
    Description
    -----------
    Saves in a FASTA file orthologs and their corresponding protein sequences.

    Parameters
    ----------
    fastaFinalDict
        dict, keys are orthologs and values are protein sequences.
    processType
         str, defines the processing type (manual, semi-automatic or automatic). 

    Returns
    ----------
    None
    """
    while outputName == None:
        print("\nPlease choose a code name for your results\n")
        outputName = input()
    outputPath = 'results/final/'+outputName
    lock = False
    while lock == False:
        if not os.path.exists(outputPath):
            os.makedirs(outputPath)
            lock = True
        else: 
            print("\nThere is already an output folder with this name, please choose another one\n")
            outputName = input()
            outputPath = 'results/final/'+outputName
    with open(outputPath+'/shared_orthologs.fasta', 'w') as file:
        for k,v in fastaFinalDict.items():
            file.write(k)
            file.write(v)
    file.close()
    print("\nFasta file saved in folder "+outputName)
    return
    

def extractSharedOrthologs(matrixOrthologs,subjectAnnot,speciesList,lifestageList,processType,extractClusterList1,verifiedCluster,startingSpecies,startingLifestage):
    """
    Description
    -----------
    Extracts shared orthologs between lists of cell clusters of two subjects.

    Parameters
    ----------
    subjectAnnot
        list, manual annotation of cell clusters of a subject.
    speciesList
        list, contains species of subjects for comparison.
    lifestageList
        list, contains lifestages of subjects for comparison. 
    processType
        str, defines the processing type (manual, semi-automatic or automatic).
    extractClusterList1
        list, list that contains cell clusters number of species 1.
    verifierCluster
        list, list that contains cell clusters number of species 2.
    startingSpecies
        str, species name of subject 1.
    startingLifestage
        str, lifestage of subject 1.

    Returns
    ----------
    extractSharedOrthologsList
        list, contains shared orthologs between list of cell clusters of two subjects.
    """
    corresp = None
    if processType == "m":
        extractClusterList1,extractClusterList2 = chooseClusters(subjectAnnot,speciesList,lifestageList)
    elif processType=="a" or processType=="sa" or processType == "f":  
        extractClusterList2 = verifiedCluster
        if (speciesList[0] == startingSpecies and lifestageList[0] == startingLifestage):
            corresp = "\nCorresponding cells clusters of "+speciesList[1]+" "+lifestageList[1]+" when compared with "+speciesList[0]+" "+lifestageList[0]+" : "+', '.join([str(x) for x in extractClusterList2])
            print(corresp)
        else:
            print("\nRunning comparison between "+speciesList[1]+" "+lifestageList[1]+" and "+speciesList[0]+" "+lifestageList[0])
    try:
        extractSharedOrthologs = matrixOrthologs.iloc[extractClusterList1,extractClusterList2]
    except IndexError:
        extractSharedOrthologs = matrixOrthologs.iloc[extractClusterList2,extractClusterList1]
    if speciesList[0] != speciesList[1]:
        extractSharedOrthologsStr = ""
        for i in extractSharedOrthologs.iterrows():
            for j in i[1]:
                j = j.replace("]","")
                j = j.replace("[","")
                extractSharedOrthologsStr += j 
        extractSharedOrthologsList = re.findall(r'\(.*?\)',extractSharedOrthologsStr)
        extractSharedOrthologsList = list(set(extractSharedOrthologsList))
        extractSharedOrthologsList = [sub.replace("'",'') for sub in extractSharedOrthologsList] 
        extractSharedOrthologsList = [sub.replace('(','') for sub in extractSharedOrthologsList]
        extractSharedOrthologsList = [sub.replace(')','') for sub in extractSharedOrthologsList] 
    else:
        extractSharedOrthologsList = []
        for i in extractSharedOrthologs.iterrows():
            for j in i[1]:
                j = j.replace("]","")
                j = j.replace("[","")
                j = j.replace("'","")
                j = j.replace(" ","")
                j = j.split(",")
                extractSharedOrthologsList = extractSharedOrthologsList + j
        extractSharedOrthologsList = list(set(extractSharedOrthologsList))
    return extractSharedOrthologsList,corresp


def getNewOrthologsMatrix(speciesList,lifestageList,clusterListIntegers,processType,geneList,startingSpecies,startingLifestage,contMethod,test,pValMethod,subjectList):
    """
    Description
    -----------
    Gets new combination of subjects to extract the shared orthologs between their equivalent cell clusters.

    Parameters
    ----------
    speciesList
        list, contains species of subjects for comparison.
    lifestageList
        list, contains lifestages of subjects for comparison. 
    clusterListIntegers
        list, list that contains cell clusters number of subject 1. 
    processType
        str, defines the processing type (manual, semi-automatic or automatic).
    geneList
        list, list of marker genes belonging to the cell type of interest of the SPS.
    startingSpecies
        str, species name of subject 1.
    startingLifestage
        str, lifestage of subject 1.

    Returns
    ----------
    newSharedOrthologsSubject
        list, contains shared orthologs between list of cell clusters of two subjects.
    newSpeciesList
        list, contains species that have been used for comparison.
    """
    if processType == "m":
        clusterListIntegers = None
        verifiedCluster = None
        newSpeciesList = param.chooseParametersSpecies(subjectList)
        newLifestageList = param.chooseParametersLifestage(subjectList)
        try: 
            matrixOrthologs = pd.read_csv('results/id/'+newSpeciesList[0]+'_'+newLifestageList[0]+'_'+newSpeciesList[1]+'_'+newLifestageList[1]+"_orthologsId.csv",index_col=0) 
        except FileNotFoundError:
            matrixOrthologs = pd.read_csv('results/id/'+newSpeciesList[1]+'_'+newLifestageList[1]+'_'+newSpeciesList[0]+'_'+newLifestageList[0]+"_orthologsId.csv",index_col=0) 
        matrixOrthologs,newSubjectAnnot = applyAnnotation(matrixOrthologs,newSpeciesList,newLifestageList) 
    elif processType == "a" or processType=="sa" or processType == "f":  
        newSubjectAnnot = None 
        newSp,newLifestage = param.chooseNewSpeciesLifestage(subjectList)
        try:
            if contMethod == "orthopairsBased":
                statMatrix = pd.read_csv('results/stat/'+speciesList[0]+'_'+lifestageList[0]+'_'+newSp+'_'+newLifestage+"_"+test+"_test_1000g_"+contMethod+".csv",index_col=0)
            elif contMethod == "genomeBased":
                statMatrix = pd.read_csv('results/stat/'+speciesList[0]+'_'+lifestageList[0]+'_'+newSp+'_'+newLifestage+"_"+test+"_test_1000g_"+contMethod+"_"+pValMethod+".csv",index_col=0)
            matrixOrthologs = pd.read_csv('results/id/'+speciesList[0]+'_'+lifestageList[0]+'_'+newSp+'_'+newLifestage+"_orthologsId.csv",index_col=0)
        except FileNotFoundError:
            if contMethod == "orthopairsBased":
                statMatrix = pd.read_csv('results/stat/'+newSp+'_'+newLifestage+'_'+speciesList[0]+'_'+lifestageList[0]+"_"+test+"_test_1000g_"+contMethod+".csv",index_col=0) 
            elif contMethod == "genomeBased":
                statMatrix = pd.read_csv('results/stat/'+newSp+'_'+newLifestage+'_'+speciesList[0]+'_'+lifestageList[0]+"_"+test+"_test_1000g_"+contMethod+"_"+pValMethod+".csv",index_col=0) 
            matrixOrthologs = pd.read_csv('results/id/'+newSp+'_'+newLifestage+'_'+speciesList[0]+'_'+lifestageList[0]+"_orthologsId.csv",index_col=0)
        newSpeciesList = [speciesList[0],newSp]
        newLifestageList = [lifestageList[0],newLifestage]
        verifiedCluster = verifyClusterIdentity(statMatrix,clusterListIntegers,newSpeciesList,newLifestageList,geneList)
    newSharedOrthologsSubject,corresp = extractSharedOrthologs(matrixOrthologs,newSubjectAnnot,newSpeciesList,newLifestageList,processType,clusterListIntegers,verifiedCluster,startingSpecies,startingLifestage)
    return newSharedOrthologsSubject,newSpeciesList


def parseSharedOrthologsList(sharedOrthologsList):
    """
    Description
    -----------
    Parses list of tupples of orthologs into list of single elements.

    Parameters
    ----------
    sharedOrthologsList
        list, list of tuples.
    Returns
    ----------
    sharedOrthologsList
        list, list of str.
    """
    delList = []
    addList = []
    for i in range(len(sharedOrthologsList)):
        if ',' in sharedOrthologsList[i]:
            g1 = sharedOrthologsList[i].split(',')[0]
            g2 = sharedOrthologsList[i].split(',')[1]
            delList.append(i) ; addList.append(g1.strip()) ; addList.append(g2.strip())
    for i in sorted(delList, reverse=True):
        del sharedOrthologsList[i]
    sharedOrthologsList = sharedOrthologsList + addList
    return sharedOrthologsList


def continueOrthologsComparison(speciesList,lifestageList,clusterListIntegers,sharedOrthologsList,dfListSharedOrthologs,listOfUsedSpecies,processType,geneList,startingSpecies,startingLifestage,contMethod,test,pValMethod,subjectList):
    """
    Description
    -----------
    Takes previous shared orthologs results and does another comparison of shared orthologs between different subjects.

    Parameters
    ----------
    speciesList
        list, contains species of subjects for comparison.
    lifestageList
        list, contains lifestages of subjects for comparison. 
    clusterListIntegers
        list, list that contains cell clusters number of subject 1.    
    sharedOrthologsList
        list, list of shared orthologs between species.
    dfListSharedOrthologs
        list, contains lists of shared orthologs between two compared subjects.
    listOfUsedSpecies
        list, list of species already used for comparison.

    Returns
    ----------
    sharedOrthologsList
        list, list of shared orthologs between species.
    dfListSharedOrthologs
        list, contains lists of shared orthologs between two compared subjects.
    listOfUsedSpecies
        list, list of species already used for comparison.  
    """
    newOrthologsList,newSpeciesList = getNewOrthologsMatrix(speciesList,lifestageList,clusterListIntegers,processType,geneList,startingSpecies,startingLifestage,contMethod,test,pValMethod,subjectList)
    for i in newSpeciesList:
        if i not in listOfUsedSpecies:
            listOfUsedSpecies.append(i) 
    newOrthologsList = parseSharedOrthologsList(newOrthologsList) 
    dfListSharedOrthologs.append(newOrthologsList)
    print("\n")
    print(dfListSharedOrthologs)
    sharedOrthologsList = sharedOrthologsList + newOrthologsList 
    print("\n")
    print(len(sharedOrthologsList))
    res = ""
    while (res != "y" or res !="n"):
        print("\nDo you want to continue with another subject (y / n)?")
        res = input()
        if res == "n":
            return sharedOrthologsList,dfListSharedOrthologs,listOfUsedSpecies
        elif res == "y":
            sharedOrthologsList,dfListSharedOrthologs,listOfUsedSpecies = continueOrthologsComparison(speciesList,lifestageList,clusterListIntegers,sharedOrthologsList,dfListSharedOrthologs,listOfUsedSpecies,processType,geneList,startingSpecies,startingLifestage,contMethod,test,pValMethod,subjectList) 


def trimOrthologsMatricesManual(speciesList,lifestageList,processType,outputName,subjectList,geneSource,orthoMode,contMethod):
    """
    Description
    -----------
    From the results of shared orthologs between cell clusters of subjects, allows the user to manually choose which cell clusters he wants to compare to extract shared orthologs, at the condition that they are already annotated.

    Parameters
    ----------
    speciesList
        list, contains species of subjects for comparison.
    lifestageList
        list, contains lifestages of subjects for comparison. 
    processType
        str, defines the processing type (manual, semi-automatic or automatic). 

    Returns
    ----------
    None
    """
    geneList = None
    clusterListIntegers = None
    verifiedCluster = None
    startingSpecies = None
    startingLifestage = None
    #contMethod = None
    test = None
    pValMethod = None
    correspList = []
    listOfUsedSpecies = speciesList.copy()
    if speciesList[0] != speciesList[1]:
        try:
            matrixOrthologs = pd.read_csv('results/id/'+geneSource+'/'+orthoMode+'/'+contMethod+'/'+speciesList[0]+'_'+lifestageList[0]+'_'+speciesList[1]+'_'+lifestageList[1]+"_orthologsId.csv",index_col=0)
        except FileNotFoundError:
            matrixOrthologs = pd.read_csv('results/id/'+geneSource+'/'+orthoMode+'/'+contMethod+'/'+speciesList[1]+'_'+lifestageList[1]+'_'+speciesList[0]+'_'+lifestageList[0]+"_orthologsId.csv",index_col=0)
    else:
        try:
            matrixOrthologs = pd.read_csv('results/id/sameSpecies/'+contMethod+'/'+speciesList[0]+'_'+lifestageList[0]+'_'+speciesList[1]+'_'+lifestageList[1]+"_orthologsId.csv",index_col=0)
        except FileNotFoundError:
            matrixOrthologs = pd.read_csv('results/id/sameSpecies/'+speciesList[1]+'_'+lifestageList[1]+'_'+speciesList[0]+'_'+lifestageList[0]+"_orthologsId.csv",index_col=0)
    matrixOrthologs,subjectAnnot1 = applyAnnotation(matrixOrthologs,speciesList,lifestageList)
    sharedOrthologsList,corresp = extractSharedOrthologs(matrixOrthologs,subjectAnnot1,speciesList,lifestageList,processType,clusterListIntegers,verifiedCluster,startingSpecies,startingLifestage)
    correspList.append(corresp)
    dfListSharedOrthologs = []
    sharedOrthologsList = parseSharedOrthologsList(sharedOrthologsList) 
    dfListSharedOrthologs.append(sharedOrthologsList)
    res = ""
    while (res != "y" or res !="n"):
        print("\nDo you want to continue with another subject (y / n)?")
        res = input()
        if res == "n":
            sharedOrthologsList = trimDuplicateOrthologs(sharedOrthologsList)
            print(sharedOrthologsList)
            fastaFinalDict = manageProteinSequences(listOfUsedSpecies,sharedOrthologsList)
            saveProteinsFasta(fastaFinalDict,outputName)
            #subjectNameBarplot2 = subjectNameBarplot.copy()
            #horizontalBarplotGenes(nbrGenesList,subjectNameBarplot,outputName)
            #horizontalBarplotOrthologs(nbrOrthologsList,subjectNameBarplot2,sharedOrthologsListConcat,outputName)
            #logOutput(correspList,processType,contMethod,test,pValMethod,iterationNbr,SPS,outputName)
            #chordDiagram(subjectOrderChordDiagram,dfListSharedOrthologs)
            return
        elif res == "y":
            sharedOrthologsList,dflistOfUsedSpecies,listOfUsedSpecies = continueOrthologsComparison(speciesList,lifestageList,clusterListIntegers,sharedOrthologsList,dfListSharedOrthologs,listOfUsedSpecies,processType,geneList,startingSpecies,startingLifestage,contMethod,test,pValMethod,subjectList)
            sharedOrthologsList = trimDuplicateOrthologs(sharedOrthologsList)
            print(sharedOrthologsList)
            fastaFinalDict = manageProteinSequences(listOfUsedSpecies,sharedOrthologsList)
            saveProteinsFasta(fastaFinalDict,outputName)
            #subjectNameBarplot2 = subjectNameBarplot.copy()
            #horizontalBarplotGenes(nbrGenesList,subjectNameBarplot,outputName)
            #horizontalBarplotOrthologs(nbrOrthologsList,subjectNameBarplot2,sharedOrthologsListConcat,outputName)
            #logOutput(correspList,processType,contMethod,test,pValMethod,iterationNbr,SPS,outputName)
            #chordDiagram(subjectOrderChordDiagram,dfListSharedOrthologs)
            return 


def extractClusterLowerPvalue(statMatrix,inputClusterList):
    """
    Description
    -----------
    From a list of cell clusters of subject 1 and statistical results between subject 1 and subject 2, extracts cells clusters of subject 2 where the p-value is the lowest compared to input cell clusters of subject 1.

    Parameters
    ----------
    statMatrix
        pandas.core.frame.DataFrame, dataframe of p-value between every combination of cell clusters of two subjects.
    inputClusterList
        list, list of cell clusters number belonging to subject 1 that will be examined.

    Returns
    ----------
    outputClusterList
        list, list of cell clusters number of subject 2 that might correspond to same cell type as cells clusters number of subject 1.
    """
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        statMatrix = np.log10(statMatrix)
    outputClusterList = []
    try:
        for i in inputClusterList:
            rowStatTest = statMatrix.iloc[int(i)]
    except IndexError:
        statMatrix = statMatrix.T
    for i in inputClusterList:
        rowStat = statMatrix.iloc[int(i)]
        rowStat = rowStat.abs()
        rowThreshold = rowStat[rowStat != 0]
        rowThreshold = rowThreshold[rowThreshold != np.inf] 
        rowThreshold = (rowThreshold.max() - rowThreshold.mean())
        rowStatCured = rowStat[rowStat >= rowThreshold]
        for j in rowStatCured:
            outputClust = statMatrix.columns[(statMatrix == -j).iloc[int(i)]]
            outputClusterList.append(int(outputClust[0]))
    outputClusterList = list(set(outputClusterList))
    return outputClusterList


def verifyClusterIdentity(statMatrix,inputClusterList,speciesList,lifestageList,geneList,orthoMode,geneSource):
    """
    Description
    -----------
    Determines which cell clusters of subject 2 belong to the same cell type as cell clusters of subject 1.
    It uses both statistical results between subjects and marker genes of SPS as a double threshold.

    Parameters
    ----------
    statMatrix
        pandas.core.frame.DataFrame, dataframe of p-value between every combination of cell clusters of two subjects.
    inputClusterList
        list, list of cell clusters number belonging to subject 1 that will be examined.
    speciesList
        list, contains species of subjects for comparison.
    lifestageList
        list, contains lifestages of subjects for comparison.
    geneList
        list, list of marker genes belonging to the cell type of interest of the SPS. 

    Returns
    ----------
    verifiedCluster
        list, list of cell clusters number of subject 2 that correspond to same cell type as cells clusters number of subject 1.
    """
    verifiedCluster = []
    outputClusterList = extractClusterLowerPvalue(statMatrix,inputClusterList)
    scanpyOutput = basic.getAllMarkers(speciesList[1],lifestageList[1]) 
    if speciesList[0] != speciesList[1]:
        targetGenes = []
        RBHgenesTupleList = basic.getTupleList(speciesList,orthoMode,geneSource) 
        for i in RBHgenesTupleList:
            if orthoMode == "one2one":
                if i[0] in geneList:
                    targetGenes.append(i[1])
                elif i[1] in geneList:
                    targetGenes.append(i[0])
            elif orthoMode == "many2many":
                intersectSp1 = [item for item in i[0] if item in geneList]
                intersectSp2 = [item for item in i[1] if item in geneList]
                if len(intersectSp1) != 0:
                    targetGenes = targetGenes + i[1]
                elif len(intersectSp2) != 0:
                    targetGenes = targetGenes + i[0] 
                #print(targetGenes)
        for i in outputClusterList:
            for k,v in scanpyOutput.items():
                sharedGenesList = []
                if k == i:
                    for j in v:
                        if j in targetGenes:
                            sharedGenesList.append(j)
                if len(sharedGenesList) != 0:
                    verifiedCluster.append(k)
    else:
        for i in outputClusterList:
            for k,v in scanpyOutput.items():
                sharedGenesList = []
                if k == i:
                    for j in v:
                        if j in geneList:
                            sharedGenesList.append(j)
                if len(sharedGenesList) != 0:
                    verifiedCluster.append(k)
    verifiedCluster = list(set(verifiedCluster))
    return verifiedCluster


def trimOrthologsMatricesSemiAutomatic(contMethod,test,pValMethod,processType,SPS,outputName,subjectList,geneSource,orthoMode):
    """
    Description
    -----------
    From the results of shared orthologs and statistical significance between cell clusters of subjects, semi-automatically determines clusters that belong to the same cell type and extract their shared orthologs.
    Semi-automatically means that user still can choose which couple of subjects he wants to compare, but does not need clusters to be already annotated.
    This mode requires the use of a SPS (Starting Point Subject).

    Parameters
    ----------
    contMethod
        str, contains method used for building contingency matrix.
    test
        str, chosen statistical test. 
    processType
        str, defines the processing type (manual, semi-automatic or automatic). 

    Returns
    ----------
    None
    """
    subjectAnnot = None 
    clusterList,geneList,startingSpecies,startingLifestage = automaticIdentityStartingPoint(SPS)
    newSp,newLifestage = param.chooseNewSpeciesLifestage(subjectList)
    print(newSp)
    speciesList = [startingSpecies,newSp] ; lifestageList = [startingLifestage,newLifestage]
    listOfUsedSpecies = speciesList.copy()
    correspList = []
    try:
        if contMethod == "orthopairsBased":
            statMatrix = pd.read_csv('results/stat/'+geneSource+'/'+orthoMode+'/'+contMethod+'/'+speciesList[0]+'_'+lifestageList[0]+'_'+speciesList[1]+'_'+lifestageList[1]+"_"+test+"_test_1000g.csv",index_col=0)
        elif contMethod == "genomeBased":
            statMatrix = pd.read_csv('results/stat/'+geneSource+'/'+orthoMode+'/'+contMethod+'/'+speciesList[0]+'_'+lifestageList[0]+'_'+speciesList[1]+'_'+lifestageList[1]+"_"+test+"_test_1000g_"+pValMethod+".csv",index_col=0)
        matrixOrthologs = pd.read_csv('results/id/'+geneSource+'/'+orthoMode+'/'+contMethod+'/'+speciesList[0]+'_'+lifestageList[0]+'_'+speciesList[1]+'_'+lifestageList[1]+"_orthologsId.csv",index_col=0)
    except FileNotFoundError:
            try:
                if contMethod == "orthopairsBased":
                    statMatrix = pd.read_csv('results/stat/'+geneSource+'/'+orthoMode+'/'+contMethod+'/'+speciesList[1]+'_'+lifestageList[1]+'_'+speciesList[0]+'_'+lifestageList[0]+"_"+test+"_test_1000g.csv",index_col=0) 
                elif contMethod == "genomeBased":
                    statMatrix = pd.read_csv('results/stat/'+geneSource+'/'+orthoMode+'/'+contMethod+'/'+speciesList[1]+'_'+lifestageList[1]+'_'+speciesList[0]+'_'+lifestageList[0]+"_"+test+"_test_1000g_"+pValMethod+".csv",index_col=0)  
                matrixOrthologs = pd.read_csv('results/id/'+geneSource+'/'+orthoMode+'/'+contMethod+'/'+speciesList[1]+'_'+lifestageList[1]+'_'+speciesList[0]+'_'+lifestageList[0]+"_orthologsId.csv",index_col=0)
            except FileNotFoundError:
                #print("\nOne of the subject does not exit")
                return
    verifiedCluster = verifyClusterIdentity(statMatrix,clusterList,speciesList,lifestageList,geneList,orthoMode,geneSource)
    clusterListIntegers = list(map(int, clusterList))
    sharedOrthologsList,corresp = extractSharedOrthologs(matrixOrthologs,subjectAnnot,speciesList,lifestageList,processType,clusterListIntegers,verifiedCluster,startingSpecies,startingLifestage)
    correspList.append(corresp)
    dfListSharedOrthologs = []
    sharedOrthologsList = parseSharedOrthologsList(sharedOrthologsList) 
    dfListSharedOrthologs.append(sharedOrthologsList)
    res = ""
    while (res != "y" or res !="n"):
        print("\nDo you want to continue with another subject (y / n)?")
        res = input()
        if res == "n":
            sharedOrthologsList = trimDuplicateOrthologs(sharedOrthologsList)
            fastaFinalDict = manageProteinSequences(listOfUsedSpecies,sharedOrthologsList)
            saveProteinsFasta(fastaFinalDict,outputName)
            return
        elif res == "y":
            sharedOrthologsList,dfListSharedOrthologs,listOfUsedSpecies = continueOrthologsComparison(speciesList,lifestageList,clusterListIntegers,sharedOrthologsList,dfListSharedOrthologs,listOfUsedSpecies,processType,geneList,startingSpecies,startingLifestage,contMethod,test,pValMethod,subjectList)
            sharedOrthologsList = finalOrthologsTrimming(dfListSharedOrthologs,sharedOrthologsList)
            fastaFinalDict = manageProteinSequences(listOfUsedSpecies,sharedOrthologsList)
            saveProteinsFasta(fastaFinalDict,outputName)
            return 

    
def finalOrthologsTrimming(dfListSharedOrthologs,sharedOrthologsList):
    """
    Description
    -----------
    Uses the collection of lists that contain shared orthologs between two subjects to extract the shared orthologs in common between every subject.

    Parameters
    ----------
    dfListSharedOrthologs
        list, contains lists of shared orthologs between two compared subjects.
    sharedOrthologsList
        list, contains all shared orthologs that are in common between every subjects.
    
    Returns
    ----------
    sharedOrthologsListPosttrim
        list, contains all shared orthologs that are in common between every subjects. 
    """
    regexGeneList = [] 
    for i in dfListSharedOrthologs:
        for j in i:
             regexGene = re.split(r"([0-9])",j)[0]
             regexGeneList.append(regexGene)
    regexGeneList = list(set(regexGeneList)) # Get list of identifiable species genes that were used
    regexOccurenceList = []
    for i in regexGeneList:
        countOccurence = 0
        for j in dfListSharedOrthologs:
            for k in j:
                if i in k:
                    countOccurence +=1
                    break
        regexOccurenceList.append(countOccurence)
    test = []
    for i in dfListSharedOrthologs:
        test = test + i
    indexListRemoval = []
    c = 0
    for i in sharedOrthologsList:
        regexTag = re.split(r"([0-9])",i)[0]
        indx = [i for i, s in enumerate(regexGeneList) if regexTag in s][0]
        comparisonNumber = regexOccurenceList[indx]
        geneOccurence = sharedOrthologsList.count(i)
        if comparisonNumber != geneOccurence:
            indexListRemoval.append(c)
        c += 1
    for i in sorted(indexListRemoval, reverse=True):
        del sharedOrthologsList[i]
    sharedOrthologsListPosttrim = list(set(sharedOrthologsList))
    return sharedOrthologsListPosttrim


def fullyAutomaticProcess(contMethod,test,pValMethod,processType,iterationNbr,SPS,outputName,subjectList,geneSource,orthoMode):
    """
    Description
    -----------
    From the results of shared orthologs and statistical significance between cell clusters of subjects, semi-automatically determines clusters that belong to the same cell type and extract their shared orthologs.
    Automatically means that every comparisons between every subjects will be done without user intervention.
    This mode requires the use of a SPS (Starting Point Subject).

    Parameters
    ----------
    contMethod
        str, contains method used for building contingency matrix.
    test
        str, chosen statistical test. 
    processType
        str, defines the processing type (manual, semi-automatic or automatic). 

    Returns
    ----------
    None
    """
    subjectAnnot = None
    startingClusterList,startingGeneList,startingSpecies,startingLifestage = automaticIdentityStartingPoint(SPS)
    speciesList,lifestageList = param.automaticSpeciesLifestage(subjectList)
    dfListSharedOrthologs = []
    sharedOrthologsListConcat = []
    nbrGenesList = []
    nbrOrthologsList = []
    subjectNameBarplot = []
    nbrStartingGenes = getGenesMetrics(startingSpecies,startingLifestage,startingClusterList)
    nbrGenesList.append(nbrStartingGenes)
    subjectNameBarplot.append(startingSpecies+" "+startingLifestage)
    subjectOrderChordDiagram = []
    listSpeciesLifestage = [] # For comparison not involving starting species
    listCluster = [] # For comparison not involving starting species
    correspList = []
    c = 0
    for species in speciesList:
        for lifestage in lifestageList:
            if (species != startingSpecies or lifestage != startingLifestage):
                statMatrix = ""
                comparedSpecies = [startingSpecies,species] ; comparedLifestage = [startingLifestage,lifestage]
                try:
                    if startingSpecies == species:
                        statMatrix = pd.read_csv('results/stat/sameSpecies/'+startingSpecies+'_'+startingLifestage+'_'+species+'_'+lifestage+'_'+test+"_test_1000g.csv",index_col=0)
                        matrixOrthologs = pd.read_csv('results/id/sameSpecies/'+startingSpecies+'_'+startingLifestage+'_'+species+'_'+lifestage+"_orthologsId.csv",index_col=0)
                    else:
                        if contMethod == "orthopairsBased":
                            statMatrix = pd.read_csv('results/stat/'+geneSource+'/'+orthoMode+'/'+contMethod+'/'+startingSpecies+'_'+startingLifestage+'_'+species+'_'+lifestage+'_'+test+"_test_1000g.csv",index_col=0)
                        elif contMethod == "genomeBased":  
                            statMatrix = pd.read_csv('results/stat/'+geneSource+'/'+orthoMode+'/'+contMethod+'/'+startingSpecies+'_'+startingLifestage+'_'+species+'_'+lifestage+'_'+test+"_test_1000g_"+pValMethod+".csv",index_col=0) 
                        matrixOrthologs = pd.read_csv('results/id/'+geneSource+'/'+orthoMode+'/'+contMethod+'/'+startingSpecies+'_'+startingLifestage+'_'+species+'_'+lifestage+"_orthologsId.csv",index_col=0) 
                except FileNotFoundError:
                    try:
                        if startingSpecies == species:
                            statMatrix = pd.read_csv('results/stat/sameSpecies/'+species+'_'+lifestage+'_'+startingSpecies+'_'+startingLifestage+'_'+test+"_test_1000g.csv",index_col=0)
                            matrixOrthologs = pd.read_csv('results/id/sameSpecies/'+species+'_'+lifestage+'_'+startingSpecies+'_'+startingLifestage+"_orthologsId.csv",index_col=0) 
                        else:
                            if contMethod == "orthopairsBased": 
                                statMatrix = pd.read_csv('results/stat/'+geneSource+'/'+orthoMode+'/'+contMethod+'/'+species+'_'+lifestage+'_'+startingSpecies+'_'+startingLifestage+'_'+test+"_test_1000g.csv",index_col=0) 
                            elif contMethod == "genomeBased":
                                statMatrix = pd.read_csv('results/stat/'+geneSource+'/'+orthoMode+'/'+contMethod+'/'+species+'_'+lifestage+'_'+startingSpecies+'_'+startingLifestage+'_'+test+"_test_1000g_"+pValMethod+".csv",index_col=0) 
                            matrixOrthologs = pd.read_csv('results/id/'+geneSource+'/'+orthoMode+'/'+contMethod+'/'+species+'_'+lifestage+'_'+startingSpecies+'_'+startingLifestage+"_orthologsId.csv",index_col=0) 
                    except FileNotFoundError:
                        print("\nSubject "+species+" "+lifestage+" does not exit")
                        pass 
                if len(statMatrix) != 0:
                    verifiedCluster = verifyClusterIdentity(statMatrix,startingClusterList,comparedSpecies,comparedLifestage,startingGeneList,orthoMode,geneSource)
                    clusterListIntegers = list(map(int,startingClusterList))
                    sharedOrthologsList,corresp = extractSharedOrthologs(matrixOrthologs,subjectAnnot,comparedSpecies,comparedLifestage,processType,clusterListIntegers,verifiedCluster,startingSpecies,startingLifestage)
                    correspList.append(corresp)
                    sharedOrthologsList = parseSharedOrthologsList(sharedOrthologsList) 
                    sharedOrthologsListConcat = sharedOrthologsListConcat + sharedOrthologsList
                    dfListSharedOrthologs.append(sharedOrthologsList)
                    listSpeciesLifestage.append(species+"_"+lifestage) ;  listCluster.append(verifiedCluster)
                    nbrGenes = getGenesMetrics(species,lifestage,verifiedCluster)
                    nbrGenesList.append(nbrGenes)
                    subjectNameBarplot.append(species+" "+lifestage)
                    subjectOrder = (comparedSpecies[0]+'_'+comparedLifestage[0],comparedSpecies[1]+'_'+comparedLifestage[1]) 
                    subjectOrderChordDiagram.append(subjectOrder)
                    c+=1
    for i in range(len(listSpeciesLifestage)):
        for j in range(len(listSpeciesLifestage)):
            if listSpeciesLifestage[i] != listSpeciesLifestage[j]:
                statMatrix = ""
                sp1 = listSpeciesLifestage[i].split('_')[0] ; sp2 = listSpeciesLifestage[j].split('_')[0]
                lf1 = listSpeciesLifestage[i].split('_')[1] ; lf2 = listSpeciesLifestage[j].split('_')[1] 
                comparedSpecies = [sp1,sp2] ; comparedLifestage = [lf1,lf2]
                try:
                    if listSpeciesLifestage[i] == listSpeciesLifestage[j]:
                        statMatrix = pd.read_csv('results/stat/sameSpecies/'+listSpeciesLifestage[i]+'_'+listSpeciesLifestage[j]+'_'+test+"_test_1000g.csv",index_col=0)
                        matrixOrthologs = pd.read_csv('results/id/sameSpecies/'+listSpeciesLifestage[i]+'_'+listSpeciesLifestage[j]+"_orthologsId.csv",index_col=0)
                    else:
                        if contMethod == "orthopairsBased":  
                            statMatrix = pd.read_csv('results/stat/'+geneSource+'/'+orthoMode+'/'+contMethod+'/'+listSpeciesLifestage[i]+'_'+listSpeciesLifestage[j]+'_'+test+"_test_1000g.csv",index_col=0)
                        elif contMethod == "genomeBased":
                            statMatrix = pd.read_csv('results/stat/'+geneSource+'/'+orthoMode+'/'+contMethod+'/'+listSpeciesLifestage[i]+'_'+listSpeciesLifestage[j]+'_'+test+"_test_1000g_"+pValMethod+".csv",index_col=0)
                        matrixOrthologs = pd.read_csv('results/id/'+geneSource+'/'+orthoMode+'/'+contMethod+'/'+listSpeciesLifestage[i]+'_'+listSpeciesLifestage[j]+"_orthologsId.csv",index_col=0)
                except FileNotFoundError:
                    pass
                if len(statMatrix) != 0: 
                    sharedOrthologsList,corresp = extractSharedOrthologs(matrixOrthologs,subjectAnnot,comparedSpecies,comparedLifestage,processType,listCluster[i],listCluster[j],startingSpecies,startingLifestage)
                    correspList.append(corresp)
                    sharedOrthologsList = parseSharedOrthologsList(sharedOrthologsList) 
                    sharedOrthologsListConcat = sharedOrthologsListConcat + sharedOrthologsList
                    dfListSharedOrthologs.append(sharedOrthologsList)
                    subjectOrder = (comparedSpecies[0]+'_'+comparedLifestage[0],comparedSpecies[1]+'_'+comparedLifestage[1]) 
                    subjectOrderChordDiagram.append(subjectOrder)
    nbrStartingOrthologs = getOrthologsMetrics(dfListSharedOrthologs,startingSpecies,startingLifestage,startingClusterList)
    nbrOrthologsList.append(nbrStartingOrthologs)
    for species in speciesList:
        for lifestage in lifestageList:
            if (species != startingSpecies or lifestage != startingLifestage):
                statMatrix = ""
                comparedSpecies = [startingSpecies,species] ; comparedLifestage = [startingLifestage,lifestage]
                try:
                    if startingSpecies == species:
                        statMatrix = pd.read_csv('results/stat/sameSpecies/'+startingSpecies+'_'+startingLifestage+'_'+species+'_'+lifestage+'_'+test+"_test_1000g.csv",index_col=0)
                        matrixOrthologs = pd.read_csv('results/id/sameSpecies/'+startingSpecies+'_'+startingLifestage+'_'+species+'_'+lifestage+"_orthologsId.csv",index_col=0) 
                    else:
                        if contMethod == "orthopairsBased":  
                            statMatrix = pd.read_csv('results/stat/'+geneSource+'/'+orthoMode+'/'+contMethod+'/'+startingSpecies+'_'+startingLifestage+'_'+species+'_'+lifestage+'_'+test+"_test_1000g.csv",index_col=0)
                        elif contMethod == "genomeBased":
                            statMatrix = pd.read_csv('results/stat/'+geneSource+'/'+orthoMode+'/'+contMethod+'/'+startingSpecies+'_'+startingLifestage+'_'+species+'_'+lifestage+'_'+test+"_test_1000g_"+pValMethod+".csv",index_col=0) 
                        matrixOrthologs = pd.read_csv('results/id/'+geneSource+'/'+orthoMode+'/'+contMethod+'/'+startingSpecies+'_'+startingLifestage+'_'+species+'_'+lifestage+"_orthologsId.csv",index_col=0)
                except FileNotFoundError:
                    try:
                        if startingSpecies == species:
                            statMatrix = pd.read_csv('results/stat/sameSpecies/'+species+'_'+lifestage+'_'+startingSpecies+'_'+startingLifestage+'_'+test+"_test_1000g.csv",index_col=0)
                            matrixOrthologs = pd.read_csv('results/id/sameSpecies/'+species+'_'+lifestage+'_'+startingSpecies+'_'+startingLifestage+"_orthologsId.csv",index_col=0)
                        else:
                            if contMethod == "orthopairsBased":  
                                statMatrix = pd.read_csv('results/stat/'+geneSource+'/'+orthoMode+'/'+contMethod+'/'+species+'_'+lifestage+'_'+startingSpecies+'_'+startingLifestage+'_'+test+"_test_1000g.csv",index_col=0) 
                            elif contMethod == "genomeBased":
                                statMatrix = pd.read_csv('results/stat/'+geneSource+'/'+orthoMode+'/'+contMethod+'/'+species+'_'+lifestage+'_'+startingSpecies+'_'+startingLifestage+'_'+test+"_test_1000g_"+pValMethod+".csv",index_col=0)  
                            matrixOrthologs = pd.read_csv('results/id/'+geneSource+'/'+orthoMode+'/'+contMethod+'/'+species+'_'+lifestage+'_'+startingSpecies+'_'+startingLifestage+"_orthologsId.csv",index_col=0)
                    except FileNotFoundError:
                        #print("\nSubject "+species+" "+lifestage+" does not exit")
                        pass
                if len(statMatrix) != 0: 
                    verifiedCluster = verifyClusterIdentity(statMatrix,startingClusterList,comparedSpecies,comparedLifestage,startingGeneList,orthoMode,geneSource)
                    nbrOrthologs = getOrthologsMetrics(dfListSharedOrthologs,species,lifestage,verifiedCluster)
                    nbrOrthologsList.append(nbrOrthologs)
    sharedOrthologsListConcat = finalOrthologsTrimming(dfListSharedOrthologs,sharedOrthologsListConcat) 
    print("\nList of pan-cnidarian orthologs for targeted cell type :\n")
    print(sharedOrthologsListConcat)
    print("\nLength of genes list :\n")
    print(len(sharedOrthologsListConcat))
    fastaFinalDict = manageProteinSequences(speciesList,sharedOrthologsListConcat)
    saveProteinsFasta(fastaFinalDict,outputName)
    subjectNameBarplot2 = subjectNameBarplot.copy()
    horizontalBarplotGenes(nbrGenesList,subjectNameBarplot,outputName)
    horizontalBarplotOrthologs(nbrOrthologsList,subjectNameBarplot2,sharedOrthologsListConcat,outputName)
    logOutput(correspList,processType,contMethod,test,pValMethod,iterationNbr,SPS,outputName,dfListSharedOrthologs,subjectOrderChordDiagram)
    pfamAnnotPCO(outputName)
    genesPfamDict = linkPfamGenes(outputName)
    pfamCodeList = retrieveTranscriptionFactorsPfam()
    sequencesTF(outputName,genesPfamDict,pfamCodeList)
    getPfamAnnotTF(outputName)
    chordDiagram(subjectOrderChordDiagram,dfListSharedOrthologs,outputName)
    return 


def horizontalBarplotGenes(nbrGeneList,subjectNameBarplot,outputName):
    """
    Description
    -----------
    Builds and saves horizontal barplot from number of genes belonging to groups, where each group are each subjects and the total of subjects.

    Parameters
    ----------
    nbrGeneList
        list, list of number of genes for each group.
    subjectNameBarplot
        list, list of name of each group.

    Returns
    ----------
    None    
    """
    total = sum(nbrGeneList)
    nbrGeneList.append(total)
    subjectNameBarplot.append('Total genes')    
    df = pd.DataFrame({'Groups':subjectNameBarplot,'Genes':nbrGeneList})
    df = df.sort_values(by=['Genes'])
    fig = plt.figure(figsize=(15, 10))
    plt.barh(y=df.Groups, width=df.Genes)
    plt.title('Genes barplot')
    outputPath = 'results/final/'+outputName
    fig.savefig(outputPath+'/genes_barplot.png')
    plt.close(fig)


def horizontalBarplotOrthologs(nbrOrthologsList,subjectNameBarplot,sharedOrthologsListConcat,outputName):
    """
    Description
    -----------
    Builds and saves horizontal barplot from number of shared orthologs belonging to groups, where each group are each subjects and the total of subjects.

    Parameters
    ----------
    nbrOrthologsList
        list, list of number of shared orthologs for each group.
    subjectNameBarplot
        list, list of name of each group.
    sharedOrthologsListConcat
        list, list of every shared orthologs between every subjects concatenated.

    Returns
    ----------
    None    
    """
    total = sum(nbrOrthologsList)
    nbrOrthologsList.append(total)
    nbrOrthologsList.append(len(sharedOrthologsListConcat))
    subjectNameBarplot.append('Total orthologs')   
    subjectNameBarplot.append('Pan-cnidarian')
    df = pd.DataFrame({'Groups':subjectNameBarplot,'Genes':nbrOrthologsList})
    df = df.sort_values(by=['Genes'])
    colorList=['#595757','#A0DB8E','#ffe599','#B7E4F1','#EEB1C7','#A7A0DA','#E7E7E7']
    fig = plt.figure(figsize=(15, 10))
    plt.barh(y=df.Groups, width=df.Genes,color=colorList,linewidth=0.9,edgecolor='black')
    plt.title('Orthologs barplot')
    outputPath = 'results/final/'+outputName
    fig.savefig(outputPath+'/orthologs_barplot.png')
    plt.close(fig)


def getGenesMetrics(species,lifestage,clusterList): # total genes / genes per subject / total orthologs / orthologs per subject / pancnidarian orthologs
    genesDict = basic.getAllMarkers(species,lifestage)
    nbrGenes = 0
    for i in clusterList:
        if type(i) == str:
            nbrGenes += len(genesDict[int(i)])
        else:
            nbrGenes += len(genesDict[i])
    return nbrGenes


def getOrthologsMetrics(dfListSharedOrthologs,species,lifestage,clusterList):
    """
    Description
    -----------
    Counts number of orthologs belonging to specific cell clusters of a subject in order to builds horizontal barplots. 

    Parameters
    ----------
    dfListSharedOrthologs
        list, contains lists of shared orthologs between two compared subjects.
    species
        str, species name of subject.
    lifestage
        str, lifestage name of subject.
    clusterList
        list, list of number of cell clusters.

    Returns
    ----------
    lenNbrOrthologs
        int, number of orthologs belonging to specific cell clusters of a subject.    
    """
    genesDict = basic.getAllMarkers(species,lifestage)
    genes = []
    for i in clusterList:
        if type(i) == str:
            genes = genes + genesDict[int(i)]
        else:
            genes = genes + genesDict[i]
    allOrthologs = []
    for i in dfListSharedOrthologs:
        allOrthologs = allOrthologs + i
    allOrthologs = list(set(allOrthologs))
    nbrOrthologs = [value for value in genes if value in allOrthologs]
    lenNbrOrthologs = len(nbrOrthologs)
    return lenNbrOrthologs


def chordDiagram(subjectOrderChordDiagram,dfListSharedOrthologs,outputName):
    """
    Description
    -----------
    Builds and saves a chord diagram from the differents lists of shared orthologs between every subjects.

    Parameters
    ----------
    subjectOrderChordDiagram
        list, list of subjects which were used for comparisons.
    dfListSharedOrthologs
        list, contains lists of shared orthologs between two compared subjects.

    Returns
    ----------
    None
    """
    source = []
    target = []
    weights = []
    for i in range(len(subjectOrderChordDiagram)):    
       source.append(subjectOrderChordDiagram[i][0]) 
       target.append(subjectOrderChordDiagram[i][1])
       weights.append(len(dfListSharedOrthologs[i]))
    df = pd.DataFrame(data=np.c_[source, target, weights], columns=['source','target','weight'])
    d3 = D3Blocks(frame=False)
    d3.chord(df, color='#2D3031', opacity='source', cmap='Set2')
    initColor = Color(rgb=[247, 177, 79])
    for i in range(len(subjectOrderChordDiagram)):
        randomColor = Color(hsv=[RANDOM, initColor.hsv.s, 98])
        d3.node_properties.get(subjectOrderChordDiagram[i][0])['color']=randomColor.hex
        #d3.node_properties.get(i)['color']='#ffe599'
        #d3.node_properties.get(i)['color']='#EEB1C7'
        #d3.node_properties.get(i)['color']='#A7A0DA'
        #d3.node_properties.get(i)['color']='#B7E4F1'
    #d3.filepath = "result/final/"+outputName+"/chordDiagram.html"
    #img = urllib.request.urlretrieve(example_url, "PHI.svg")
    
    d3.show()


def logOutput(correspList,processType,contMethod,test,pValMethod,iterationNbr,SPS,outputName,dfListSharedOrthologs,subjectOrderChordDiagram):
    with open("results/final/"+outputName+"/log_results.txt","w") as f:
        f.write("Parameters :\n") 
        f.write("Processing : "+processType+" | Contingency matrix method : "+contMethod+" | Statistical test : "+test+" | Combining p-value method : "+pValMethod+" | Starting point subject : "+SPS+" | Number of iteration : "+str(iterationNbr)+"\n")
        for i in correspList:
            if i != None:
                f.write(i)
        f.write("\n\nShared orthologs between subjects:\n")
        for i in range(len(dfListSharedOrthologs)):
            sp1 = subjectOrderChordDiagram[i][0].split('_')[0] ; lf1 = subjectOrderChordDiagram[i][0].split('_')[1] 
            sp2 = subjectOrderChordDiagram[i][1].split('_')[0] ; lf2 = subjectOrderChordDiagram[i][1].split('_')[1] 
            f.write("\nSubjects compared : "+sp1 +" "+lf1+" & "+sp2+" "+lf2)
            if sp1 != sp2:
                f.write("\nShared orthologs :\n"+', '.join(dfListSharedOrthologs[i])+"\n")
            else:
                f.write("\nShared genes :\n"+', '.join(dfListSharedOrthologs[i])+"\n") 
        f.close()


def retrieveTranscriptionFactorsPfam():
    pfamCodeList = []
    with open("input/transcriptionFactors.txt","r") as f:
        l = f.readline()
        while l != "":
            if l[0] != "#":
                pfamCode = l.split(" ")[1]
                pfamCode = pfamCode.split(" ")[0]
                if "PF" in pfamCode:
                    pfamCodeList.append(pfamCode)
            l = f.readline()
    f.close()
    return pfamCodeList
            

def pfamAnnotPCO(outputName):
    if os.path.exists("results/final/"+outputName+"/shared_orthologs.fasta") == False:
        print("PCO were not retrieved, impossible to process Pfam annotation.\n")
    else:
        if os.path.exists("results/fina/"+outputName+"/shared_orthologs_pfam") == False: 
            print("Currently processing Pfam annotation of PCO.\n")
            subprocess.call(['hmmsearch','--cut_ga','-o',"results/final/"+outputName+"/shared_orthologs_pfam",'--domtblout',"results/final/"+outputName+"/shared_orthologs_dt_out",
                     "input/pfam/Pfam-A.hmm","results/final/"+outputName+"/shared_orthologs.fasta"])
        print("Pfam annotation of PCO has been processed.\n")
    return

def linkPfamGenes(outputName):
    genesPfamDict = {}
    with open("results/final/"+outputName+"/shared_orthologs_dt_out","r") as f:
        l = f.readline()
        while l != "":
            if l[0] != "#":
                tempList = l.split(" ")
                strippedList = list(filter(None, tempList))
                gene = strippedList[0]
                pfam = strippedList[4].split('.')[0]
                genesPfamDict[gene] = pfam
            l = f.readline()
    f.close()
    return genesPfamDict


def sequencesTF(outputName,genesPfamDict,pfamCodeList):
    tfDict = {}
    with open("results/final/"+outputName+"/shared_orthologs.fasta","r") as f:
        l = f.readline()
        seqFlag = False
        prot = ""
        while l != "":
            if l[0] == ">":
                gene = l.split(">")[1]
                gene = gene.replace("\n","")
                if gene in genesPfamDict:
                    pfam = genesPfamDict[gene]
                    if pfam in pfamCodeList:
                        seqFlag = True
            else:
                if seqFlag == True:
                    while l[0] != ">":
                        prot += l
                        l = f.readline()
                    if "*" in prot:
                        prot = prot.replace("*","")
                    tfDict[gene] = prot
                    prot = ""
                    seqFlag = False
            l = f.readline()
    with open("results/final/"+outputName+"/transcription_factors.fasta", 'w') as f:
        for k,v in tfDict.items():
            f.write('>' + k + '\n')
            f.write(v)
    f.close()
    return


def getPfamAnnotTF(outputName):
    tfList = []
    with open("results/final/"+outputName+"/transcription_factors.fasta","r") as f:
        l = f.readline()
        while l != "":
            if l[0] == ">":
                gene = l.replace(">","")
                gene = gene.replace("\n","")
                tfList.append(gene)
            l = f.readline()
    f.close()
    writeList = []
    with open("results/final/"+outputName+"/shared_orthologs_dt_out","r") as f:
        l = f.readline()
        while l != "":
            if l[0] != "#":
                tempList = l.split(" ")
                strippedList = list(filter(None, tempList))
                gene = strippedList[0]
                if gene in tfList:
                    writeList.append(l)
            else:
                writeList.append(l)
            l = f.readline()
    f.close()
    with open("results/final/"+outputName+"/transcription_factors_dt_out","w") as f:
        f.writelines(writeList)
    f.close()
    return



