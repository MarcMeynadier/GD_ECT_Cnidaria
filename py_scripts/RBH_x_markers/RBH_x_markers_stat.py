#-----------------------------------------------------------------------------------------#
#                       Importation of external libraries and modules                     #
#-----------------------------------------------------------------------------------------#

import pandas as pd
import statistics
import math
import numpy as np
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
from functools import reduce
import scipy.stats as stats
import RBH_x_markers_basic as basic
import os

def contingencyMatrices(RBHtupleList,markersSubjectList,spList,lifestageList,contMethod,orthoMode,geneSource):
    """
    Description
    -----------
    Builds contingency matrices for statistical test.

    Parameters
    ----------
    RBHgenesTupleList
        list, contains tuples of one-to-one orthologs between two species.
    markersDictSpList
        list, contains dictionnary for each species where keys of the dictionnaries are cluster number and values are list of marker genes. 
    spList
        list, contains species of subjects for comparison.
    lifestageList
        list, contains lifestages of subjects for comparison. 
    contMethod
        str, contains method used for building contingency matrix.

    Returns
    ----------
    None
    """ 
    if RBHtupleList == "sameSpecies":
        resPath = 'results/stat/sameSpecies'
    else:
        if geneSource == "RBH":
            resPath = 'results/stat/RBH/'+orthoMode+'/'+contMethod 
        elif geneSource == "orthofinder":
            resPath = 'results/stat/orthofinder/'+orthoMode+'/'+contMethod 
    if not os.path.exists(resPath):
        os.makedirs(resPath)
    if RBHtupleList == "sameSpecies":
        species = None
        matrixList = basic.genesCounting(RBHtupleList,markersSubjectList,orthoMode,contMethod,species,spList)
        with open(resPath+'/'+spList[0]+'_'+lifestageList[0]+'_'+spList[1]+'_'+lifestageList[1]+'_markers_matrix_list_1000g.csv', "wb") as f:   
            pickle.dump(matrixList, f)
    else:
        if contMethod == "orthopairsBased":
            species = None
            matrixList = basic.genesCounting(RBHtupleList,markersSubjectList,orthoMode,contMethod,species,spList)
            with open(resPath+'/'+spList[0]+'_'+lifestageList[0]+'_'+spList[1]+'_'+lifestageList[1]+'_markers_matrix_list_1000g.csv', "wb") as f:   
                pickle.dump(matrixList, f)
        elif contMethod == "genomeBased":
            matrixListSp1 = basic.genesCounting(RBHtupleList,markersSubjectList,orthoMode,contMethod,spList[0],spList)
            with open(resPath+'/'+spList[0]+'_'+lifestageList[0]+'_'+spList[1]+'_'+lifestageList[1]+'_markers_matrix_list_1000g_'+spList[0]+'.csv', "wb") as f:   
                pickle.dump(matrixListSp1, f)
            matrixListSp2 = basic.genesCounting(RBHtupleList,markersSubjectList,orthoMode,contMethod,spList[1],spList)  
            with open(resPath+'/'+spList[0]+'_'+lifestageList[0]+'_'+spList[1]+'_'+lifestageList[1]+'_markers_matrix_list_1000g_'+spList[1]+'.csv', "wb") as f:   
                pickle.dump(matrixListSp2, f) 
    print("Contingency matrix processed for "+spList[0]+' '+lifestageList[0]+' & '+spList[1]+' '+lifestageList[1]) 
    return


def combinePvalues(matrixSp1,matrixSp2,pValMethod):
    """
    Description
    -----------
    Takes two p-values in input to give one p-value in output.

    Parameters
    ----------
    matrixSp1,
        pandas.core.frame.DataFrame, dataframe of p-values between clusters of species 1 and species 2, with contingency matrix builds by using genomeBased method with genome of species 1.
    matrixSp2,
        pandas.core.frame.DataFrame, dataframe of p-values between clusters of species 1 and species 2, with contingency matrix builds by using genomeBased method with genome of species 2.

    Returns
    ----------
    outputMatrix,
        pandas.core.frame.DataFrame, dataframe of combined p-values between the two dataframes in input.
    """
    matrixSp2T = matrixSp2.transpose()
    dfSp1 = matrixSp1.values.tolist()
    dfSp2 = matrixSp2T.values.tolist()
    outputMatrix = []
    for i in range(len(dfSp1)):
        outputRow = []
        for j in range(len(dfSp1[i])):
            pval = stats.combine_pvalues([dfSp1[i][j],dfSp2[i][j]],method=pValMethod)
            pvalReturn = pval[1]
            if pvalReturn == 0.0:
                pvalReturn = 1e-293 
            outputRow.append(pvalReturn)
        outputMatrix.append(outputRow)
    outputMatrix = pd.DataFrame(outputMatrix)
    return outputMatrix


def inferentialTest(RBHtupleList,speciesList,lifestageList,contMethod,test,pValMethod,orthoMode,geneSource):
    """
    Description
    -----------
    Calculates p-values between clusters from contingency matrix by applying statistical test chosen in parameter.

    Parameters
    ----------
    speciesList
        list, contains species of subjects for comparison.
    lifestageList
        list, contains lifestages of subjects for comparison. 
    contMethod
        str, contains method used for building contingency matrix.
    test
        str, chosen statistical test.

    Returns
    ----------
    None
    """
    listMatrices = []
    testMatrixList = []
    if RBHtupleList != "sameSpecies":
        resPath = 'results/stat/'+geneSource+'/'+orthoMode+'/'+contMethod+'/' 
    else:
        resPath = 'results/stat/sameSpecies/'
    if contMethod == "orthopairsBased" or RBHtupleList == "sameSpecies":
        testMatrix = []
        with open(resPath+speciesList[0]+'_'+lifestageList[0]+'_'+speciesList[1]+'_'+lifestageList[1]+'_markers_matrix_list_1000g.csv', "rb") as f:   
            matrixList = pickle.load(f)
        listMatrices.append(matrixList)
        for i in range(len(matrixList[0])):
            listZero = []
            for j in range(len(matrixList[0][i])):
                listZero.append(0)
            testMatrix.append(listZero)
        testMatrixList.append(testMatrix)
    elif contMethod == "genomeBased":
        with open(resPath+speciesList[0]+'_'+lifestageList[0]+'_'+speciesList[1]+'_'+lifestageList[1]+'_markers_matrix_list_1000g_'+speciesList[0]+'.csv', "rb") as f:   
            matrixList = pickle.load(f)
        with open(resPath+speciesList[0]+'_'+lifestageList[0]+'_'+speciesList[1]+'_'+lifestageList[1]+'_markers_matrix_list_1000g_'+speciesList[1]+'.csv', "rb") as f:   
            matrixList2 = pickle.load(f) 
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
                if test == "Fisher":
                    pStat = stats.fisher_exact(testDf,alternative='greater')
                    pvalReturn = pStat[1]
                    if pvalReturn == 0.0:
                        pvalReturn = 1e-293 
                    testMatrixList[m][i][j] = pvalReturn
                elif test == "Barnard":
                    pStat = stats.barnard_exact(testDf)
                    pvalReturn = pStat.pvalue
                    if pvalReturn == 0.0:
                        pvalReturn = 1e-293 
                    testMatrixList[m][i][j] = pvalReturn
                elif test =="Chi2":
                    pStat = stats.chi2_contingency(testDf)
                    pvalReturn = pStat[1]
                    if pvalReturn == 0.0:
                        pvalReturn = 1e-293 
                    testMatrixList[m][i][j] = pvalReturn
                elif test == "Boschloo":
                    pStat = stats.boschloo_exact(testDf)
                    pvalReturn = pStat.pvalue
                    if pvalReturn == 0.0:
                        pvalReturn = 1e-293 
                    testMatrixList[m][i][j] = pvalReturn
        outputTestMatrix = pd.DataFrame(testMatrixList[m])
        outputMatricesList.append(outputTestMatrix)
    if contMethod == "orthopairsBased":
        resultMatrix = outputMatricesList[0]
    elif contMethod == "genomeBased":
        if speciesList[0] != speciesList[1]:
            resultMatrix = combinePvalues(outputMatricesList[0],outputMatricesList[1],pValMethod)
        else:
            resultMatrix = outputMatricesList[0] 
    resultMatrix = resultMatrix.values.tolist()
    flattenMatrix = reduce(lambda a, b: a + b, resultMatrix) 
    pAdjustedMatrix = sm.stats.multipletests(pvals=flattenMatrix,method="bonferroni")
    colLen = len(resultMatrix[0])
    rowLen = len(resultMatrix)
    pAdjustedMatrix = pd.DataFrame(np.array(pAdjustedMatrix[1]).reshape(rowLen,colLen))
    if contMethod == "genomeBased":
        pAdjustedMatrix.to_csv(resPath+speciesList[0]+'_'+lifestageList[0]+'_'+speciesList[1]+'_'+lifestageList[1]+'_'+test+"_test_1000g_"+pValMethod+".csv")   
    elif contMethod == "orthopairsBased" or RBHtupleList == "sameSpecies":
        pAdjustedMatrix.to_csv(resPath+speciesList[0]+'_'+lifestageList[0]+'_'+speciesList[1]+'_'+lifestageList[1]+'_'+test+"_test_1000g.csv")
    heatmapGeneration(RBHtupleList,speciesList,lifestageList,contMethod,test,pValMethod,orthoMode,geneSource)
    print("Inferential test processed for "+speciesList[0]+' '+lifestageList[0]+' & '+speciesList[1]+' '+lifestageList[1])
    return 

# Bootstrap test

def basicStats(matrix):
    """
    Description
    -----------
    Calculates means and standard deviation 

    Parameters
    ----------
    matrixSp1,
        pandas.core.frame.DataFrame, dataframe of p-values between clusters of species 1 and species 2, with contingency matrix builds by using genomeBased method with genome of species 1.

    matrixSp2,
        pandas.core.frame.DataFrame, dataframe of p-values between clusters of species 1 and species 2, with contingency matrix builds by using genomeBased method with genome of species 2.

    Returns
    ----------
    outputMatrix,
        pandas.core.frame.DataFrame, dataframe of combined p-values between the two dataframes in input.
    """ 
    headersList = matrix.columns.tolist()
    valuesList = []
    for i in headersList:
        columnList = matrix[i].tolist()
        for j in columnList:
            valuesList.append(j)
    mean = statistics.mean(valuesList)
    SD = statistics.pstdev(valuesList)
    return mean,SD

def zTest(meanSample,meanPop,popSD,sqrtNbrPop):
    """
    Description
    -----------
    Calculates Z-test to test the mean of the distribution. 

    Parameters
    ----------
    meanSample
        float, mean of sample 
    meanPop
        float, mean of population.
    popSD
        float, standard deviation of population.
    sqrtNbrPop
        float, square root of number of values in population.

    Returns
    ----------
    test,
        float, p-value of the Z-test.
    """ 
    test = (meanSample - meanPop) / (popSD/sqrtNbrPop)
    return test

def randomMatricesGenerator(markersDictList,RBHdfTupleList,reciprocality):
    """
    Description
    -----------
    Calls shuffleTuple() to shuffle the pairs of one-to-one orthologs to calculate a random contingency matrix. 

    Parameters
    ----------
    markersDictSpList
        list, contains dictionnary for each species where keys of the dictionnaries are cluster number and values are list of marker genes.
    RBHgenesTupleList
        list, contains tuples of one-to-one orthologs between two species.
    reciprocality,
        str, determines which case of the contingency matrix is calculated. 
    
    Returns
    ----------
    randomMatrix,
       pandas.core.frame.DataFrame, dataframe that contains random values of one of the case of contingency matrix depending on the reciprocality value. 
    """ 
    shuffleTupleList = basic.shuffleTuple(RBHdfTupleList)
    randomMatrix = basic.countOrthologsPairs(shuffleTupleList,markersDictList,reciprocality)
    return randomMatrix

def bootstrapTest(markersDictList,RBHdfTupleList,speciesList,lifestageList,iterationNbr,orthoMode):
    """
    Description
    -----------
    Applies an achieved significance level (ASL) bootstrap on contingency matrix for testing statistical significancy between expression profiles of cell clusters. 

    Parameters
    ----------
    markersDictSpList
        list, contains dictionnary for each species where keys of the dictionnaries are cluster number and values are list of marker genes.
    RBHgenesTupleList
        list, contains tuples of one-to-one orthologs between two species. 
    spList
        list, contains species of subjects for comparison.
    lifestageList
        list, contains lifestages of subjects for comparison.  
    iterationNbr
        int, number of iteration for bootstrap test.

    Returns
    ----------
    None
    """ 
    with open('results/stat/'+orthoMode+'/'+speciesList[0]+'_'+lifestageList[0]+'_'+speciesList[1]+'_'+lifestageList[1]+'_markers_matrix_list_1000g_orthopairsBased.csv', "rb") as f:   
        matrixList = pickle.load(f)
    originalMatrix = matrixList[0]
    originalMatrix = pd.DataFrame(originalMatrix)
    sqrtNbrPopOriginal = math.sqrt(originalMatrix.to_numpy().sum())
    originalMatrixMean,originalMatrixSD = basicStats(originalMatrix)
    originalMatrixTest = []
    for i in originalMatrix.iterrows():
        testList = []
        for j in i[1].tolist():
            test = zTest(j,originalMatrixMean,originalMatrixSD,sqrtNbrPopOriginal)
            testList.append(test)
        originalMatrixTest.append(testList)
    originalMatrixTest = pd.DataFrame(originalMatrixTest) 
    simMatricesTestList = []
    if iterationNbr == None:
        print("\nPlease choose a number of iteration for ASL bootstrap test")
        iterationNbr = int(input())
    for i in range(iterationNbr):
        simMatrix = randomMatricesGenerator(markersDictList,RBHdfTupleList,"yes")
        sqrtNbrPopSim = math.sqrt(simMatrix.to_numpy().sum())
        simMatrixMean,simMatrixSD = basicStats(simMatrix)
        simMatrixTest = []
        for j in simMatrix.iterrows():
            testList = []
            for k in j[1].tolist():
                test = zTest(k,simMatrixMean,simMatrixSD,sqrtNbrPopSim)
                testList.append(test)
            simMatrixTest.append(testList)
        simMatrixTest = pd.DataFrame(simMatrixTest)
        simMatricesTestList.append(simMatrixTest)
    ASLmatrix = originalMatrix.copy()
    for col in ASLmatrix.columns:
        ASLmatrix[col].values[:] = 0
    for i in simMatricesTestList:
        diffDf = originalMatrixTest.subtract(i)
        rowIndex = 0
        for j in diffDf.iterrows():
            colIndex = 0
            for k in j[1].tolist():
                if k <= 0:
                    ASLmatrix.loc[rowIndex,colIndex] += 1
                colIndex += 1
            rowIndex += 1
    ASLmatrix = ASLmatrix.div(iterationNbr)
    ASLmatrix.to_csv('results/stat/'+speciesList[0]+'_'+lifestageList[0]+'_'+speciesList[1]+'_'+lifestageList[1]+"_ASL_bootstrap_"+str(iterationNbr)+"_iterations.csv")
    return


def heatmapGeneration(RBHtupleList,listSp,listLifestage,contMethod,test,pValMethod,orthoMode,geneSource):
    if RBHtupleList == "sameSpecies":
        resPath = 'results/stat/sameSpecies/heatmap/'
    else:
        resPath = 'results/stat/'+geneSource+'/'+orthoMode+'/'+contMethod+'/heatmap/'
    if not os.path.exists(resPath):
        os.makedirs(resPath)
    if RBHtupleList == "sameSpecies":
        matrix = pd.read_csv('results/stat/sameSpecies/'+listSp[0]+"_"+listLifestage[0]+"_"+listSp[1]+"_"+listLifestage[1]+"_"+test+"_test_1000g.csv").drop(columns=['Unnamed: 0'])
    else:
        if contMethod == "orthopairsBased":
            matrix = pd.read_csv("results/stat/"+geneSource+'/'+orthoMode+'/'+contMethod+'/'+listSp[0]+"_"+listLifestage[0]+"_"+listSp[1]+"_"+listLifestage[1]+"_"+test+"_test_1000g.csv").drop(columns=['Unnamed: 0'])
        elif contMethod == "genomeBased":
            matrix = pd.read_csv('results/stat/'+geneSource+'/'+orthoMode+'/'+contMethod+'/'+listSp[0]+'_'+listLifestage[0]+'_'+listSp[1]+'_'+listLifestage[1]+'_'+test+"_test_1000g_"+pValMethod+".csv").drop(columns=['Unnamed: 0']) 
    matrix = np.log10(matrix) * -1 
    matrix = pd.DataFrame.transpose(matrix)
    plt.figure(figsize=(11,9))
    sns.heatmap(data=matrix,cmap="coolwarm",square=True,linewidths=0.3,linecolor="black",xticklabels=True,yticklabels=True,cbar_kws={"shrink": 0.5})
    plt.title('Cell clusters equivalence between '+listSp[0]+" "+listLifestage[0]+" and "+listSp[1]+" "+listLifestage[1], weight='bold',pad=20)
    plt.xlabel('Cell clusters of '+listSp[0]+" "+listLifestage[0],labelpad=20)
    plt.ylabel('Cell clusters of '+listSp[1]+" "+listLifestage[1],labelpad=20)
    plt.savefig(resPath+listSp[0]+"_"+listLifestage[0]+"_"+listSp[1]+"_"+listLifestage[1]+"_"+test+"_test_1000g.png") 