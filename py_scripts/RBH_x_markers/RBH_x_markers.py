"""
RBH_x_markers

Five modules are available in this program :
1 : Main module that user can use to navigate across the program.
2 : Parameters module that contains functions relative to menu choice and display, as well as parameters settings.
3 : Basic module that contains functions useful to both statistical and functional module.
4 : Statistical module that contains functions relative to statistical tests to assess expression proximity between cell clusters.
5 : Functional module that contains functions useful to retrieve the identity of orthologs shared between cell clusters.

Marc Meynadier
"""

#-----------------------------------------------------------------------------------------#
#                       Importation of external libraries and modules                     #
#-----------------------------------------------------------------------------------------#

import sys
import argparse
from os.path import exists
import RBH_x_markers_stat as stat
import RBH_x_markers_identity as id
import RBH_x_markers_basic as basic
import RBH_x_markers_parameters as param


def autoFullRun(args):
    """
    Description
    -----------
    Launches the program in autorun to process everything. 

    Parameters
    ----------
    args
        Command-line arguments

    Returns
    -------
    None
    """
    contMethod = args.contMethod ; test = args.test ; pValMethod = args.pValMethod ; processType = args.processType ; iterationNbr = args.iterationNbr ; SPS = args.SPS ; outputName = args.outputName ; subjectList = args.subjectList
    if SPS == None:
        print("\nStarting point subject is mandatory for automatic processing")
        SPS = param.selectStartingPointSubject()
    spList,lifestageList = param.automaticSpeciesLifestage(subjectList)
    unique_combinations = []
    for i in range(len(spList)):
        for j in range(len(lifestageList)):
            unique_combinations.append((spList[i],lifestageList[j]))
    used_combinations = []
    not_existing = []
    for i in unique_combinations:
        used_combinations.append(i)
        for j in unique_combinations: 
            if j not in used_combinations:
                skipFlag = False
                sp1 = i[0] ; sp2 = j[0]
                lf1 = i[1] ; lf2 = j[1]
                try:
                    markersSubject1 = basic.getAllMarkers(sp1,lf1)
                except FileNotFoundError:
                    not_existing.append((sp1,lf1))
                    skipFlag = True
                    #break
                try:
                    markersSubject2 = basic.getAllMarkers(sp2,lf2)
                except FileNotFoundError:
                    skipFlag = True
                    not_existing.append((sp2,lf2)) 
                    #break
                print("\nCurrently processing "+sp1+" "+lf1+" & "+sp2+" "+lf2+"\n")
                if skipFlag == False:
                    RBHdfTupleList = basic.getRBHtupleList([sp1,sp2])
                    if contMethod == "orthopairsBased":
                        if not exists('results/stat/'+sp1+'_'+lf1+'_'+sp2+'_'+lf2+'_markers_matrix_list_1000g_'+contMethod+'.csv'):
                            stat.contingencyMatrices(RBHdfTupleList,[markersSubject1,markersSubject2],[sp1,sp2],[lf1,lf2],contMethod)  
                    elif contMethod == "genomeBased":
                        if not exists('results/stat/'+sp1+'_'+lf1+'_'+sp2+'_'+lf2+'_markers_matrix_list_1000g_'+contMethod+'_'+sp1+'.csv'):
                            stat.contingencyMatrices(RBHdfTupleList,[markersSubject1,markersSubject2],[sp1,sp2],[lf1,lf2],contMethod)   
                    if contMethod == "orthopairsBased": 
                        if not exists('results/stat/'+sp1+'_'+lf1+'_'+sp2+'_'+lf2+'_'+test+"_test_1000g_"+contMethod+".csv"):
                            stat.inferentialTest([sp1,sp2],[lf1,lf2],contMethod,test,pValMethod)
                    elif contMethod == "genomeBased":
                        if not exists('results/stat/'+sp1+'_'+lf1+'_'+sp2+'_'+lf2+'_'+test+"_test_1000g_"+contMethod+"_"+pValMethod+".csv"):
                            stat.inferentialTest([sp1,sp2],[lf1,lf2],contMethod,test,pValMethod) 
                    if not exists('results/id/'+sp1+'_'+lf1+'_'+sp2+'_'+lf2+"_orthologsId.csv"):  
                        id.getOrthologsPairs(RBHdfTupleList,[markersSubject1,markersSubject2],[sp1,sp2],[lf1,lf2])
    id.fullyAutomaticProcess(contMethod,test,pValMethod,processType,iterationNbr,SPS,outputName,subjectList)
                    

def main(args):
    """
    Description
    -----------
    Launches the program and allows the user to navigate across the different modules. 

    Parameters
    ----------
    args
        Command-line arguments

    Returns
    -------
    None
    """ 
    contMethod = args.contMethod ; test = args.test ; pValMethod = args.pValMethod ; processType = args.processType ; iterationNbr = args.iterationNbr ; SPS = args.SPS ; outputName = args.outputName ; subjectList = args.subjectList
    spList = [] ; lifestageList = []
    if processType == "f":
        autoFullRun(args)
        exit()
    while True:
        param.menuDisplay()
        while True:
            try:
                answer = int(input())
            except ValueError:
                print("\nYou must indicate an integer value ranging from 1 to 4\n")
                continue
            break
        if answer == 1:
            spList,lifestageList,contMethod,test,pValMethod,iterationNbr,SPS,processType = param.menuParameters(spList,lifestageList,contMethod,test,pValMethod,iterationNbr,SPS,processType,subjectList)
            if len(spList) != 0:
                RBHdfTupleList = basic.getRBHtupleList(spList)
        elif answer == 2:
            markersSubjectList = []
            try:
                markersSubject1 = basic.getAllMarkers(spList[0],lifestageList[0])
                markersSubject2 = basic.getAllMarkers(spList[1],lifestageList[1])
            except (UnboundLocalError,IndexError):
                print("\nChoosing species and lifestages is a mandatory step")
                spList,lifestageList,contMethod,test,pValMethod,iterationNbr,SPS,processType = param.menuParameters(spList,lifestageList,contMethod,test,pValMethod,iterationNbr,SPS,processType,subjectList)
                RBHdfTupleList = basic.getRBHtupleList(spList)
                markersSubject1 = basic.getAllMarkers(spList[0],lifestageList[0])
                markersSubject2 = basic.getAllMarkers(spList[1],lifestageList[1])
            markersSubjectList.append(markersSubject1)
            markersSubjectList.append(markersSubject2)
            param.menuStatistics(RBHdfTupleList,markersSubjectList,spList,lifestageList,contMethod,test,pValMethod,iterationNbr)     
        elif answer == 3:
            markersSubjectList = []
            if processType == "m":
                try: 
                    markersSubject1 = basic.getAllMarkers(spList[0],lifestageList[0])
                    markersSubject2 = basic.getAllMarkers(spList[1],lifestageList[1])
                except (UnboundLocalError,IndexError):
                    print("\nChoosing species and lifestages is a mandatory step")
                    spList,lifestageList,contMethod,test,pValMethod,iterationNbr,SPS,processType = param.menuParameters(spList,lifestageList,contMethod,test,pValMethod,iterationNbr,SPS,processType,subjectList)
                    RBHdfTupleList = basic.getRBHtupleList(spList)
                    markersSubject1 = basic.getAllMarkers(spList[0],lifestageList[0])
                    markersSubject2 = basic.getAllMarkers(spList[1],lifestageList[1])
                markersSubjectList.append(markersSubject1)
                markersSubjectList.append(markersSubject2)
                param.menuIdentity(RBHdfTupleList,markersSubjectList,spList,lifestageList,contMethod,test,pValMethod,processType,iterationNbr,SPS,outputName,subjectList)
            elif processType == "a" or processType=="sa":
                RBHdfTupleList = None
                param.menuIdentity(RBHdfTupleList,markersSubjectList,spList,lifestageList,contMethod,test,pValMethod,processType,iterationNbr,SPS,outputName,subjectList)
        elif answer == 4:
            sys.exit(0)
    

# Command line args

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Equivalencing cell types between species & lifestages")
    parser.add_argument('-p','--process',action='store',dest='processType',help='Processing mode - f (full-automatic) / a (automatic) / sa (semi-automatic) / m (manual)') 
    parser.add_argument('-l','--list',action='store',dest='subjectList',help='Text file species and lifestages - speciesLifestage.txt')
    parser.add_argument('-c','--contingency',action='store',dest='contMethod',help='Method for contingency table in statistical module - orthopairsBased / genomeBased')
    parser.add_argument('-t','--test',action='store',dest='test',help='Statistical test in statistical module - Fisher / Barnard / Chi2 / Boschloo')
    parser.add_argument('-v','--pval',action='store',dest='pValMethod',help='Statistical method for combining p-values - tippett / fisher / pearson / stouffer / mudholkar_george')
    parser.add_argument('-i','--iteration',action='store',dest='iterationNbr',help='Iteration number for bootstrap method')
    parser.add_argument('-s','--starting',action='store',dest='SPS',help='Text file of starting point subject - e.g. automaticIdentityCnidocytes.txt')
    parser.add_argument('-o','--output',action='store',dest='outputName',help='Folder name of output results')

    args = parser.parse_args()

    if args.contMethod == None:
        args.contMethod = "orthopairsBased"
    if args.test == None:
        args.test = "Fisher"
    if args.pValMethod == None:
        args.pValMethod = "tippett"
    if args.processType == None:
        args.processType = "m"
    if args.processType == "m":
        args.SPS = None
    if args.subjectList == None:
        args.subjectList = param.selectSubjectList()

print(args)
main(args)