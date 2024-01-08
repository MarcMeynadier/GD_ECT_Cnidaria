import os

# Module importation
import RBH_x_markers_stat as stat
import RBH_x_markers_identity as id

def menuParametersDisplay():
    """
    Description
    -----------
    Displays the parameters module menu.

    Parameters
    ----------
    None
    
    Returns
    -------
    None
    """ 
        
    print("\n")
    print("----------------------------------------------------")
    print("|                                                  |")
    print("|                Parameters menu                   |")
    print("|                                                  |")
    print("|                                                  |") 
    print("|             Species setting : 1                  |")
    print("|                                                  |")
    print("|             Lifestage setting : 2                |")
    print("|                                                  |") 
    print("|             Contingency table method : 3         |")
    print("|                                                  |")
    print("|             Statistical test setting : 4         |")
    print("|                                                  |")
    print("|             Combine p-value method setting : 5   |")
    print("|                                                  |")
    print("|             Iteration number setting : 6         |")
    print("|                                                  |")
    print("|             Process type identity setting : 7    |")
    print("|                                                  |")
    print("|             Help / Guideline : 8                 |") 
    print("|                                                  |")
    print("|             Return : 9                           |")
    print("|                                                  |")
    print("----------------------------------------------------")
    print("\n")
    return 

def menuParameters(spList,lifestageList,contMethod,test,pValMethod,iterationNbr,SPS,processType,subjectList):
    """
    Description
    -----------
    Displays the parameters module menu while letting the user define the different parameters.

    Parameters
    ----------
    spList
        list, empty at first call, designed to contain species of subjects for comparison.
    lifestageList
        list, empty at first call, designed to contain lifestages of subjects for comparison.
    contMethod
        str, contains method used for building contingency matrix.
    test
        str, chosen statistical test. 
    pValMethod
        str, chosen method for combining p-values in genomeBased method.
    iterationNbr
        int, number of iteration for bootstrap test.
    processType
        str, defines the processing type (manual, semi-automatic or automatic).

    Returns
    ----------
    spList
        list, contains species of subjects for comparison.
    lifestageList
        list, contains lifestages of subjects for comparison. 
    contMethod
        str, contains method used for building contingency matrix.
    test
        str, chosen statistical test.
    pValMethod
        str, chosen method for combining p-values in genomeBased method. 
    iterationNbr
        int, number of iteration for bootstrap test.
    processType
        str, defines the processing type (manual, semi-automatic or automatic).
    """
        
    while True:
        menuParametersDisplay()
        while True:
            try:
                answer = int(input())
            except ValueError:
                print("\nYou must indicate an integer value ranging from 1 to 9.\n")
                continue
            break
        if answer == 1:
            sp1,sp2 = chooseParametersSpecies(subjectList)
            spList = [sp1,sp2]
        elif answer == 2:
            sp1Lifestage,sp2Lifestage = chooseParametersLifestage(subjectList)
            lifestageList = [sp1Lifestage,sp2Lifestage]
        elif answer == 3:
            contMethod = chooseContingencyMethod()
        elif answer == 4:
            test = chooseParametersTest()
        elif answer == 5:
            pValMethod = choosePValMethod()
        elif answer == 6:
            iterationNbr = chooseParametersIteration()
        elif answer == 7:
            processType = processTypeId()
            if (processType == "a" or processType == "sa"):
                SPS = selectStartingPointSubject()
        elif answer == 8:
            guideline() 
        elif answer == 9:
            if len(spList) == 0:
                answer2=""
                while (answer2 != "y" or answer2 !="n"):
                    print("\nChoosing species is a mandatory step for statistical and manual identity modules.")
                    print("Are you sure you want to skip this step ? (y / n)")
                    answer2=input()
                    if answer2=="y":
                        spList = [] 
                    elif answer2=="n":    
                        sp1,sp2 = chooseParametersSpecies(subjectList)
                        spList = [sp1,sp2]
                    print(spList,lifestageList,contMethod,test,pValMethod,iterationNbr,processType)
                    return spList,lifestageList,contMethod,test,pValMethod,iterationNbr,SPS,processType
            if len(lifestageList) == 0:
                answer2=""
                while (answer2 != "y" or answer2 !="n"):
                    print("\nChoosing lifestages is a mandatory step for statistical and manual identity modules.")
                    print("Are you sure you want to skip this step ? (y / n)")
                    answer2=input()
                    if answer2=="y":
                        lifestageList = [] 
                    elif answer2=="n":    
                        sp1Lifestage,sp2Lifestage = chooseParametersLifestage(subjectList)
                        lifestageList = [sp1Lifestage,sp2Lifestage]
                    print(spList,lifestageList,test,pValMethod,iterationNbr,processType)
                    return spList,lifestageList,contMethod,test,pValMethod,iterationNbr,SPS,processType
            print(spList,lifestageList,contMethod,test,pValMethod,iterationNbr,processType)
            return spList,lifestageList,contMethod,test,pValMethod,iterationNbr,SPS,processType


def automaticSpeciesLifestage(subjectList):
    """
    Description
    -----------
    Retrieves from the input folder the txt file used for automatically use data of subjects.

    Parameters
    ----------
    None

    Returns
    ----------
    speciesList
        list, contains species of subjects for comparison.
    lifestageList
        list, contains lifestages of subjects for comparison. 
    """ 
    try:
        
        with open('input/subjectList/'+subjectList) as f:
            speciesList = f.readline()
            speciesList = speciesList.replace('\n','')
            speciesList=speciesList.split(',')
            lifestageList = f.readline()
            lifestageList = lifestageList.replace('\n','')
            lifestageList = lifestageList.split(',')
    except FileNotFoundError:
            print("\nFile "+subjectList+" is not available.\nPlease choose a new one.")
            subjectList = selectSubjectList()
            with open('input/subjectList/'+subjectList) as f: 
                speciesList = f.readline()
                speciesList = speciesList.replace('\n','')
                speciesList=speciesList.split(',')
                lifestageList = f.readline()
                lifestageList = lifestageList.replace('\n','')
                lifestageList = lifestageList.split(',')
            subjectList = subjectList
    return speciesList,lifestageList


def chooseParametersSpecies(subjectList):
    """
    Description
    -----------
    Allows the user to manually choose which species will be used for comparison.

    Parameters
    ----------
    speciesList
        list, contains species of subjects for comparison. 

    Returns
    ----------
    choiceSp1
        str, manually choose species 1 for comparison.
    choiceSp2 
        str, manually choose species 2 for comparison.
    """
    speciesList,lifestageList = automaticSpeciesLifestage(subjectList) 
    print("\nSpecies list : "+' - '.join(speciesList))
    lock = False
    lock2 = False
    while lock==False:
        print("\nSpecies 1: ") 
        choiceSp1 = input()
        if choiceSp1 in speciesList:
            lock=True
        else:
            print("\nPlease choose a species in the list")
    while lock2==False:
        print("\nSpecies 2: ") 
        choiceSp2 = input()
        if choiceSp2 in speciesList:
            lock2=True
        else:
            print("\nPlease choose a species in the list") 
    print("\nSpecies selected : "+choiceSp1+" & "+choiceSp2)
    return choiceSp1,choiceSp2


def chooseParametersLifestage(subjectList):
    """
    Description
    -----------
    Allows the user to manually choose which lifestage of the species will be used for comparison.

    Parameters
    ----------
    lifestageList
        list, contains species of subjects for comparison. 

    Returns
    ----------
    choiceSp1
        str, manually choose lifestage of species 1 for comparison.
    choiceSp2 
        str, manually choose lifestage of species 2 for comparison.
    """
    speciesList,lifestageList = automaticSpeciesLifestage(subjectList)  
    print("\nLifestage list : "+' - '.join(lifestageList))
    lock = False
    lock2 = False
    while lock==False:
        print("\nLifestage of species 1: ") 
        choiceSp1 = input()
        if choiceSp1 in lifestageList:
            lock=True
        else:
            print("\nPlease choose a lifestage in the list")
    while lock2==False:
        print("\nLifestage of species 2: ") 
        choiceSp2 = input()
        if choiceSp2 in lifestageList:
            lock2=True
        else:
            print("\nPlease choose a lifestage in the list") 
    print("\nLifestage selected : "+choiceSp1+" & "+choiceSp2)
    return choiceSp1,choiceSp2


def chooseNewSpeciesLifestage(subjectList):
    """
    Description
    -----------
    Allows the user to manually choose which lifestage and species will be used for further comparison.

    Parameters
    ----------
    None

    Returns
    ----------
    choiceSp
        str, manually choose species for comparison.
    choiceSLifestage 
        str, manually choose lifestage of species for comparison.
    """
    speciesList,lifestageList = automaticSpeciesLifestage(subjectList) 
    lock = False
    lock2 = False
    print("Species list : "+', '.join(speciesList))
    while lock==False:
        print("\nSpecies: ") 
        choiceSp1 = input()
        if choiceSp1 in speciesList:
            lock=True
        else:
            print("\nPlease choose a species in the list")
    print("Lifestage list : "+', '.join(lifestageList))
    while lock2==False:
        print("\nLifestage: ") 
        choiceSp2 = input()
        if choiceSp2 in lifestageList:
            lock2=True
        else:
            print("\nPlease choose a lifestage in the list") 
    print("\nLifestage selected : "+choiceSp1+" & "+choiceSp2)
    return choiceSp1,choiceSp2


def chooseContingencyMethod():
    """
    Description
    -----------
    Allows the user to manually choose which method will be used for building contingency matrix.

    Parameters
    ----------
    None

    Returns
    ----------
    method
        str, method for building contingency matrix.
    """ 
    methodList = ["orthopairsBased","genomeBased"]
    print("\nMethod list : "+' - '.join(methodList))
    lock = False
    while lock==False:
        print("\nContingency table building method :") 
        method = input()
        if method in methodList:
            lock=True
        else:
            print("\nPlease choose a method in the list")
    print("\nMethod selected : "+method)
    return method 


def chooseParametersTest():
    """
    Description
    -----------
    Allows the user to manually choose which statistical test will be used for testing equivalence between cell clusters.

    Parameters
    ----------
    None
    
    Returns
    -------
    test
        str, chosen statistical test.
    """

    testList=['Fisher','Chi2','Barnard','Boschloo']
    print("\nTest list : "+' - '.join(testList))
    lock = False
    while lock==False:
        print("\nStatistical test :") 
        test = input()
        if test in testList:
            lock=True
        else:
            print("\nPlease choose a statistical test in the list")
    print("\nTest selected : "+test)
    return test


def choosePValMethod():
    """
    Description
    -----------
    Allows the user to manually choose which statistical method for combining p-values from the statistical results of two different species with genomeBased method.

    Parameters
    ----------
    None
    
    Returns
    -------
    pValMethod
        str, chosen method for combining p-values in genomeBased method. 
    """
    pValMethodList=['tippett','fisher','pearson','stouffer','mudholkar_george']
    print("\nMethod list : "+' - '.join(pValMethodList))
    lock = False
    while lock==False:
        print("\nStatistical method :") 
        pValMethod = input()
        if pValMethod in pValMethodList:
            lock=True
        else:
            print("\nPlease choose a statistical method in the list")
    print("\nMethod selected : "+pValMethod)
    return pValMethod


def chooseParametersIteration():
    """
    Description
    -----------
    Let the user chooses the number of iterations he wants for bootstrap method.

    Parameters
    ----------
    None
    
    Returns
    -------
    iterationNbr
        int, chosen number of iterations.
    """

    print("\nIteration settings for bootstrap and iterative methods :")
    try:
        iterationNbr = int(input())
    except:
        print("Please enter an integer number")
    return iterationNbr


def processTypeId():
    """
    Description
    -----------
    Let the user chooses the type of processing for functional module.

    Parameters
    ----------
    None
    
    Returns
    -------
    processType
        str, chosen type of processing.
    """

    lock = False
    print("Do you want to run the identity orthologs module manually, semi_automatically or automatically ? (m / sa / a)")
    while lock==False:
        processType = input()
        if (processType == "m" or processType == "sa" or processType == "a"):
            lock=True
        else:
            print("\nPlease choose between manual, semi-automatic and automatic (m / sa / a)")
    return processType 


def selectStartingPointSubject():
    filesList = []
    filesList2 = []
    for i in os.listdir("input/SPS/"):
        if i.endswith(".txt"):
            filesList2.append(i)
            filesList.append(i.split(".")[0])
    print("\nSelect which file you want as starting point subject :\n")
    for i in range(len(filesList)):
        filesList[i] = str(i+1) + " - " + filesList[i]
    print(' | '.join(filesList))
    filesList.append("dummy")
    lock = False
    while lock == False:
        choice = int(input())
        if (choice-1) in list(range(len(filesList))) and choice != 0:
            try:
                SPS = filesList2[choice-1]
                lock = True
            except IndexError:
                print("You have to choose a number in the list, or select 0 to return to menu")
        elif choice == 0:
            return
        else:
            print("You have to choose a number in the list, or select 0 to return to menu")
    print("\nStarting point subject file chosen : "+SPS)
    return SPS


def selectSubjectList():
    filesList = []
    filesList2 = []
    for i in os.listdir("input/subjectList/"):
        if i.endswith(".txt"):
            filesList2.append(i)
            filesList.append(i.split(".")[0])
    print("\nSelect which file you want as subjects list :\n")
    for i in range(len(filesList)):
        filesList[i] = str(i+1) + " - " + filesList[i]
    print(' | '.join(filesList))
    filesList.append("dummy")
    lock = False
    while lock == False:
        choice = int(input())
        if (choice-1) in list(range(len(filesList))) and choice != 0:
            try:
                subjectList = filesList2[choice-1]
                lock = True
            except IndexError:
                print("You have to choose a number in the list, or select 0 to return to menu")
        elif choice == 0:
            return
        else:
            print("You have to choose a number in the list, or select 0 to return to menu")
    print("\nSubjects list file chosen : "+subjectList)
    return subjectList
    

def menuDisplay():
    """
    Description
    -----------
    Displays the main menu.

    Parameters
    ----------
    None
    
    Returns
    -------
    None
    """ 

    print("\n")
    print("----------------------------------------------------")
    print("|                                                  |")
    print("|                RBH_x_markers                     |")
    print("|                                                  |")
    print("|                                                  |") 
    print("|       Parameters (species - lifestage) : 1       |")
    print("|                                                  |")
    print("|       Statistics module : 2                      |")
    print("|                                                  |")
    print("|       Identity module : 3                        |")
    print("|                                                  |")
    print("|       Exit : 4                                   |")
    print("|                                                  |")
    print("----------------------------------------------------")
    print("\n")
    return


def guideline():
    """
    Description
    -----------
    Displays the guideline of the program.

    Parameters
    ----------
    None
    
    Returns
    -------
    None
    """ 

    print("\nThis program allows the user to compare the number of orthologs shared between cell clusters of two subjects.")
    print("Each subject is composed of a species and a lifestage. Choosing a subject is mandatory (parameters menu - option 1 & 2).")
    print("If clusters are already annotated, it is possible to target them directly (main menu - option 3 - not working yet).") 
    print("If not, all clusters must be selected (main menu - option 2).") 
    print("Two main types of statistical test are available, inferential and bootstrap.")
    print("The type of inferential test (Fisher's exact test, Squared Chi2, Barnard and Boschloo) can be chosen in the parameters (parameters menu - option 3 - Fisher by default).")
    print("The number of iteration for bootstraping test can also be chosen (parameters menu - option 4 - 999 iterations by default).\n")
    print("For the inferential test, a mandatory list of matrices for the test has to be build first (comparison menu - option 1)")
    print("Because this step can take quite long, the list is then saved as a CSV format.") 
    print("This list of matrices is used by the inferential test, which gives as an output a matrix in CSV format (comparison menu - option 2)")
    print("The list of matrices is also used for the bootstrap test (comparison menu - option 3).")
    print("Depending on the number of iterations, this test can take a consequent amount of time.")


# Part 2 : Global markers comparison 

# Inferential matrices preparation :

def menuStatisticsDisplay():
    """
    Description
    -----------
    Displays the statistical module menu.

    Parameters
    ----------
    None
    
    Returns
    -------
    None
    """ 

    print("\n")
    print("----------------------------------------------------")
    print("|                                                  |")
    print("|            Statistics module menu                |")
    print("|                                                  |")
    print("|                                                  |")
    print("|             infTestMatrices : 1                  |")
    print("|                                                  |")
    print("|             inferentialTest : 2                  |")
    print("|                                                  |")
    print("|             Bootstrap test : 3                   |")
    print("|                                                  |")
    print("|             Return : 4                           |") 
    print("|                                                  |")
    print("----------------------------------------------------")
    print("\n")
    return 


def menuStatistics(RBHgenesTupleList,markersDictSpList,speciesList,lifestageList,contMethod,test,pValMethod,iterationNbr,orthoMode,geneSource):
    """
    Description
    -----------
    Allows user to choose action in statistical module.

    Parameters
    ----------
    spList
        list, empty at first call, designed to contain species of subjects for comparison.
    lifestageList
        list, empty at first call, designed to contain lifestages of subjects for comparison. 
    iterationNbr
        int, number of iteration for bootstrap test.
    processType
        str, defines the processing type (manual, semi-automatic or automatic).

    Returns
    ----------
    None
    """
    while True: 
        menuStatisticsDisplay()
        while True:
            try:
                answer = int(input())
            except ValueError:
                print("\nYou must indicate an integer value ranging from 1 to 4\n")
                continue
            break
        if answer == 1:
            stat.contingencyMatrices(RBHgenesTupleList,markersDictSpList,speciesList,lifestageList,contMethod,orthoMode,geneSource)  
        elif answer == 2:
            stat.inferentialTest(RBHgenesTupleList,speciesList,lifestageList,contMethod,test,pValMethod,orthoMode,geneSource)
        elif answer == 3:
            stat.bootstrapTest(markersDictSpList,RBHgenesTupleList,speciesList,lifestageList,iterationNbr,orthoMode)
        elif answer == 4:
            return    
        

def menuIdentityDisplay():
    """
    Description
    -----------
    Displays the identity module menu.

    Parameters
    ----------
    None
    
    Returns
    -------
    None
    """ 

    print("\n")
    print("----------------------------------------------------")
    print("|                                                  |")
    print("|              Identity module menu                |")
    print("|                                                  |")
    print("|                                                  |")
    print("|             Orthologs pairs : 1                  |")
    print("|                                                  |") 
    print("|             Trim orthologs matrices : 2          |")
    print("|                                                  |")
    print("|             Return : 3                           |") 
    print("|                                                  |")
    print("----------------------------------------------------")
    print("\n")
    return

def menuIdentity(RBHgenesTupleList,markersDictSpList,speciesList,lifestageList,contMethod,test,pValMethod,processType,iterationNbr,SPS,outputName,subjectList,geneSource,orthoMode):
    """
    Description
    -----------
    Allows user to choose action in identity module.

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
    
    while True: 
        menuIdentityDisplay()
        while True:
            try:
                answer = int(input())
            except ValueError:
                print("\nYou must indicate an integer value ranging from 1 to 3\n")
                continue
            break
        if answer == 1:
            id.getOrthologsPairs(RBHgenesTupleList,markersDictSpList,speciesList,lifestageList,geneSource,orthoMode,contMethod)
        elif answer == 2:
            if processType=="m":
                id.trimOrthologsMatricesManual(speciesList,lifestageList,processType,outputName,subjectList,geneSource,orthoMode,contMethod)
            elif processType=="sa":
                id.trimOrthologsMatricesSemiAutomatic(contMethod,test,pValMethod,processType,SPS,outputName,subjectList,geneSource,orthoMode)
            elif processType=="a":
                id.fullyAutomaticProcess(contMethod,test,pValMethod,processType,iterationNbr,SPS,outputName,subjectList,geneSource,orthoMode)
        elif answer == 3:
            return   
        


