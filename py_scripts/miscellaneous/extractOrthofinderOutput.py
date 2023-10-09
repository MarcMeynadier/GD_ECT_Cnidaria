"""
Extract orthogroups, orthologs and paralogs from OrthoFinder outputs
Marc Meynadier
"""

# Packages and parameters
import sys
import os
import pandas as pd
from glob import glob
import regex as re

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', None)

# Functions module 1

def getOrthoFile(path,fastaNames):
    path += 'Orthologues/'
    fastaNameParsed = []
    for i in fastaNames:
        i = i.replace('.fasta','')
        fastaNameParsed.append(i) 
    dfList = []
    for i in fastaNameParsed:
        speciesPath = path + 'Orthologues_' + i
        df = pd.read_csv(glob(speciesPath+"/*.tsv")[0], sep="\t")
        dfList.append(df)
    return dfList

def getParaFile(path):
   path += 'Gene_Duplication_Events/' 
   df = pd.read_csv(glob(path+"*.tsv")[0], sep="\t")
   return df 

def extractGenes(col): 
    genesL = []
    for i in col:
        genesL.append(i)
    genesStr = ' '.join(genesL) 
    genesStr = genesStr.replace(',','') 
    genesList = genesStr.split() 
    return genesList

def extractOrthologs(orthofile):
    genesColumnSpecies1 = orthofile[orthofile.columns[1]]
    genesColumnSpecies2 = orthofile[orthofile.columns[2]]
    genesSpecies1 = extractGenes(genesColumnSpecies1)
    genesSpecies2 = extractGenes(genesColumnSpecies2)
    return genesSpecies1,genesSpecies2

def extractParalogs(parafile):
    print(parafile)

# Functions module 2

def common(a,b): 
    c = [value for value in a if value in b] 
    return c

def extractBiomarkersMatches(listMarkers,orthofiles,parafile):
    orthoBiomarkerMatchesAll = []
    orthoHeadersAll = []
    for i in orthofiles:
        orthoHeaders = list(i.columns)
        orthoHeadersAll.append(orthoHeaders)
        orthoBiomarkerMatches = pd.DataFrame()
        for j in listMarkers:
            index = 0
            for k in i.iterrows():
                if j in str(k[1]):
                    orthoBiomarkerMatches = pd.concat([orthoBiomarkerMatches,i.loc[index]])
                index += 1
        orthoBiomarkerMatchesAll.append(orthoBiomarkerMatches)
    paraBiomarkerMatches = pd.DataFrame()
    for i in listMarkers:
        index = 0
        for j in parafile.iterrows():
            if i in str(j[1]):
                paraBiomarkerMatches = pd.concat([paraBiomarkerMatches, parafile.loc[index]]) 
            index +=1
    return orthoHeadersAll,orthoBiomarkerMatchesAll,paraBiomarkerMatches

def metadataBiomarkersMatches(listMarkers,orthoMarkerMatchesAll,paraMarkerMatches,orthoHeadersAll):
    biomarkerList = []
    infoList = []
    potentialBiomarker = ""
    infoStr = ""
    for i in paraMarkerMatches.iterrows():
        infoStr += str(i)
        infoStr += '\n'
        if 'Genes 1' in str(i):
            for j in listMarkers:
                if j in str(i):
                    potentialBiomarker = j
        if 'Genes 2' in str(i):
            for j in listMarkers:
                if j in str(i):
                    potentialBiomarker = j  
            biomarkerList.append(potentialBiomarker)
            potentialBiomarker = ""
            infoList.append(infoStr)
            infoStr = ""
    zipList = list(zip(biomarkerList,infoList))
    paraDf = pd.DataFrame(zipList,columns=['Biomarkers','Informations'])
    orthoDfList = []
    for i in range(len(orthoMarkerMatchesAll)):
        biomarkerList = []
        infoList = []
        potentialBiomarker = ""
        infoStr = ""
        for j in orthoMarkerMatchesAll[i].iterrows():
            infoStr += str(j)
            infoStr += '\n'
            if orthoHeadersAll[i][1] in str(j):
                for k in listMarkers:
                    if k in str(j):
                        potentialBiomarker = k
            if orthoHeadersAll[i][2] in str(j):
                for k in listMarkers:
                    if k in str(j):
                        potentialBiomarker = k 
                biomarkerList.append(potentialBiomarker)
                potentialBiomarker = ""
                infoList.append(infoStr)
                infoStr = ""
        zipList = list(zip(biomarkerList,infoList))
        orthoDf =  pd.DataFrame(zipList,columns=['Biomarkers','Informations'])
        orthoDfList.append(orthoDf)
    return orthoDfList,paraDf

def saveDf(eooPath,orthoDfList,speciesList,paraDf):
    if os.path.exists(eooPath) == False:
        os.mkdir(eooPath)
    for i in range(len(orthoDfList)):
        orthoDfList[i].to_csv(eooPath+speciesList[i]+'_orthoBiomarkers.csv',index=False)
    paraDf.to_csv(eooPath+'paraBiomarkers.csv',index=False) 

# Functions module 3

def getFastaNames(path):
    with open(path+'Log.txt','r') as f:
        lines = f.readlines()
        lines = ''.join(lines)
        fastaFiles = re.findall(r'\d: (.*)\n',lines)
    f.close()
    return fastaFiles
        
def parseBiomatches(path,speciesList,fastaNames):
    orthoDict1 = {}
    orthoDict2 = {}
    paraDict = {}
    fastaNameParsed = []
    for i in fastaNames:
        i = i.replace('.fasta','')
        fastaNameParsed.append(i)
    paraName = ['Genes 1','Genes 2']
    ortho1 = pd.read_csv(path + speciesList[0]+'_orthoBiomarkers.csv')
    ortho2 = pd.read_csv(path + speciesList[1]+'_orthoBiomarkers.csv')
    para = pd.read_csv(path + 'paraBiomarkers.csv')
    for index, row in ortho1.iterrows():
        info = row['Informations']
        biomatches = []
        for j in fastaNameParsed:
            regex = r'{0}(.*)\nName:'.format(j)
            markers = re.findall(regex,info)
            markers = ''.join(markers)
            markers = markers.split(' 0 ')[1]
            markers = markers.strip()
            if 'XLOC' in markers:
                markers = markers.replace('type_','type:')
                markers = markers.replace('len_','len:')
                markers = markers.replace('_+_','(+)')
                markers = markers.replace('_-_','(-)')
                markers = re.sub(r'(TCONS_\d+)_',r'\1:',markers)
            markers = markers.split(',')
            for k in markers:
                k = '>' + k
                k = k.replace('> ','>')
                biomatches.append(k)
        biomarker = row['Biomarkers']
        orthoDict1[biomarker] = biomatches
    for index, row in ortho2.iterrows():
        info = row['Informations']
        biomatches = []
        for j in fastaNameParsed:
            regex = r'{0}(.*)\nName:'.format(j)
            markers = re.findall(regex,info)
            markers = ''.join(markers)
            markers = markers.split(' 0 ')[1]
            markers = markers.strip()
            if 'XLOC' in markers:
                markers = markers.replace('type_','type:')
                markers = markers.replace('len_','len:')
                markers = markers.replace('_+_','(+)')
                markers = markers.replace('_-_','(-)')
                markers = re.sub(r'(TCONS_\d+)_',r'\1:',markers)
            markers = markers.split(',')
            for k in markers:
                k = '>' + k
                k = k.replace('> ','>')
                biomatches.append(k)
        biomarker = row['Biomarkers']
        orthoDict2[biomarker] = biomatches
    for index, row in para.iterrows():
        info = row['Informations']
        biomatches = []
        for j in paraName:
            regex = r'{0}(.*)\nName:'.format(j)
            markers = re.findall(regex,info)
            markers = ''.join(markers)
            markers = markers.split(' 0 ')[1]
            markers = markers.strip()
            if 'XLOC' in markers:
                markers = markers.replace('type_','type:')
                markers = markers.replace('len_','len:')
                markers = markers.replace('_+_','(+)')
                markers = markers.replace('_-_','(-)')
                markers = re.sub(r'(TCONS_\d+)_',r'\1:',markers)
            markers = markers.split(',')
            for k in markers:
                k = '>' + k
                k = k.replace('> ','>')
                biomatches.append(k)
        biomarker = row['Biomarkers']
        paraDict[biomarker] = biomatches
    return orthoDict1,orthoDict2,paraDict
            
def prepareFastaBiomatches(orthoFinderPath,eooPath,speciesList,orthoDict1,orthoDict2,paraDict,fastaNames):
    speciesFastaList = []
    for i in fastaNames:
        with open(orthoFinderPath + i,'r') as f:
            fasta = f.readlines()
            fasta = ''.join(fasta)
            speciesFastaList.append(fasta)
    fastaConcatenate = '\n'.join(speciesFastaList)       
    if os.path.exists(eooPath+speciesList[0]+'_orthoBiomatches') == False:
        os.mkdir(eooPath+speciesList[0]+'_orthoBiomatches')
    if os.path.exists(eooPath+speciesList[1]+'_orthoBiomatches') == False:
        os.mkdir(eooPath+speciesList[1]+'_orthoBiomatches')
    if os.path.exists(eooPath+'paraBiomatches') == False:
        os.mkdir(eooPath+'paraBiomatches')  
    for k1,v1 in orthoDict1.items():
        fastaDict = {}
        for i in v1:
            fastaSplit = fastaConcatenate.split(i)[1]
            fastaSplit = fastaSplit.split('>',1)[0]
            fastaSplit = fastaSplit.split('\n',1)[1]
            fastaDict[i] = fastaSplit
            filename = i.split('>')[1]
            filename = filename.split(' ')[0]
            with open(eooPath+speciesList[0]+'_orthoBiomatches/'+filename+'.fasta', 'w') as file:
                for k2,v2 in fastaDict.items():
                    file.write(k2 + '\n')
                    file.write(v2)
    for k1,v1 in orthoDict2.items():
        fastaDict = {}
        for i in v1:
            fastaSplit = fastaConcatenate.split(i)[1]
            fastaSplit = fastaSplit.split('>',1)[0]
            fastaSplit = fastaSplit.split('\n',1)[1]
            fastaDict[i] = fastaSplit
            filename = i.split('>')[1]
            filename = filename.split(' ')[0]
            with open(eooPath+speciesList[1]+'_orthoBiomatches/'+filename+'.fasta', 'w') as file:
                for k2,v2 in fastaDict.items():
                    file.write(k2 + '\n')
                    file.write(v2)
    fastaNameParsed = []
    for i in fastaNames:
        i = i.replace('.fasta','')
        fastaNameParsed.append(i)   
    for k1,v1 in paraDict.items():
        fastaDict = {}
        for i in v1:
            for j in fastaNameParsed:
                if j in i:
                    i = i.replace(j,'')
                    i = i.replace('_','',1)
            fastaSplit = fastaConcatenate.split(i)[1]
            fastaSplit = fastaSplit.split('>',1)[0]
            fastaSplit = fastaSplit.split('\n',1)[1]
            fastaDict[i] = fastaSplit
            filename = i.split('>')[1]
            filename = filename.split(' ')[0]
            with open(eooPath+'paraBiomatches/'+filename+'.fasta', 'w') as file:
                for k2,v2 in fastaDict.items():
                    file.write(k2 + '\n')
                    file.write(v2) 

            
# Menu

def menuDisplay():
    print("\n")
    print("----------------------------------------------------------------")
    print("|                                                              |")
    print("|                  Extract OrthoFinder Output                  |")
    print("|                                                              |")
    print("|      Extract raw list of orthologs / paralogs : 1            |") 
    print("|                                                              |")
    print("|      Extract orthologs and paralogs from biomatches : 2      |")
    print("|                                                              |")
    print("|      Retrieve protein sequences of biomatches : 3            |")
    print("|                                                              |")
    print("|      Exit : 4                                                |")
    print("|                                                              |")
    print("----------------------------------------------------------------")
    print("\n")
    return

def mainMenu():
    with open("markerList.txt") as f:
        listMarkers = [x.rstrip() for x in f] 
    with open("speciesList.txt") as f:
        speciesList = [x.rstrip() for x in f]
    speciesStr = '_'.join(speciesList)
    path = '../../../species/crossSpecies/'+speciesStr+'/'
    orthoFinderPath = path + 'OrthoFinder/'
    eooPath = path +'EOO/'
    fastaNames = getFastaNames(orthoFinderPath)
    orthofiles = getOrthoFile(orthoFinderPath,fastaNames)
    parafile = getParaFile(orthoFinderPath) 
    while True:
        menuDisplay()
        while True:
            try:
                answer = int(input())
            except ValueError:
                print("\nYou must indicate an integer value ranging from 1 to 4\n")
                continue
            break 
        if answer == 1:
            orthoList1, orthoList2 = extractOrthologs(orthofiles)
        elif answer == 2:
            orthoHeadersAll,orthoMarkerMatchesAll,paraMarkerMatches = extractBiomarkersMatches(listMarkers,orthofiles,parafile)
            orthoDfList,paraDf = metadataBiomarkersMatches(listMarkers,orthoMarkerMatchesAll,paraMarkerMatches,orthoHeadersAll)
            saveDf(eooPath,orthoDfList,speciesList,paraDf)
        elif answer == 3: 
            orthoDict1,orthoDict2,paraDict = parseBiomatches(eooPath,speciesList,fastaNames)
            prepareFastaBiomatches(orthoFinderPath,eooPath,speciesList,orthoDict1,orthoDict2,paraDict,fastaNames)
        elif answer == 4:
            sys.exit(0)

# Program activation

mainMenu()