# Folder architecture for SRA sequences
# Marc Meynadier

# Packages

import pandas as pd
import getopt
import glob
import sys
import os 

# Functions

def argsParser():
    geo = None
    argv = sys.argv[1:]
    try:
        opts, args = getopt.getopt(argv, "g", 
                                   ["geo ="])   
    except:
        print("Error")
        sys.exit(0)
    for opt, arg in opts:
        if opt in ['-g', '--geo']:
            geo = arg
    return geo

def getCSV(geo):
    script_dir = os.path.dirname(__file__)
    path = os.path.join(script_dir,'../../SRA_files')
    os.chdir(path)
    files = glob.glob('*'+geo+'*.csv')
    return files

def getGsm(files):
    gsmList = []
    for i in range(len(files)):
        df = pd.read_csv(files[i],error_bad_lines=False,sep=',')
        gsmList.append(df.iloc[0])
    realGsmList = []
    for i in gsmList:
        temporaryGsmList = [] 
        for j in i:
            if "GSM" in str(j) and str(j) != "GSM":
                temporaryGsmList.append(j)
        realGsmList.append(temporaryGsmList)
    return realGsmList   

def getMetadata(files):
    orgList = [] ; paperList = [] ; geoList = []
    for i in files:
        split1 = i.split('_',1)[0]
        orgList.append(split1)
        split2 = i.split('_',1)[1]
        split3 = split2.split('_',1)[1]
        split2 = split2.split('_',1)[0]
        split2 = split2.lower()
        paperList.append(split2)
        split3 = split3.split('_',1)[0]
        geoList.append(split3)
    return orgList,paperList,geoList

def createFolders(orgList,paperList,geoList,files,gsmList):  
    for i in range(len(orgList)):
        script_dir = os.path.dirname(__file__)
        checkPointPath = os.path.join(script_dir,'../../SRA_files')
        os.chdir(checkPointPath)
        while os.path.exists(checkPointPath+'/'+orgList[i]+'/'+paperList[i]+'/'+geoList[i]) == False:
            os.chdir(checkPointPath)
            if os.path.exists(checkPointPath+'/'+orgList[i]) == False:
                os.mkdir(orgList[i])
            else:
                path = checkPointPath
                path += '/' ; path += orgList[i]
                os.chdir(path)
                if os.path.exists(path+'/'+paperList[i]) == False:
                    os.mkdir(paperList[i])
                else:
                    path += '/' ; path += paperList[i]
                    os.chdir(path)
                    if os.path.exists(path+'/'+geoList[i]) == False:
                        os.mkdir(geoList[i])
                    else:
                        break
        os.chdir(checkPointPath+'/'+orgList[i]+'/'+paperList[i]+'/'+geoList[i])              
        createTxtSRA(files,gsmList,i,geoList)

def createTxtSRA(files,gsmList,count,geo):
    pathSRA = '../../../'
    df = pd.read_csv(pathSRA+files[count],error_bad_lines=False,sep=',')
    cols = [0,1] ; rows = [0,1]
    df = df.drop(df.columns[cols],axis=1,inplace=False)
    df = df.drop(index=rows)
    gsmCount = 0
    for j in range(len(gsmList[count])):
        if gsmCount == 0:
            SRA = df[geo[count]].values.astype(str)
        else : 
            SRA = df[geo[count]+'.'+str(gsmCount)].values.astype(str)
        SRAtxt = ''
        for k in SRA:
            if k != 'nan':
                SRAtxt += k
                SRAtxt += '\n'         
        gsmCount += 1
        fileSRA = open(gsmList[count][j]+'.txt','w')
        fileSRA.write(SRAtxt)
        fileSRA.close()
            
        
# Main

def main():
    geo = argsParser()
    files = getCSV(geo)
    gsmList = getGsm(files)
    orgList,paperList,geoList = getMetadata(files)
    createFolders(orgList,paperList,geoList,files,gsmList)

main()