# Webscrapping SRA sequences
# Marc Meynadier

# Packages

from bs4 import BeautifulSoup as bs
import getopt
import requests
import pandas as pd
import sys
import re 

# Functions

def argsParser():
    geo = None
    organism = None
    paper = None
    argv = sys.argv[1:]
    try:
        opts, args = getopt.getopt(argv, "g:o:p:", 
                                   ["geo =",
                                    "organism =",
                                    "paper ="])   
    except:
        print("Error")
        sys.exit(0)
    for opt, arg in opts:
        if opt in ['-g', '--geo']:
            geo = arg
        elif opt in ['-o', '--organism']:
            organism = arg
        elif opt in ['-p', '--paper']:
            paper = arg
    return geo,organism,paper 

def getGeo():
    print("\nEnter GEO accession number :\n")
    geo = input()
    return geo

def geoUrl(geo):
    url = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='
    url += geo
    return url

def fetchURL(url):
    response = requests.get(url)
    if response.status_code!=200:
        print("Error fetching page")
        sys.exit(0)
    else:
        content=response.content
    return content

def getContentGsm(content):
    soup = bs(content,'lxml')
    table = soup.find_all('table', attrs={'style':'position:relative;top:-5px;left:-5px'})
    tableStr = str(table)
    return tableStr

def getGsmTissue(contentGsm):
    gsmBs4 = re.findall(r'(">GSM)(.*?)(</a></td>)', contentGsm)
    tissueBs4 = re.findall(r'(\n<td valign="top")(.*?)(</td>)\n', contentGsm)
    tissueBs4.pop(0)
    gsmList = [] ; tissueList = []
    for i in gsmBs4:
        if i:
            i = str(i)
            gsmList.append(i)
    for i in tissueBs4:
        if i:
            i = str(i)
            tissueList.append(i)
    gsmFinalList = [] ; tissueFinalList = [] 
    for i in range(len(gsmList)):
        gsm = "GSM"
        gsmSplit = gsmList[i].split(", '",1)[1]
        gsmSplit = gsmSplit.split("'",1)[0]
        gsm += gsmSplit
        gsmFinalList.append(gsm)
        tissue = tissueList[i].split("'>",1)[1]
        tissue = tissue.split("',",1)[0]
        tissueFinalList.append(tissue)
    return gsmFinalList,tissueFinalList

def getSRX(gsmList):
    addressListSRX = []
    for i in gsmList:
        urlGsm = geoUrl(i)
        content = fetchURL(urlGsm)
        SRXcode = re.findall(r'(sra\?term=)(.*?)(</a>)',str(content)) 
        strSRX = str(SRXcode) 
        strSRX = strSRX.split("', '",1)[1]
        strSRX = strSRX.split('">',1)[0] 
        addressSRX = "https://www.ncbi.nlm.nih.gov/sra?term=" 
        addressSRX += strSRX
        addressListSRX.append(addressSRX)
    return addressListSRX

def getSRA(srxList):
    sraFinalList = []
    for i in srxList:
        content = fetchURL(i)
        sra = re.findall(r'(Traces\?run=)(.*?)(">)',str(content))
        sraList = []
        for j in sra:
            j = str(j)
            j = j.split("', '",1)[1]
            j = j.split("',",1)[0]
            sraList.append(j)
        sraFinalList.append(sraList)
    maxSra = max(len(i) for i in sraFinalList)
    for i in sraFinalList:
        while len(i) < maxSra:
            i.append('NA')
    return sraFinalList

"""
def getReadLength(sraFinalList):
    urlBeginning = 'https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc='
    urlEnd = '&display=metadata'
    readLengthList = []
    for i in sraFinalList:
        url = urlBeginning
        url += i[0]
        url += urlEnd
        response  = requests.get(url)
        data = response.json()
        data = json_normalize(data)
        print(data)
        #readLength = re.findall(r'(average length: )(.*?)(, )',str(selectDiv))
        #readLengthList.append(readLength)
    print(readLengthList)
"""

def buildOutputTable(gsmList,tissueList,sraList,geo):
    geoList = []
    infColumn = ['GSM','Sample']
    maxLenSra = 0
    for i in range(len(sraList)):
        geoList.append(geo)
        lenSra = len(sraList[i])
        if lenSra > maxLenSra:
           maxLenSra = lenSra 
    for i in range(maxLenSra):
        infColumn.append('SRA')
    dic = {k: v for k, v in zip(gsmList,sraList)}
    df = pd.DataFrame(dic)
    df.columns = geoList
    df.loc[-1] = tissueList
    df.index = df.index + 1
    df.loc[-1] = gsmList
    df.index = df.index + 1
    df.sort_index(inplace=True) 
    df.insert(0, 'GEO', infColumn)
    return df

# Main

def main():
    geo,organism,paper = argsParser()
    urlGeo = geoUrl(geo)  
    content = fetchURL(urlGeo)
    contentGsm = getContentGsm(content)
    gsmList,tissueList = getGsmTissue(contentGsm)
    addressListSRX = getSRX(gsmList)
    sraList = getSRA(addressListSRX)
    #getReadLength(sraList)
    df = buildOutputTable(gsmList,tissueList,sraList,geo)
    print(df)
    outputPath = '../../SRA_files/'
    df.to_csv(outputPath+organism+'_'+paper+'_'+geo+'_'+'SRA_info.csv',index = True,encoding='utf-8')

main()