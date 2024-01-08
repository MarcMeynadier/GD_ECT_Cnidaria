import argparse
import subprocess,os
import numpy as np
import regex as re
import pandas as pd
import os.path
from copy import deepcopy
from itertools import chain
from nltk import flatten
from Bio.Blast.Applications import NcbiblastpCommandline
import sys
from commonFunctions import *
from smallPicture import *
from bigPicture import *
from taxoModule import *
from scanpyLink import *



def main(args):
  if not os.path.exists("output/"+args.output):
    os.makedirs("output/"+args.output) 
  """
  if os.path.exists("output/"+args.output+"/blast/"+args.query.split(".")[0]+"BlastOut") == False:
    if os.path.exists("output/"+args.output+"/blast/") == False:
      os.makedirs("output/"+args.output+"/blast/")
    blastQuery(args.query,args.output)
  print("\nBlast of "+args.query.split(".")[0]+" has been processed.\n")
  """
  if os.path.exists("output/pfam/metazoanHmmOutput") == False:
    if os.path.exists("output/pfam/") == False:
      os.makedirs("output/pfam/")
    subprocess.call(['hmmsearch','--cut_ga','-o',"output/pfam/metazoanHmmOutput",'--domtblout',"output/pfam/metazoanDtOut",
                     "input/pfam/"+args.hmm_lib,"input/proteins/metazoanProteins.fasta"])
  print("Pfam annotation of metazoan sequences database has been processed.\n")

  domainQuery = getQueryDomains(args.query)
  if len(domainQuery) == 0:
    print("Domains were not retrieved correctly, please check if gene exist in Pfam output. Exit.")
    sys.exit()
  targetGenes = retrieveGenesPfam(domainQuery)
  if len(targetGenes) == 0:
    print("No other other genes with same domains has been retrieved. Exit.")
    sys.exit()
  geneList,protList = retrieveFastaSeq("input/proteins/metazoanProteins.fasta",targetGenes)
  if args.query != None:
    args.domain = args.query.split(".")[0]
  dtOut = "metazoanDtOut"
  # This block is for getting only domains of transcripts, running MSA and building tree
  if not os.path.exists("output/"+args.output+"/MSA"):
    os.makedirs("output/"+args.output+"/MSA")  
  if os.path.exists("output/"+args.output+"/MSA/"+args.domain+"Domains.fasta") == False:
    for i in range(len(geneList)):
      listDomain = multiDomain(geneList[i],"output/pfam/"+dtOut)
      if len(listDomain) != 0:
        flag = verifyDomains(listDomain,domainQuery)
        #flag = True
        if flag == False:
          pass
        else:
          domainDic = trimList(listDomain)
          orderList = reorderDomainDict(domainDic)
          #print(orderList)
          concatenateDomain(args.output,geneList[i],protList[i],orderList,domainDic,args.domain)
  if os.path.exists("output/"+args.output+"/MSA/"+args.domain+"MafftMSA.fasta") == False:
    mafftMSA(args.output,args.domain) 
  if os.path.exists("output/"+args.output+"/MSA/"+args.domain+"Tree.tree") == False: 
    treeBuild(args.output,args.domain)
  if os.path.exists("output/"+args.output+"/MSA/"+args.domain+"Tree.pdf") == False:
    treeGenerator(args.output,args.domain)  

  if os.path.exists("output/metazoanTaxo.txt") == False:
    dictMetazoan = parseMetazoanList()
    createTaxoFile(dictMetazoan)
  
  dictTaxo = getDictTaxoFile()
  familiesList = retrieveFamilies(args.output,args.query,dictTaxo,"Cnidaria")
  paralogsFamilies = getParalogs(familiesList) 
  #eligibleFiles = scanpyEligibleSpecies()
  scExpr(args.output,dictTaxo,paralogsFamilies,args.query.split(".")[0])



def getDictTaxoFile():
  dictTaxo = {}
  with open("output/metazoanTaxo.txt","r") as f:
    l = f.readline() 
    while l != "":
      l = l.replace("\n","")
      abbr = l.split(":")[0] 
      lineage = l.split(":")[1]
      dictTaxo[abbr] = lineage.split(",")
      l = f.readline()
  f.close()
  return dictTaxo


def retrieveFamilies(output,query,dictTaxo,taxaLimit):
  query = query.split(".")[0]
  with open("output/"+output+"/MSA/"+query+"Tree.tree") as f:
    tree = f.readline()
  f.close()
  nodesList = tree.split(',')
  spList = []
  for i in nodesList:
    sp = i.split("|")[0]
    if "(" in sp:
      sp = sp.replace("(","")
    if sp not in spList:
      spList.append(sp)
  familiesList = []
  family = []
  trapFlag = False
  for i in nodesList:
    gene = i.replace('(',"")
    gene = gene.split(":")[0]
    lineage = dictTaxo[gene.split("|")[0]]
    if taxaLimit in lineage:
      trapFlag = True
    else:
      trapFlag = False
      if len(family) != 0:
        familiesList.append(family)
        family = []
    if trapFlag == True:
      family.append(gene)
  return familiesList
      
    
def getParalogs(familiesList):
  paralogsFamilies = []
  for i in familiesList:
    familyGene = [j.split("|")[0] for j in i]
    idx = np.where(pd.DataFrame(familyGene).duplicated(keep=False))[0].tolist()
    paralogs = [i[j] for j in idx]
    if len(paralogs) != 0:
      paralogsFamilies.append(paralogs)
  return paralogsFamilies
      

def verifyDomains(listDomain,domainQuery):
  verifList = []
  for i in listDomain:
    tempList = i.split(" ")
    strippedList = list(filter(None, tempList))
    cVal = float(strippedList[11])
    if cVal < 10e-14:
      domain = strippedList[3]
      verifList.append(domain)
  flag = False
  if verifList== domainQuery:
    flag = True
  return flag


def getQueryDomains(query):
  with open("input/queries/"+query,"r") as f:
    l = f.readline()
    queryName = l.split(" ")[0]
    queryName = queryName.replace(">","")
    queryName = queryName.replace("\n","") 
  f.close()
  listDomain = []
  with open("output/pfam/metazoanDtOut","r") as f:
    l = f.readline()
    while l != "":
      if l[0] != "#":
        geneName = l.split(" ")[0]
        if geneName == queryName:
          tempList = l.split(" ")
          strippedList = list(filter(None, tempList))
          cVal = float(strippedList[11])
          if cVal < 10e-14:
            domain = strippedList[3]
            listDomain.append(domain)
      l = f.readline()
  f.close()
  return listDomain


def retrieveGenesBlast(output,query):
  with open("input/metazoanList.txt","r") as f:
    spList = f.readline()
  f.close()
  spList = spList.replace("\n","")
  spList = spList.split(",")
  targetGenesList = []
  with open("output/"+output+"/blast/"+query.split(".")[0]+"BlastOut","r") as f:
    l = f.readline()
    while l != "":
      target = l.split("\t")[1]
      targetSp = target.split("|")[0]
      if targetSp in spList:
        targetGene = target.split("\t")[0]
        targetGenesList.append(targetGene)
        l = f.readline()
      else:
        l = f.readline()
  f.close()
  return targetGenesList


def retrieveGenesPfam(domainQuery):
  dicGeneDomain = {}
  with open("output/pfam/metazoanDtOut","r") as f:
    l = f.readline()
    while l != "":
      if l[0] != "#":
        geneName = l.split(" ")[0]
        domainList = []
        try:
          outputSplit = l.split(" PF")[1]
          splitStr = ''.join(outputSplit)
          splitStr = splitStr.split(" ")
          splitStr = list(filter(None, splitStr))
          evalue = splitStr[2]
          domain = l.split("PF")[0]
          domain = domain.split(" ")
          domain = list(filter(None, domain))
          domain = domain[-1]
          if float(evalue) < 10e-14:
            try:
              dicGeneDomain[geneName].append(domain)
            except KeyError:
              dicGeneDomain[geneName] = domainList
              dicGeneDomain[geneName].append(domain)
          l = f.readline()
        except IndexError:
          l = f.readline()
      else:
        l = f.readline()
  f.close()
  targetGenes = []
  for k,v in dicGeneDomain.items():
    if v == domainQuery:
      targetGenes.append(k) 
  return targetGenes
  


def concatProteom(spList,input):
  with open('output/'+input+'/blast/'+'_'.join(spList)+'_proteins.fasta', 'w') as outfile:
    for fname in spList:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)
  outfile.close()
 

def blastQuery(query,output):
  blastp = '/Users/mmeynadier/Documents/PhD/scripts/tools/ncbi-blast/bin/blastp' 
  queryPath = 'input/queries/'+query
  dbPath = '/Users/mmeynadier/Documents/PhD/scripts/GD_ECT_Cnidaria/py_scripts/phyloMarkers/input/proteins/metazoanProteins.fasta'
  queryName = query.split('.')[0]
  outPath = "/Users/mmeynadier/Documents/PhD/scripts/GD_ECT_Cnidaria/py_scripts/phyloMarkers/output/"+output+"/blast/"+queryName+"BlastOut"
  cmd = str(NcbiblastpCommandline(cmd=blastp,query=queryPath,subject=dbPath,evalue=1e-14,out=outPath,outfmt="6 qseqid sseqid pident qcovs qlen slen length bitscore evalue"))
  os.system(cmd)


def pcoOrthogroups():
  genePCO,protPCO = retrieveFastaSeq("/Users/mmeynadier/Documents/PhD/scripts/GD_ECT_Cnidaria/py_scripts/RBH_x_markers/results/final/orthopairsCnidocytes/shared_orthologs.fasta")
  curedPCO = []
  for i in genePCO:
    i = i.replace(">","")
    i = i.replace("\n","")
    curedPCO.append(i)
  spList = retrieveSpList()
  orthoGroupsList = []
  usedCombination = []
  for i in spList:
    for j in spList:
      if (i != j and (i,j) not in usedCombination and (j,i) not in usedCombination):
        orthoGroup1 = orthogroupFromBlast(i,j,curedPCO)
        orthoGroup2 = orthogroupFromBlast(j,i,curedPCO)
        for k in orthoGroup1:
          orthoGroup2.append(k) 
        orthoGroups = mimomu(orthoGroup2)
        for k in orthoGroups:
          orthoGroupsList.append(k)
        usedCombination.append((i,j))
        usedCombination.append((j,i))
  orthoGroupsFinal = mimomu(orthoGroupsList)   
  return orthoGroupsFinal


def retrieveSpList():
  with open("input/speciesLifestage.txt","r") as f:
    sp = f.readline()
    sp = sp.replace("\n","")
    sp = sp.split(",")
  f.close()
  return sp

def pcoOrthologsTranscripts(sp,pco,output):
  pco = flatten(pco)
  transcrList = []
  transcrReturnList = []
  seqList = []
  for i in sp:
    with open("input/proteins/longest_isoform_"+i+"Proteins.fasta","r") as f:
      l = f.readline()
      while l != "":
        if l[0] == ">":
          gene = l.replace(">","") ; gene = gene.replace("\n","")
          if gene in pco:
            transcrReturnList.append(gene)
            transcrList.append(l)
            l = f.readline()
            prot = ""
            while l[0] != ">":
              prot += l
              l = f.readline()
            seqList.append(prot)
          else:
            l = f.readline()
        else:
          l = f.readline()
  f.close()
  transcriptDict = dict(zip(transcrList,seqList))
  with open("output/"+output+"/transcriptsFasta/pcoOrthologsTranscripts.fasta", 'w') as file:
      for k,v in transcriptDict.items():
        file.write(k)  
        file.write(v)
  file.close()
  return transcrReturnList,seqList


def orthogroupFromBlast(sp1,sp2,pcoList):
  threshold = float(10**-5)
  listOrthoPairs = []
  with open("/Users/mmeynadier/Documents/PhD/scripts/GD_ECT_Cnidaria/py_scripts/RBH_x_markers/input/RBH/"+sp1+'_'+sp2+'.txt',"r") as f:
    l = f.readline()
    while l != "":
      lSplit = l.split("\t")
      eVal = lSplit[-1]
      eVal = eVal.replace("\n","")
      if lSplit[0] == "qseqid":
         l = f.readline()
      else:
        if lSplit[0] in pcoList or lSplit[1] in pcoList:
          eVal = float(eVal)
          if eVal > threshold or eVal == 0.0:
            l = f.readline()
          else:
            listOrthoPairs.append([lSplit[0],lSplit[1]])  
            l = f.readline()
        else:
          l = f.readline()
  f.close()
  return listOrthoPairs


def getSubjectInfo(logPCO):
  subjectList = []
  clusterList = []
  with open(logPCO,"r") as f:
    l = f.readline()
    while l != "":
      if "Corresponding " in l:
        split = l.split("of ")[1]
        subj = split.split(" when")[0]
        cluster = split.split(" : ")[1]
        cluster = cluster.strip()
        cluster = cluster.replace(" ","")
        subjectList.append(subj)
        clusterList.append(cluster)
        l = f.readline()
      else:
        l = f.readline()
  f.close()
  subjectInfo = dict(zip(subjectList,clusterList))
  return subjectInfo


def selectOrthologs(logPCO):
  with open(logPCO,"r") as f:
    l = f.readline()
    subjectPairsList = []
    orthologsList = []
    while l != "":
      l = f.readline()
      if "Subjects compared " in l:
        split = l.split(" : ")[1]
        sp1 = split.split(" ")[0]
        if l.count(sp1) == 1:
          split = split.replace("\n","")
          subjectPairsList.append(split)
          ortho = ""
          l = f.readline()
          l = f.readline()
          while len(l.strip()) != 0 :
            ortho += l
            l = f.readline()
          ortho = ortho.replace("\n","")
          ortho = ortho.replace(" ","")
          ortho = ortho.split(",")
          orthologsList.append(ortho)
  subjectPairsOrthologs = dict(zip(subjectPairsList,orthologsList))
  return subjectPairsOrthologs 
  

def mimomu(l):
  l = deepcopy(l)
  s = set(chain.from_iterable(l))
  for i in s:
    components = [x for x in l if i in x]
    for j in components:
      l.remove(j)
    l += [list(set(chain.from_iterable(components)))]
  return l


def reconstructOrthologs(subjectPairsOrthologs):
  orthogroupsList = []
  c = 0
  initPairs = []
  for k,v in subjectPairsOrthologs.items():
    for i in v:
      initPairs.append(i)
      c += 1
      if c % 2 == 0:
        orthogroupsList.append(initPairs)
        initPairs = []
  orthogroups = mimomu(orthogroupsList)
  return orthogroups


def multiDomain(geneName,dt_out):
  listDomain = []
  if " " in geneName:
    geneCode = geneName.split(' ')[0]
  elif "\n" in geneName:
    geneCode = geneName.split('\n')[0]
  else:
    geneCode = geneName
  with open(dt_out,"r") as f:
    l = f.readline()
    while l != "":
      if l[0] != "#" and geneCode in l:
        listDomain.append(l)
        l = f.readline()
      else:
        l = f.readline()
  f.close()
  return listDomain


def trimList(listDomain):
  domainDic = {}
  for i in listDomain:
    outputSplit = i.split(" PF")[1]
    splitStr = ''.join(outputSplit)
    splitStr = splitStr.split(" ")
    pfamAcces = splitStr[0]
    pfamAcces = "PF" + pfamAcces
    splitStr = list(filter(None, splitStr))
    evalue = splitStr[2]
    begin = splitStr[15]
    end = splitStr[16]
    if float(evalue) < 10e-14:
      if pfamAcces in domainDic:
        c = 1
        pfamAcces = pfamAcces + "_" + str(c)
        while pfamAcces in domainDic:
          c += 1
          pfamAcces = pfamAcces.split("_")[0]
          pfamAcces = pfamAcces + "_" + str(c)
      domainDic[pfamAcces] = (begin,end)
  return domainDic


def reorderDomainDict(domainDic):
  orderList = []
  for k,v in domainDic.items():
    orderList.append((k,v[0]))
  orderList = sorted(orderList, key=lambda tup: int(tup[1]))
  orderReturnList = []
  for i in orderList:
    orderReturnList.append(i[0]) 
  return orderReturnList


def concatenateDomain(output,gene,prot,orderList,domainDic,domain):
  domainSeq=""
  gene = ">"+gene
  for i in orderList:
    begin = domainDic[i][0]
    end = domainDic[i][1]
    seq = prot[int(begin):int(end)]
    seq = seq.replace("\n","")
    domainSeq+=seq
  domainSeq = "\n".join(re.findall("(?s).{,64}", domainSeq))[:-1]
  
  with open("output/"+output+"/MSA/"+domain+"Domains.fasta","a") as f:
    if os.stat("output/"+output+"/MSA/"+domain+"Domains.fasta").st_size != 0:
      f.write("\n")
    f.write(gene)
    f.write("\n")
    f.write(domainSeq)
    f.write("\n")
  f.close()


def mafftMSA(output,domain):
  maftPath = "/Users/mmeynadier/Documents/PhD/scripts/tools/mafft-mac/mafft.bat"
  cmd = maftPath + " output/"+output+"/MSA/"+domain+"Domains.fasta > output/"+output+"/MSA/"+domain+"MafftMSA.fasta"
  os.system(cmd)
  

def treeBuild(output,domain):
  fasttreePath = "/Users/mmeynadier/Documents/PhD/scripts/tools/FastTree/FastTree"
  cmd = fasttreePath + " output/"+output+"/MSA/"+domain+"MafftMSA.fasta > output/"+output+"/MSA/"+domain+"Tree.tree" 
  os.system(cmd)


def treeGenerator(output,domain):
  figtreePath = "java -jar /Users/mmeynadier/Documents/PhD/scripts/tools/FigTree_v1.4.4/lib/figtree.jar"
  cmd = figtreePath + " -width 5000 -height 5000 -graphic PDF output/"+output+"/MSA/"+domain+"Tree.tree output/"+output+"/MSA/"+domain+"Tree.pdf"
  os.system(cmd)


def getPanCnidarianOrthologs(pco):
  pathPCO = "/Users/mmeynadier/Documents/PhD/scripts/GD_ECT_Cnidaria/py_scripts/RBH_x_markers/results/final/"+pco+"/shared_orthologs.fasta"
  listPCO = []
  with open(pathPCO,"r") as f:
    line = f.readline()
    while line != "":
      if line[0] == ">":
        line = line.split(">")[1]
        line = line.split("\n")[0]
        listPCO.append(line) 
        line = f.readline()
        while line != "" and line[0] != ">":
          line = f.readline()
  f.close()
  return listPCO 






if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = '')
    parser.add_argument('-p','--pco',action='store',dest='pco',help='list of potential pan-cnidarian orthologs')
    parser.add_argument('-l','--lib',action='store',dest='hmm_lib',help='Pfam-A.hmm library file',default='Pfam-A.hmm')
    parser.add_argument('-o','--output',action='store',dest='output',help='Output name',default='')
    parser.add_argument('-m','--method',action='store',dest='method',help='Method name',default='')
    parser.add_argument('-d','--domain',action='store',dest='domain',help='Domain name',default='') 
    parser.add_argument('-q','--query',action='store',dest='query',help='Protein query name',default='')  
    args = parser.parse_args()

    if (args.hmm_lib == None) and (args.output == None):
      print()
      print("--lib <list> and --output <output> are needed.")
      print()
      exit()

    print(args)
    main(args)