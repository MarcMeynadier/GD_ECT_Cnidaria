#!/usr/bin/env python

###########################################
# Phylogenetic module
###########################################

import os
import requests
import regex as re
from bs4 import BeautifulSoup
import numpy as np
import pandas as pd


def processPhylogenetic(output,domainOut,msaInput,alignment):
  tree = treeBuild(output,domainOut,msaInput,alignment) # tree in Newick
  treeGenerator("output/"+output+"/"+alignment+"/"+domainOut +"/MSA","Tree") # tree in pdf
  if os.path.exists("metazoanTaxo.txt") == False:
    dictMetazoan = parseMetazoanList()
    createTaxoFile(dictMetazoan)
  dictTaxo = getDictTaxoFile()
  familiesList = retrieveFamilies(tree,dictTaxo,"Cnidaria")
  #paralogsFamilies = getParalogs(familiesList)
  print(familiesList)
  return dictTaxo,familiesList


def treeBuild(output,domainOut,filtered_fn,alignment):
  output_file = "output/"+output+"/"+alignment+"/"+domainOut +"/MSA/Tree.tree" 
  fasttreePath = "/Users/mmeynadier/Documents/PhD/scripts/tools/FastTree/FastTree"
  cmd = fasttreePath + " " +filtered_fn+" > "+output_file 
  os.system(cmd)
  return output_file


def treeGenerator(path,fileName):
    with open(path + "/" + fileName+ ".tree","r") as f:
        l = f.readline()
    f.close()
    countGene = l.count('|')
    if countGene < 30:
        dim = 300
    elif countGene > 30:
        dim = 300 + 10 * countGene
    figtreePath = "java -jar /Users/mmeynadier/Documents/PhD/scripts/tools/FigTree_v1.4.4/lib/figtree.jar"
    cmd = figtreePath + " -width " + str(dim) + " -height " + str(dim) + " -graphic PDF "+ path + "/" + fileName+ ".tree " + path + "/" + fileName +".pdf"
    os.system(cmd)


def parseMetazoanList():
  with open("input/metazoanList.txt","r") as f:
    l = f.readline()
    abbrSpecies = ""
    wholeSpecies = ""
    while l != "":
      if l[0] == "#":
        if "'" in l:
          abbrSpecies += l
        else:
          wholeSpecies += l
        l = f.readline()
      else:
        l = f.readline()
  f.close()
  abbrSpecies = abbrSpecies.replace("#","") 
  abbrSpecies = abbrSpecies.replace("'","") 
  abbrSpecies = abbrSpecies.replace(" ","")
  abbrSpecies = abbrSpecies.replace("\n","")
  abbrSpecies = abbrSpecies.split(",")
  wholeSpecies = wholeSpecies.replace("#","")
  wholeSpecies = wholeSpecies.replace("\n","")
  wholeSpecies = wholeSpecies.split(",") 
  wholeSpecies = [x.strip() for x in wholeSpecies if x.strip()]
  dictMetazoan = {}
  for i in range(len(abbrSpecies)):
    dictMetazoan[abbrSpecies[i]] = wholeSpecies[i]
  return dictMetazoan


def createTaxoFile(dictMetazoan):
  listTaxo = []
  listSpecies = []
  for k,v in dictMetazoan.items():
    listSpecies.append(k)
    genus = v.split(" ")[0]
    species = v.split(" ")[1] 
    taxo = fetchTaxonomy(genus,species)
    listTaxo.append(taxo)
  with open("output/metazoanTaxo.txt","w") as f:
    for i in range(len(listSpecies)):
      f.write(listSpecies[i])
      f.write(":")
      f.write(",".join(listTaxo[i]))
      f.write("\n")
  f.close()


def fetchTaxonomy(genus,species):
  if species != 'sp':
    url="https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Undef&name="+genus+"+"+species+"&lvl=0&srchmode=1&keep=1&unlock"
  else:
    url="https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Undef&name="+genus+"+&lvl=0&srchmode=1&keep=1&unlock"
  request = requests.get(url)
  soup = BeautifulSoup(request.text,"lxml")
  soup_str = str(soup)
  try:
    parsedTrunk = soup_str.split('" title="superkingdom"')[1]
    parsedTrunk = parsedTrunk.split("/a></dd>")[0]
    taxo = re.findall(r'">(.*?)<',parsedTrunk)
    if genus not in taxo:
      taxo.append(genus)
    return taxo
  except IndexError:
    newUrl = fetchTaxonomyById(genus,soup_str)
    request = requests.get(newUrl)
    soup = BeautifulSoup(request.text,"lxml")
    soup_str = str(soup)
    try:
      parsedTrunk = soup_str.split('" title="superkingdom"')[1]
      parsedTrunk = parsedTrunk.split("/a></dd>")[0]
      taxo = re.findall(r'">(.*?)<',parsedTrunk)
      if genus not in taxo:
        taxo.append(genus)
      return taxo
    except IndexError:
      print("Taxonomy fetching did not work, genus might be wrong.")
      return


def fetchTaxonomyById(genus,soup_str):
  motif = r'<a(.*?)</a>'
  parseList = re.findall(motif, soup_str)
  for i in parseList:
    if '<strong>'+genus+' &lt;'+genus+'&gt;</strong>' in i:
      url = i.split('href="')[1]
      url = url.split('">')[0]
      url = url.replace("amp;","")
      url = "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/" + url
  if url:
    return url
  else:
    print("Taxonomy fetching did not work, genus might be wrong.")
    return
  

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


def retrieveFamilies(tree,dictTaxo,taxaLimit):
  with open(tree) as f:
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
    try:
      lineage = dictTaxo[gene.split("|")[0]]
    except KeyError:
      pass
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
  #print(familiesList)
  for i in familiesList:
    familyGene = [j.split("|")[0] for j in i]
    #print(familyGene)
    idx = np.where(pd.DataFrame(familyGene).duplicated(keep=False))[0].tolist()
    #print(idx)
    paralogs = [i[j] for j in idx]
    if len(paralogs) != 0:
      paralogsFamilies.append(paralogs)
  return paralogsFamilies