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
from ete3 import Tree

def processPhylogenetic(output,domainOut,msaInput,alignment):
  tree = treeBuild(output,domainOut,msaInput,alignment) # tree in Newick
  treeGenerator("output/"+output+"/"+alignment+"/"+domainOut +"/MSA","Tree") # tree in pdf
  if os.path.exists("output/metazoanTaxo.txt") == False:
    dictMetazoan = parseMetazoanList()
    createTaxoFile(dictMetazoan)
  dictTaxo = getDictTaxoFile()
  familiesList = retrieveFamilies(tree,dictTaxo,"Cnidaria")
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
    taxoCode = set()
    for code, categories in dictTaxo.items():
        if taxaLimit in categories:
            taxoCode.add(code)
    t = Tree(tree)
    paralogsNestList = []
    paralogsList = []
    for node in t.traverse("postorder"):
        gene = node.name
        if any(not c.isspace() for c in gene):
            gene = gene.split('|')[0]
            if gene in taxoCode:
                paralogsList.append(node.name)
            else:
                if len(paralogsList) != 0:
                    paralogsNestList.append(paralogsList) 
                paralogsList = []
    paralogsNestList.append(paralogsList)
    iterativeFlag1 = True
    iterativeFlag2 = True
    banLeafList = []
    banCount = 0
    while iterativeFlag1 == True:
        paralogsNestList,iterativeFlag1,banLeafList,banCount = testMonophyly(tree,paralogsNestList,taxoCode,banLeafList,banCount)
    while iterativeFlag2 == True:
      paralogsNestList,iterativeFlag2 = fusionner_listes_avec_communs(paralogsNestList)
    paralogsNestList = [subList for subList in paralogsNestList if subList]
    return paralogsNestList


def testMonophyly(tree,paralogsNestList,taxoCode,banLeafList,banCount):
    t = Tree(tree)
    iterativeFLag1 = False
    for i in paralogsNestList:
          listIndex = paralogsNestList.index(i)
          for j in range(len(i) - 1):
              for k in range(j + 1, len(i)):
                  commonAncestor = t.get_common_ancestor(i[j],i[k])
                  caLeaves = commonAncestor.get_leaves()
                  caLeaves = [str(element) for element in caLeaves]
                  caLeavesBanVerif = [element.strip('\n-') for element in caLeaves] 
                  caLeaves = [element.split('|')[0].strip('\n-') for element in caLeaves]
                  if set(caLeaves) - taxoCode:  
                    intersection = set(banLeafList) & set(caLeavesBanVerif)
                    if len(intersection) == 0:
                      iterativeFLag1 = True  
                      cutIndex = i.index(i[k-1])
                      cutList1 = i[:cutIndex+1]
                      cutList2 = i[cutIndex+1:]                    
                      if cutList1 in paralogsNestList and cutList2 in paralogsNestList:
                        banCount += 1 
                        if banCount > len(i):
                          ban = [item for item in caLeavesBanVerif if item.split('|')[0] not in taxoCode]
                          banLeafList = list(set(banLeafList) | set(ban))
                        break
                      else:
                        banCount = 0
                        paralogsNestList[listIndex] = cutList1 
                        paralogsNestList.insert(listIndex + 1,cutList2) 
                        break
    return paralogsNestList,iterativeFLag1,banLeafList,banCount


def fusionner_listes_avec_communs(liste):
    i = 0
    iterativeFlag2 = False
    while i < len(liste):
        j = i + 1
        while j < len(liste):
            if set(liste[i]) & set(liste[j]):
                iterativeFlag2 = True
                liste[i] = list(set(liste[i]) | set(liste[j]))  
                del liste[j]
                j -= 1  
            j += 1
        i += 1
    return liste,iterativeFlag2