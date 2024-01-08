import regex as re
from domainRetrieve import *
import requests
from bs4 import BeautifulSoup





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




