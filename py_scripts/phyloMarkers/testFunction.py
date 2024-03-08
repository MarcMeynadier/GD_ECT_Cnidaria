import argparse
import subprocess,os
import numpy as np
import regex as re
import pandas as pd
import os.path
import ast
from ete3 import Tree

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


def newick_to_nested_list(newick_string):
    tree = []
    stack = []
    current_node = tree
    gene_name = ''

    for char in newick_string:
        if char == '(':
            new_node = []
            current_node.append(new_node)
            stack.append(current_node)
            current_node = new_node
        elif char == ',':
            if gene_name:
                current_node.append(gene_name)
                gene_name = ''
            current_node = stack[-1]
        elif char == ')':
            if gene_name:
                current_node.append(gene_name)
                gene_name = ''
            stack.pop()
            if stack:
                current_node = stack[-1]
        elif char == ';':
            if gene_name:
                current_node.append(gene_name)
                gene_name = ''
            break
        else:
            gene_name += char
    return tree


def is_liste_imbriquee(liste):
    for element in liste:
        if isinstance(element, list):
            return True
    return False


def verify_codes(liste, taxoCode):
    codeList = []
    for i in liste:
        code = i.split("|")[0]
        codeList.append(code)
    return set(codeList).issubset(taxoCode)

def remove_duplicates_nested(input_list):
    # Vérifier si l'élément est une liste
    if isinstance(input_list, list):
        # Utiliser une liste temporaire pour stocker les éléments uniques
        temp_list = []
        for item in input_list:
            # Appeler récursivement la fonction pour les sous-listes
            item = remove_duplicates_nested(item)
            # Vérifier si l'élément n'est pas déjà dans la liste temporaire
            if item not in temp_list:
                temp_list.append(item)
        return temp_list
    else:
        return input_list


def iterate_nested_lists(nestedList,taxoCode,fusionList,fusionNestedList):
    
    for item in nestedList:
        nestedVerif = is_liste_imbriquee(item)
        if nestedVerif == True:
            for item in nestedList:
                if isinstance(item, list):
                    #print(item)
                    fusionList,fusionNestedList = iterate_nested_lists(item,taxoCode,fusionList,fusionNestedList)
        else:
            if isinstance(item, str):
                item = item.split()
            taxoVerif = verify_codes(item,taxoCode)
            
            if taxoVerif == True:
                fusionList += item
            else:
                if len(fusionList) != 0 :
                    verifIntersect1 = set(fusionList)
                    if len(fusionNestedList) == 0:
                        fusionNestedList.append(fusionList)
                    verifIntersect2 = set(fusionNestedList[-1])
                    #print(verifIntersect1) 
                    if verifIntersect1 & verifIntersect2:
                        #print(fusionList)
                        fusionNestedList[-1] = fusionList
                    else:
                        if fusionList in fusionNestedList:
                            pass
                        else:
                            fusionNestedList.append(fusionList)
                #print(fusionList)
                fusionList = []
    #print(fusionNestedList)
    return fusionList,fusionNestedList


def fusionner_listes(nestedList,dictionnaire,taxaLimit):
    fusionList = []
    fusionNestedList = []
    # Créer un ensemble de tous les codes associés à "Cnidaria"
    taxoCode = set()
    for code, categories in dictionnaire.items():
        if taxaLimit in categories:
            taxoCode.add(code)
    for item in nestedList:
            #print(item)
            if isinstance(item, list):
                fusionList,fusionNestedList = iterate_nested_lists(item,taxoCode,fusionList,fusionNestedList)
            else:
                if isinstance(item, str):
                    item = item.split()
                taxoVerif = verify_codes(item,taxoCode)      
                if taxoVerif == True:
                    fusionList += item
                else:
                    if len(fusionList) != 0 :
                        verifIntersect1 = set(fusionList)
                        if len(fusionNestedList) == 0:
                            fusionNestedList.append(fusionList)
                        verifIntersect2 = set(fusionNestedList[-1])
                        if verifIntersect1 & verifIntersect2:
                            fusionNestedList[-1] = fusionList
                        else:
                            if fusionList in fusionNestedList:
                                pass
                            else:
                                fusionNestedList.append(fusionList)
                    fusionList = []           
    taxoVerif = verify_codes(fusionList,taxoCode)
    if taxoVerif == True:
        fusionNestedList.append(fusionList)
    newFusionNestedList = remove_duplicates_nested(fusionNestedList)
    return newFusionNestedList


def retrieveFamilies(tree,dictTaxo,taxaLimit):
    with open(tree) as f:
        chaine = f.readline()
    f.close()
    regex_pattern1 = r':(.*?)(?=\)|,)'
    regex_pattern2 = r'\)(.*?)(?=\)|,)'
    chaine = re.sub(regex_pattern1, '', chaine)
    chaine = re.sub(regex_pattern2, ")", chaine)
    nested_list = newick_to_nested_list(chaine)
    nested_list = flatten_single_element_lists(nested_list)
    paralogsFamilies = fusionner_listes(nested_list,dictTaxo,taxaLimit)
    return paralogsFamilies


def flatten_single_element_lists(nested_list):
    for i in range(len(nested_list)):
        if isinstance(nested_list[i], list) and len(nested_list[i]) == 1:
            nested_list[i] = nested_list[i][0]
        elif isinstance(nested_list[i], list):
            nested_list[i] = flatten_single_element_lists(nested_list[i])
    return nested_list










def nombre_de_listes_imbriquees(lst):
    count = 0
    for elem in lst:
        if isinstance(elem, list):
            count += 1  # Incrémente le compteur si l'élément est une liste
            count += nombre_de_listes_imbriquees(elem)  # Appel récursif pour les listes imbriquées
    return count


def retrieveFamilies2(tree,dictTaxo,taxaLimit):
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
    iterativeFlag = True
    while iterativeFlag == True:
        paralogsNestList,iterativeFlag = testMonophyly(tree,paralogsNestList,taxoCode)
    return paralogsNestList


def testMonophyly(tree,paralogsNestList,taxoCode):
    t = Tree(tree)
    iterativeFLag = False
    for i in paralogsNestList:
        listIndex = paralogsNestList.index(i)
        for j in range(len(i) - 1):
            for k in range(j + 1, len(i)):
                commonAncestor = t.get_common_ancestor(i[j],i[k])
                caLeaves = commonAncestor.get_leaves()
                caLeaves = [str(element) for element in caLeaves]
                caLeaves = [element.split('|')[0].strip('\n-') for element in caLeaves]
                if set(caLeaves) - taxoCode:
                    iterativeFLag = True 
                    cutIndex = i.index(i[k-1])
                    cutList1 = i[:cutIndex+1]
                    cutList2 = i[cutIndex+1:]                    
                    if cutList1 in paralogsNestList and cutList2 in paralogsNestList:
                        break
                    else:
                        print(commonAncestor)
                        print(cutList1)
                        print(cutList2)
                        paralogsNestList[listIndex] = cutList1 
                        paralogsNestList.insert(listIndex + 1,cutList2) 
                    break
    return paralogsNestList,iterativeFLag            
                


def createSubtree(tree,paralogsNestList,output):
    with open(tree,"r") as f:
        tree = f.readline()
    f.close()
    t = Tree(tree)
    for i in paralogsNestList:
        if len(i) > 1:
            ca = t.get_common_ancestor(i[0],i[-1])
            #ca.write(format=1, outfile=output+"new_tree.nw")
    return

tree = "output/MultidomainTest/HMM/ANF_receptor:7tm_3/MSA/Tree.tree"
t = Tree(tree)
dictTaxo = getDictTaxoFile()
paralogsNestList = retrieveFamilies2(tree,dictTaxo,"Cnidaria")
print("\n")
print(paralogsNestList)
#print(t.get_common_ancestor('Xe|LOC124452655','ast|STRG.40245'))
#createSubtree(tree,paralogsNestList)




def testMonophylyBackup(tree,paralogsNestList,taxoCode):
    t = Tree(tree)
    iterativeFLag = False
    for i in paralogsNestList:
        if len(i) > 1:
            index = paralogsNestList.index(i)
            commonAncestor = t.get_common_ancestor(i[0],i[-1])
            ancestorTree = commonAncestor.get_leaves()
            #print(ancestorTree)
            for j in ancestorTree:
                gene = str(j)
                gene = gene.split("--")[1]
                gene = gene.split("|")[0]
                if gene not in taxoCode:
                    print(commonAncestor)
                    print(gene)
                    removeGene = i[-1]
                    i.remove(removeGene)
                    paralogsNestList.insert(index + 1,removeGene.split())
                    iterativeFLag = True
                    break
    return paralogsNestList,iterativeFLag 