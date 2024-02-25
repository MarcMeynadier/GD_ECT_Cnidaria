#!/usr/bin/env python

###########################################
# Single-cell module
###########################################

import os
import regex as re
import scanpy as sc
from matplotlib import pyplot as plt
import PyPDF2
import glob

from PhylogeneticModule import treeGenerator


def scExpr(output,domainOut,dictTaxo,paralogsFamilies,alignment):
    if not os.path.exists("output/"+output+"/"+alignment+"/"+domainOut+"/scanpy"):
        os.makedirs("output/"+output+"/"+alignment+"/"+domainOut+"/scanpy") 
    with open("input/speciesLifestage.txt","r") as f:
        spList = f.readline()
        spList = spList.replace("\n","")
        spList = spList.split(",")
        lfList = f.readline()
        lfList = lfList.replace("\n","")
        lfList = lfList.split(",")
    f.close()
    pathSc = "input/scanpy/"
    countSubtrees = 0
    for i in paralogsFamilies:
        if os.path.exists("output/"+output+"/"+alignment+"/"+domainOut+"/scanpy/subtree"+str(countSubtrees)) == False:
             os.makedirs("output/"+output+"/"+alignment+"/"+domainOut+"/scanpy/subtree"+str(countSubtrees)) 
        pathOutput = "output/"+output+"/"+alignment+"/"+domainOut+"/scanpy/subtree"+str(countSubtrees)
        tree = "output/"+output+"/"+alignment+"/"+domainOut+"/MSA/Tree.tree"
        paralogsList(pathOutput,i)
        createSubtree(tree,i,pathOutput)
        treeGenerator(pathOutput,"subTree")
        for j in i:
            abbrv = j.split("|")[0]
            gene = j.split("|")[1]
            species = dictTaxo[abbrv][-1]
            for k in lfList:
                if os.path.exists(pathSc+species+"_"+k+".h5ad"):
                    scanpyFile = sc.read_h5ad(pathSc+species+"_"+k+".h5ad")
                    with plt.rc_context():  
                        try:
                            sc.pl.umap(scanpyFile,color=gene,title=j+" "+k,show=False)
                            plt.savefig("output/"+output+"/"+alignment+"/"+domainOut+"/scanpy/subtree"+str(countSubtrees)+"/"+j+"_"+k+"_UMAP.pdf", bbox_inches="tight")
                            plt.close()
                        except KeyError:
                            print(gene+" expression was not retrieved\n")
                            pass
        stitchUMAPs(pathOutput)
        countSubtrees += 1
    return


def stitchUMAPs(path):
    pdfFiles = []
    for file in glob.glob(path+"/*UMAP.pdf"):
        pdfFiles.append(file)
    if len(pdfFiles) > 1:
        merger = PyPDF2.PdfMerger()
        for i in pdfFiles: 
            merger.append(i)
        with open(path+"/allUMAPs.pdf",'wb') as f:
            merger.write(f)
            return
    else:
        return 

def createSubtree(tree,paralogsFamily,pathOutput):
    with open(tree,"r") as f:
        tree = f.readline()
    f.close()
    treeBeginning = paralogsFamily[0]
    treeEnding = paralogsFamily[-1]
    treeEnding = treeEnding.replace("|", "\\|")
    subTree = ""
    regex_pattern = r'\((.*?)' + re.escape(treeBeginning)
    subTree = re.split(regex_pattern, tree)
    if len(subTree) > 1:
        subTree =  subTree[-1].strip()
    match = re.search(treeEnding, subTree)
    if match:
        substring_after_regex = subTree[match.end():]
        first_comma_position = substring_after_regex.find(',')
        first_semicolon_position = substring_after_regex.find(';')
        first_separator_position = min(first_comma_position, first_semicolon_position)
        if first_separator_position != -1:
            text_before_separator = subTree[:match.end() + first_separator_position]
            subTree = text_before_separator.strip()
    subTree = treeBeginning + subTree
    subTree = subTree + ')'
    subTreeCount1 = subTree.count('(')
    subTreeCount2 = subTree.count(')')
    diff = max(subTreeCount2,subTreeCount1) - min(subTreeCount2,subTreeCount1)
    if subTreeCount1 > subTreeCount2:
        subTree = subTree + ')' * diff
    else:
        subTree = '(' * diff + subTree
    subTree = subTree + ';'
    with open(pathOutput+"/subTree.tree","w") as f:
        f.write(subTree)
    f.close()
    return


def paralogsList(path,paralogsList):
    with open(path+"/paralogs.txt","w") as f:
        f.write(','.join(paralogsList))
    f.close()
    return