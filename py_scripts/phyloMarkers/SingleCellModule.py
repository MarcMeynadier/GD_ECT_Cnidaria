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
from ete3 import Tree

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
        
        tree = "output/"+output+"/"+alignment+"/"+domainOut+"/MSA/Tree.tree"
        if len(i) > 1:
            if os.path.exists("output/"+output+"/"+alignment+"/"+domainOut+"/scanpy/subtree"+str(countSubtrees)) == False:
                os.makedirs("output/"+output+"/"+alignment+"/"+domainOut+"/scanpy/subtree"+str(countSubtrees)) 
                pathOutput = "output/"+output+"/"+alignment+"/"+domainOut+"/scanpy/subtree"+str(countSubtrees)
            createSubtree(tree,i,pathOutput)
            treeGenerator(pathOutput,"subTree")
            paralogsList(pathOutput,i)
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
        else:
            isolatedPath = "output/"+output+"/"+alignment+"/"+domainOut+"/scanpy/isolatedParalogs" 
            if os.path.exists(isolatedPath) == False:
                os.makedirs(isolatedPath)
            paralogsList(isolatedPath,i)
            abbrv = i[0].split("|")[0]
            gene = i[0].split("|")[1]
            species = dictTaxo[abbrv][-1]
            for j in lfList:
                if os.path.exists(pathSc+species+"_"+j+".h5ad"):
                    scanpyFile = sc.read_h5ad(pathSc+species+"_"+j+".h5ad")
                    with plt.rc_context():  
                        try:
                            sc.pl.umap(scanpyFile,color=gene,title=i[0]+" "+j,show=False)
                            plt.savefig(isolatedPath+"/"+i[0]+"_"+j+"_UMAP.pdf", bbox_inches="tight")
                            plt.close()
                        except KeyError:
                            print(gene+" expression was not retrieved\n")
                            pass
    if os.path.exists("output/"+output+"/"+alignment+"/"+domainOut+"/scanpy/isolatedParalogs/paralogs.txt"):  
        trimLastComma("output/"+output+"/"+alignment+"/"+domainOut+"/scanpy/isolatedParalogs/paralogs.txt")    
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


def createSubtree(tree,paralogsList,output):
    with open(tree,"r") as f:
        tree = f.readline()
    f.close()
    t = Tree(tree)      
    ca = t.get_common_ancestor(paralogsList[0],paralogsList[-1])
    ca.write(format=1, outfile=output+"/subTree.tree")
    return


def paralogsList(path,paralogsList):
    with open(path+"/paralogs.txt","a") as f:
        f.write(','.join(paralogsList))
        if len(paralogsList) == 1:
            f.write(',')
    f.close()
    return

def trimLastComma(path):
    with open(path, "r") as f:
        content = f.read()
    f.close()
    with open(path, "w") as f:
        f.write(content[:-1])
    f.close()