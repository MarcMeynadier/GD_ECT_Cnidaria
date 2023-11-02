import os
import regex as re
import scanpy as sc
import scanpy.external as sce
from matplotlib import pyplot as plt


def scanpyEligibleSpecies(): 
    eligibleFiles = []
    for i in os.listdir("input/scanpy/"):
        if i.endswith(".h5ad"):
            eligibleFiles.append(i)
    return eligibleFiles



def scExpr(output,dictTaxo,paralogsFamilies,query):
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
        if os.path.exists("output/"+output+"/scanpy/"+query+"_subtree"+str(countSubtrees)) == False:
             os.makedirs("output/"+output+"/scanpy/"+query+"_subtree"+str(countSubtrees))
        for j in i:
            abbrv = j.split("|")[0]
            gene = j.split("|")[1]
            species = dictTaxo[abbrv][-1]
            for k in lfList:
                if os.path.exists(pathSc+species+"_"+k+".h5ad"):
                    scanpyFile = sc.read_h5ad(pathSc+species+"_"+k+".h5ad")
                    with plt.rc_context():  
                        try:
                            sc.pl.umap(scanpyFile,color=gene,show=False)
                        except KeyError:
                            print(gene+" expression was not retrieved\n")
                        plt.close()
                        plt.savefig("output/"+output+"/scanpy/"+query+"_subtree"+str(countSubtrees)+"/"+j+"_"+k+".pdf", bbox_inches="tight")
        countSubtrees += 1
                    
            
def cellTypeAssessment(output,query):
     
















