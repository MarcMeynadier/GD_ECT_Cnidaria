import pandas as pd
import statistics
import math
import numpy as np
import pickle
from scipy.stats import shapiro 
import statsmodels.api as sm
from functools import reduce
import scipy.stats as stats
import seaborn as sns
import random
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from collections import defaultdict
import os
import sys
import re
from pheatmap import pheatmap

"""
def readPfamFile(path,file):
    with open(path+file,'r') as f:
        l = f.readline()
        while l != "":
            if l[0] == "#":
                print("test")
                l = f.readline()
                print(l)
"""
#readPfamFile('/Users/mmeynadier/Documents/PhD/species/crossSpecies/RBH_x_markers_identity/','Clytia_Nematostella_Xenia_automatic_cnidocytes.pfam')

#table = [[2, 1], [184, 316]]

a=[601, 218, 737, 247, 694, 467, 433, 203, 70, 186, 592, 552, 627, 369, 774, 570, 227, 317, 148, 215, 486, 386, 63, 169, 42, 311, 158, 220, 76, 93, 628, 41]
b=[399, 394, 263, 300, 306, 532, 452, 196, 65, 317, 408, 448, 373, 230, 226, 409, 309, 254, 249, 262, 513, 325, 186, 185, 84, 437, 157, 150, 162, 141, 372, 93]
c=[6810, 7193, 6674, 7164, 6717, 6944, 6978, 7208, 7341, 7225, 6819, 6859, 6784, 7042, 6637, 6841, 7184, 7094, 7263, 7196, 6925, 7025, 7348, 7242, 7369, 7100, 7253, 7191, 7335, 7318, 6783, 7370]
d=[18917, 18922, 19053, 19016, 19010, 18784, 18864, 19120, 19251, 18999, 18908, 18868, 18943, 19086, 19090, 18907, 19007, 19062, 19067, 19054, 18803, 18991, 19130, 19131, 19232, 18879, 19159, 19166, 19154, 19175, 18944, 19223]

pvalDict = {}

"""
for i in range(len(a)):
    table = [[a[i],b[i]],[c[i],d[i]]]
    pval = stats.fisher_exact(table,alternative='greater')
    pvalDict[i] = pval[1]
print(pvalDict)    
""" 

def combine_hex_values(d):
    """
    Description
    -----------
    

    Parameters
    ----------
    d

    Returns
    ----------
    zpad(hex(red)[2:]) + zpad(hex(green)[2:]) + zpad(hex(blue)[2:])
    """
    d_items = sorted(d.items())
    tot_weight = sum(d.values())
    red = int(sum([int(k[:2], 16)*v for k, v in d_items])/tot_weight)
    green = int(sum([int(k[2:4], 16)*v for k, v in d_items])/tot_weight)
    blue = int(sum([int(k[4:6], 16)*v for k, v in d_items])/tot_weight)
    zpad = lambda x: x if len(x)==2 else '0' + x
    return zpad(hex(red)[2:]) + zpad(hex(green)[2:]) + zpad(hex(blue)[2:])

#print(stats.fisher_exact(table))
#print(stats.fisher_exact(table,alternative='greater'))
#print(stats.fisher_exact(table,alternative='less'))

#print(stats.chi2_contingency(table,correction='True'))
#print(stats.chi2_contingency(table,correction='False'))


#print(len(nematostellaContingency[0]))

"""
print("Clytia")
for i in range(len(clytiaContingency)):
    for j in range(len(clytiaContingency[i])):
        for k in range(len(clytiaContingency[i][j])):
            if j == 16:
                if k == 9:
                    print(clytiaContingency[i][j][k]) 

print("Nematostella")
for i in range(len(nematostellaContingency)):
    for j in range(len(nematostellaContingency[i])):
        for k in range(len(nematostellaContingency[i][j])):
            if j == 0:
                if k == 0:
                    print(nematostellaContingency[i][j][k])


print(stats.fisher_exact([[105,431],[206,25705]],alternative='greater'))
print(stats.fisher_exact([[15, 498], [586, 20397]],alternative='greater'))  
"""


with open('results/stat/Clytia_adult_Clytia_larva_markers_matrix_list_1000g_genomeBased_Clytia.csv','rb') as f:
    matrixGenome = pickle.load(f)

with open('results/stat/Clytia_adult_Nematostella_adult_markers_matrix_list_1000g_orthopairsBased.csv','rb') as f:
    matrixOrtho = pickle.load(f)

print(matrixGenome)
print("\n\n\n")
#print(matrixOrtho)

print("Genome")
print(stats.fisher_exact([[68, 932], [667, 25060]],alternative='greater'))  

print("\n")
print("Ortho")
print(stats.fisher_exact([[8, 593], [45, 6765]],alternative='greater'))  

"""
listSp = ["Clytia","Nematostella"]
listLifestage=["adult","adult"]
test="Fisher"
contMethod=["genomeBased","orthopairsBased"]

genomeClytia = pd.read_csv("results/stat/"+listSp[0]+"_"+listLifestage[0]+"_"+listSp[1]+"_"+listLifestage[1]+"_"+test+"_test_1000g_"+contMethod[0]+"_tippett.csv").drop(columns=['Unnamed: 0'])
orthoClytia = pd.read_csv("results/stat/"+listSp[0]+"_"+listLifestage[0]+"_"+listSp[1]+"_"+listLifestage[1]+"_"+test+"_test_1000g_"+contMethod[1]+".csv").drop(columns=['Unnamed: 0'])

print(genomeClytia)
print("ORTHO")
print(orthoClytia)
"""
