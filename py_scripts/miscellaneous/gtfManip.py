import csv
import os
import glob

def geneIdTransform(species):
    path = "../../../../species/"+species+"/raw/"
    gtf = path+species+"Genes.gtf"
    geneProt = {}
    with open(gtf,"r") as f:
        l = f.readline()
        while l != "":
            if 'gene_id' in l and '\tCDS\t' in l:
                gene = l.split('gene_id "')[1]
                gene = gene.split('"')[0]
                prot = l.split(' ID "')[1]
                prot = prot.split('"')[0]
                prot = prot.split("-")[1]
                geneProt[gene] = prot
            l = f.readline()
    f.close()     
    return geneProt


def replaceIdSingleCell(species,geneProt):
    path = "../../../../species/"+species+"/analysis/ScanPy/"
    
    for file in os.listdir(path):
        if ".csv" in file:
            markersFile = file
    col1 = [] ; col2 = [] ; col3 = [] ; col4 = []
    with open(path+markersFile,'r',) as f:
        l = f.readline()
        while l != "":
            lineSplit = l.split("\t")
            col1.append(lineSplit[0])
            try:
                col2.append(geneProt[lineSplit[1]])
            except KeyError:
                col2.append(lineSplit[1])
            col3.append(lineSplit[2])
            col4.append(lineSplit[3])
            l = f.readline()
    f.close()
    with open(path+markersFile, 'w') as f:
        for i in range(len(col1)):
            f.write(col1[i]+"\t"+col2[i]+"\t"+col3[i]+"\t"+col4[i])
    f.close()
        

species = "Aiptasia"
geneProt = geneIdTransform(species)
replaceIdSingleCell(species,geneProt)
    

    