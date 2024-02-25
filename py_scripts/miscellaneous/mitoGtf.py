import os
import sys
import fileinput
import subprocess

def getMitoProtList(species):
    rPath = os.getcwd()
    mitoProtList = []
    with open(rPath+"/../../../../species/"+species+"/raw/"+species+"MitoGenes.txt","r") as f:
        l = f.readline()
        while l != "":
            l = l.split()
            l = [x for x in l if x != '']
            if float(l[7]) < 10e-10:
                mitoProtList.append(l[0])
            l = f.readline()
    f.close()
    return list(set(mitoProtList))

def getMitoGenesList(species,mitoProtList):
    rPath = os.getcwd()
    mitoGenesList = []
    with open(rPath+"/../../../../species/"+species+"/raw/"+species+"Genes.gtf","r") as f: 
        l = f.readline()
        while l != "":
            for i in mitoProtList:
                if i in l:
                    gene = l.split('gene_id "')[1]
                    gene = gene.split('";')[0]
                    mitoGenesList.append(gene)
            l = f.readline() 
    f.close()
    return list(set(mitoGenesList))

def replaceStarOutput(species,mitoGenesList,autorDate,accessNumber):
    rPath = os.getcwd()
    rootPath = rPath+"/../../../../species/"+species+"/analysis/STARmapping/"+autorDate+"/"+accessNumber+"/"
    for dossier_actuel, sous_dossiers, fichiers in os.walk(rootPath):
        for nom_fichier in fichiers:
            chemin_fichier = os.path.join(dossier_actuel, nom_fichier)
            if nom_fichier == "features.tsv.gz":
                subprocess.run(['gunzip', chemin_fichier], check=True)
                chemin_fichier = chemin_fichier.replace(".gz","")
                for i in mitoGenesList:
                    commande = f'awk -v old_value="{i}" -v new_value="{"MITO_"+i}" \'{{gsub(old_value, new_value)}}1\' {chemin_fichier} > temp && mv temp {chemin_fichier}'
                    subprocess.run(commande,shell=True,check=True)
                subprocess.run(['gzip', chemin_fichier], check=True)

def main():
    species = sys.argv[1]
    autorDate = sys.argv[2]
    accessNumber = sys.argv[3]
    mitoProtList = getMitoProtList(species)
    mitoGenesList = getMitoGenesList(species,mitoProtList)
    replaceStarOutput(species,mitoGenesList,autorDate,accessNumber)    

main()