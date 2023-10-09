import sys
import os
import pandas as pd
import glob

species = sys.argv[1]
autor = sys.argv[2]
geo = sys.argv[3]
gsm = sys.argv[4]

def getBiomarkers(species,autor,geo,gsm):
    os.chdir('../../../species/'+species+'/analysis/Seurat/'+autor+'/'+geo+'/'+gsm+'/')
    search = glob.glob('*_biomarkers.csv')[0]
    biomarkers = pd.read_csv(search)
    return biomarkers

def parsedProteom(species):
    gene_prot = {}
    os.chdir('../../../species/'+species+'/raw/')
    search = glob.glob('*_proteins.fasta')[0]
    with open(search, 'r') as f:
        line = f.readline()
        prot=""
        geneList = []
        protList = []
        while line != "":
            if line[0] == ">":
                gene = line.split(' ')[1] 
                gene = gene.split('=')[1]
                gene = gene.replace('\n','')
                gene = ">" + gene
                geneList.append(gene) 
                line = f.readline()
                while line != "" and line[0] != ">":
                    prot += line
                    line = f.readline()
                protList.append(prot)
                prot = ""      
    for i in range(len(geneList)):
        gene_prot[geneList[i]] = protList[i]
    return gene_prot

def trimProteom(biomarkers,proteom):
    biomarkersList = biomarkers['gene'].tolist()
    biomarkersList2 = []
    for i in biomarkersList:
        i = ">" + i
        biomarkersList2.append(i)
    trimmedProteom = {k: proteom[k] for k in biomarkersList2 if k in proteom}
    return trimmedProteom

def saveTrimmedProteom(trimmedProteom):
    os.chdir('../../../species/'+species+'/analysis/trimFastaBiomarkers/')
    with open(species+'_'+autor+'_'+geo+'_'+gsm+'.fasta', 'w') as file:
        for k,v in trimmedProteom.items():
            file.write(k + '\n')
            file.write(v)
    file.close()

def delBlank():
    write_position = 0
    with open(species+'_'+autor+'_'+geo+'_'+gsm+'.fasta', 'rb+') as f:
        for line in f:
            if line.strip():
                read_position = f.tell()
                f.seek(write_position)
                f.write(line)
                write_position = f.tell()
                f.seek(read_position)
        f.truncate(write_position)
    f.close()
        
def main():
    initPath = os.getcwd()
    biomarkers = getBiomarkers(species,autor,geo,gsm)
    os.chdir(initPath)
    proteom = parsedProteom(species)
    trimmedProteom = trimProteom(biomarkers,proteom) 
    os.chdir(initPath)
    saveTrimmedProteom(trimmedProteom)
    delBlank()
    print("Trimmed proteom with biomarkers is ready")

main()
    

