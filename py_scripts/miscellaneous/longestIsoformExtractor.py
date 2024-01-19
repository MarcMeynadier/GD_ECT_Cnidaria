import sys
import os
from glob import glob

def getProteom(path):
    proteom = glob(path + '*Proteins.fasta')[0]
    return proteom

def longestIsoformsNematostella(proteomPath,proteom):
    with open(proteom, 'r') as f:
        line = f.readline()
        prot=""
        geneList = []
        transcriptList = []
        protList = []
        while line != "":
            if line[0] == ">":
                transcript = line.split(' ')[0]
                transcript = transcript.split('>')[1]
                transcriptList.append(transcript)
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

    gene_transcript = {}
    gene_prot = {}

    for i in range(len(geneList)):
        gene_transcript[geneList[i]] = transcriptList[i]
        gene_prot[geneList[i]] = protList[i]

    for k1,v1 in gene_prot.items():
        maxLength = 0
        transcriptLength = len(v1)
        for k2,v2 in gene_transcript.items():
            if k1 == k2:
                if transcriptLength >= maxLength:
                    maxLength = transcriptLength
                else:
                    gene_prot.pop(k1,v1)

    fileName = proteom.split('/')[-1] 
    
    with open(proteomPath+'longest_isoform_'+fileName, 'w') as file:
        for k,v in gene_prot.items():
            file.write(k + '\n')
            file.write(v)
    file.close()

def longestIsoformsClytia(proteomPath,proteom):
    with open(proteom, 'r') as f:
        line = f.readline()
        prot=""
        geneList = []
        transcriptList = []
        protList = []
        while line != "":
            if line[0] == ">":
                transcript = line.split(' ')[1]
                transcript = transcript.replace('\n','')
                transcriptList.append(transcript)
                gene = line.split(' ')[0] 
                geneList.append(gene) 
                line = f.readline()
                while line != "" and line[0] != ">":
                    prot += line
                    line = f.readline()
                protList.append(prot)
                prot = ""

    gene_transcript = {}
    gene_prot = {}

    for i in range(len(geneList)):
        gene_transcript[geneList[i]] = transcriptList[i]
        gene_prot[geneList[i]] = protList[i]

    for k1,v1 in gene_prot.items():
        maxLength = 0
        transcriptLength = len(v1)
        for k2,v2 in gene_transcript.items():
            if k1 == k2:
                if transcriptLength >= maxLength:
                    maxLength = transcriptLength
                else:
                    gene_prot.pop(k1,v1)

    fileName = proteom.split('/')[-1] 
    
    with open(proteomPath+'longest_isoform_'+fileName, 'w') as file:
        for k,v in gene_prot.items():
            file.write(k + '\n')
            file.write(v)
    file.close()

def longestIsoformsXenia(proteomPath,proteom):
    with open(proteom, 'r') as f:
        line = f.readline()
        prot=""
        geneList = []
        transcriptList = []
        protList = []
        while line != "":
            if line[0] == ">":
                transcript = line.split(' ')[0]
                transcript = transcript.replace('>','')
                transcriptList.append(transcript)
                gene = line.split(' ')[1]
                gene= gene.replace('\n','') 
                gene = ">"+gene
                geneList.append(gene) 
                line = f.readline()
                while line != "" and line[0] != ">":
                    prot += line
                    line = f.readline()
                protList.append(prot)
                prot = ""
    gene_transcript = {}
    gene_prot = {}

    for i in range(len(geneList)):
        gene_transcript[geneList[i]] = transcriptList[i]
        gene_prot[geneList[i]] = protList[i]

    for k1,v1 in gene_prot.items():
        maxLength = 0
        transcriptLength = len(v1)
        for k2,v2 in gene_transcript.items():
            if k1 == k2:
                if transcriptLength >= maxLength:
                    maxLength = transcriptLength
                else:
                    gene_prot.pop(k1,v1)

    fileName = proteom.split('/')[-1] 
    
    with open(proteomPath+'longest_isoform_'+fileName, 'w') as file:
        for k,v in gene_prot.items():
            file.write(k + '\n')
            file.write(v)
    file.close()

def invertTranscriptGene(proteomPath,proteom):
    with open(proteom, 'r') as f:
        line = f.readline()
        prot=""
        geneList = []
        protList = []
        while line != "":
            if line[0] == ">":
                firstPart = line.split(' ')[0]
                secondPart = line.split(' ')[1]
                t = firstPart.split('>')[1]
                g = secondPart.split('=')[1]
                g = g.split('\n')[0]
                invertLine = ">"
                invertLine += g
                invertLine += " transcript="
                invertLine += t
                geneList.append(invertLine) 
                line = f.readline()
                while line != "" and line[0] != ">":
                    prot += line
                    line = f.readline()
                protList.append(prot)
                prot = ""

    gene_prot = {}

    for i in range(len(geneList)):
        gene_prot[geneList[i]] = protList[i] 

    fileName = proteom.split('/')[-1]  
    with open(proteomPath+'invert_t_g_'+fileName, 'w') as file:
        for k,v in gene_prot.items():
            file.write(k + '\n')
            file.write(v)
    file.close()


def longestIsoformsAurelia(proteomPath,proteom):
    with open(proteom, 'r') as f:
        line = f.readline()
        prot=""
        geneList = []
        transcriptList = []
        protList = []
        while line != "":
            if line[0] == ">":
                transcript = line.split(' ')[0]
                transcript = transcript.replace('>','')
                transcriptList.append(transcript)
                gene = transcript.split('.')[0]
                gene = ">"+gene
                geneList.append(gene) 
                line = f.readline()
                while line != "" and line[0] != ">":
                    prot += line
                    line = f.readline()
                protList.append(prot)
                prot = ""
    gene_transcript = {}
    gene_prot = {}

    for i in range(len(geneList)):
        gene_transcript[geneList[i]] = transcriptList[i]
        gene_prot[geneList[i]] = protList[i]

    for k1,v1 in gene_prot.items():
        maxLength = 0
        transcriptLength = len(v1)
        for k2,v2 in gene_transcript.items():
            if k1 == k2:
                if transcriptLength >= maxLength:
                    maxLength = transcriptLength
                else:
                    gene_prot.pop(k1,v1)

    fileName = proteom.split('/')[-1] 
    
    with open(proteomPath+'longest_isoform_'+fileName, 'w') as file:
        for k,v in gene_prot.items():
            file.write(k + '\n')
            file.write(v)
    file.close()


def longestIsoformHydractinia(proteomPath,proteom):
    with open(proteom,'r') as f:
        line = f.readline()
        prot=""
        geneList = []
        transcriptList = []
        protList = []
        while line != "":
            if line[0] == ">":
                line = line.replace("\n","")
                transcript = line.replace('>','')
                transcriptList.append(transcript)
                gene = transcript.split(".")[0]
                gene = ">"+gene
                geneList.append(gene)           
                line = f.readline()
                while line != "" and line[0] != ">":
                    prot += line
                    line = f.readline()
                protList.append(prot)
                prot = ""
    f.close()
    gene_transcript = {}
    gene_prot = {}
    for i in range(len(geneList)):
        gene_transcript[geneList[i]] = transcriptList[i]
        gene_prot[geneList[i]] = protList[i]
    for k1,v1 in gene_prot.items():
        maxLength = 0
        transcriptLength = len(v1)
        for k2,v2 in gene_transcript.items(): 
            if k1 == k2:
                if transcriptLength >= maxLength:
                    maxLength = transcriptLength
                else:
                    gene_prot.pop(k1,v1) 
    fileName = proteom.split('/')[-1] 
    
    with open(proteomPath+'longest_isoform_'+fileName, 'w') as file:
        for k,v in gene_prot.items():
            file.write(k + '\n')
            file.write(v)
    file.close() 


def longestIsoformGTF(species):
    path = "../../../../species/"+species+"/raw/"
    proteom = path+species+"Proteins.fasta"
    gtf = path+species+"Genes.gtf"
    gene_transcript = {}
    with open(gtf,'r') as f:
        l = f.readline()
        while l != "":
            try:
                gene = l.split("gene_id")[1]
                gene = gene.split(";")[0]
                gene = gene.replace('"','')
                gene = gene.strip()
                gene = ">"+gene
                transcript = l.split("protein_id")[1]
                transcript = transcript.split(";")[0]
                transcript = transcript.replace('"','')
                transcript = transcript.strip()
                gene_transcript[gene] = transcript
            except IndexError:
                pass
            l = f.readline()
    f.close()
    geneList = []
    with open(proteom,'r') as f:
        line = f.readline()
        prot=""
        geneList = []
        transcriptList = []
        protList = []
        while line != "":
            if line[0] == ">":
                line = line.replace("\n","")
                transcript = line.split(" ")[0]
                transcript = transcript.replace(">","")
                transcriptList.append(transcript)
                line = f.readline()
                while line != "" and line[0] != ">":
                    prot += line
                    line = f.readline()
                protList.append(prot)
                prot = ""

    f.close()
    transcript_prot = {}
    for i in range(len(transcriptList)):
       transcript_prot[transcriptList[i]] = protList[i]
    gene_prot = {}
    for k,v in gene_transcript.items():
        for k2,v2 in transcript_prot.items():
            if v == k2:
                gene_prot[k] = v2
    for k1,v1 in gene_prot.items():
        maxLength = 0
        transcriptLength = len(v1)
        for k2,v2 in gene_transcript.items(): 
            if k1 == k2:
                if transcriptLength >= maxLength:
                    maxLength = transcriptLength
                else:
                    gene_prot.pop(k1,v1) 
    fileName = proteom.split('/')[-1] 
    with open(path+'longest_isoform_'+fileName, 'w') as file:
        for k,v in gene_prot.items():
            file.write(k + '\n')
            file.write(v)
    file.close() 

def longestIsoformAiptasia(proteomPath,proteom):
    with open(proteom,'r') as f:
        line = f.readline()
        prot=""
        geneList = []
        transcriptList = []
        protList = []
        while line != "":
            if line[0] == ">":
                line = line.replace("\n","")
                transcript = line.replace('>','')
                transcript = transcript.split(" ")[0]
                transcript = ">"+transcript
                transcriptList.append(transcript)
                gene = transcript.split(".")[0]
                geneList.append(gene)           
                line = f.readline()
                while line != "" and line[0] != ">":
                    prot += line
                    line = f.readline()
                protList.append(prot)
                prot = ""
    f.close()
    gene_transcript = {}
    gene_prot = {}
    for i in range(len(geneList)):
        gene_transcript[geneList[i]] = transcriptList[i]
        gene_prot[geneList[i]] = protList[i]
    for k1,v1 in gene_prot.items():
        maxLength = 0
        transcriptLength = len(v1)
        for k2,v2 in gene_transcript.items(): 
            if k1 == k2:
                if transcriptLength >= maxLength:
                    maxLength = transcriptLength
                else:
                    gene_prot.pop(k1,v1) 
    fileName = proteom.split('/')[-1] 
    new_dict = {gene_transcript.get(key, key): value for key, value in gene_prot.items()}
    with open(proteomPath+'longest_isoform_'+fileName, 'w') as file:
        for k,v in new_dict.items():
            file.write(k + '\n')
            file.write(v)
    file.close() 




def main():
    species = sys.argv[1]
    longestIsoformGTF(species)

main()