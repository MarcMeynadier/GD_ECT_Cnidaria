import argparse
import subprocess,os
import numpy as np
import regex as re
import pandas as pd
import os.path
from Bio import SeqIO

def inputParsing(outputName):
  with open("input/"+outputName, 'r') as f:
    clusterList = f.readline() ; clusterList = clusterList.split("\n")[0]
    clusterList = clusterList.split(',')
    clusterList = map(int,clusterList)
    sp = f.readline() ; sp = sp.split("\n")[0]
    lf = f.readline() ; lf = lf.split("\n")[0]
    return clusterList,sp,lf


def getTranscriptCellType(sp,lf,clusterList):
    markersDictSp = {}
    pathMarkersSp = '../RBH_x_markers/input/Scanpy/'+sp+'_'+lf+'_markers_1000g.csv'
    dfSp = pd.read_csv(pathMarkersSp,sep='\t') 
    markersSp = dfSp.groupby('clus')['markerGene'].apply(list).values
    for i in range(len(markersSp)):
        markersDictSp[i] = markersSp[i]
    transcriptList = []
    for i in clusterList:
        transcriptList += markersDictSp[i]
    transcriptList = list(set(transcriptList))
    return transcriptList


def retrieveFastaSeq(fasta_file,targetGenes):
  protList = []
  geneList = []
  prot = ""
  with open(fasta_file, 'r') as f:
    line = f.readline()
    while line != "":
      if line[0] == ">":
        spGene = line.split(">")[1] ; spGene = spGene.split(" ")[0]
        if spGene in targetGenes or len(targetGenes) == 0:
          geneList.append(spGene) 
          line = f.readline()
          while line != "" and line[0] != ">":
            prot += line
            line = f.readline()
          protList.append(prot)
          prot = ""
        else:
          line = f.readline()
      else:
        line = f.readline()
  return geneList,protList
  

def transcriptSequences(transcriptList,prot,output):
  if os.path.exists("output/"+output+"/transcriptsFasta/flag.txt"):
    with open("output/"+output+"/transcriptsFasta/flag.txt","r") as f:
      l = f.readline()
      if "# Complete" in l:
          f.close()
          return
      else:
          f.close()
  geneList,seqList = retrieveFastaSeq(prot)
  transcriptSeq = []
  for i in transcriptList:
    for j in range(len(geneList)):
      geneList[j] = geneList[j].replace("\n","")
      if geneList[j] == ">"+i:
        transcriptSeq.append(seqList[j])
  transcriptDict = dict(zip(transcriptList,transcriptSeq))
  with open("output/"+output+"/transcriptsFasta/transcriptsFasta.fasta", 'a') as file:
      for k,v in transcriptDict.items():
          file.write(">"+k+'\n')
          file.write(v)
  file.close()
  return transcriptSeq


def transcriptCellTypeSubject(subjectInfo,output): 
  if os.path.exists("output/"+output+"/transcriptsFasta/flag.txt"):
    with open("output/"+output+"/transcriptsFasta/flag.txt","r") as f:
      l = f.readline()
      if "# Complete" in l:
          f.close()
          return
      else:
          f.close()
  for k,v in subjectInfo.items():
    sp = k.split(' ')[0]
    lf = k.split(' ')[1]
    clusterList = v.split(',')
    clusterList = map(int,clusterList)
    transcriptList = getTranscriptCellType(sp,lf,clusterList)
    if os.path.exists("input/proteins/longest_isoform_"+sp+"Proteins.fasta"):
      proteomPath = "input/proteins/longest_isoform_"+sp+"Proteins.fasta"
    else:
      proteomPath = "input/proteins/"+sp+"Proteins.fasta" 
    transcriptSequences(transcriptList,proteomPath,output)
  uniqueSequences(output,"output/"+output+"/transcriptsFasta/transcriptsFasta.fasta")
  with open("output/"+output+"/transcriptsFasta/flag.txt",'w') as f:
     f.write("# Complete")
     

def uniqueSequences(output,file):
  seen = []
  records = []
  for record in SeqIO.parse(file, "fasta"):
    if str(record.seq) not in seen:
      seen.append(str(record.seq))
      records.append(record)
  SeqIO.write(records, "output/"+output+"/transcriptsFasta/uniqueTranscriptsFasta.fasta", "fasta")
  longestSequences(output)
  


def longestSequences(output):
  with open("output/"+output+"/transcriptsFasta/uniqueTranscriptsFasta.fasta",'r') as f:
    line = f.readline()
    prot = ""
    geneList = []
    transcriptList = []
    protList = []
    while line != "":
      if line[0] == ">":
        line = line.replace("\n","")
        transcriptList.append(line)
        geneList.append(line)           
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
  with open("output/"+output+"/transcriptsFasta/uniqueTranscriptsFasta.fasta", 'w') as f:
    for k,v in gene_prot.items():
      f.write(k + '\n')
      f.write(v)
  f.close() 

