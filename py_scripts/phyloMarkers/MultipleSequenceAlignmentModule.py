#!/usr/bin/env python

###########################################
# Multiple Sequence Alignment Module
###########################################


import subprocess
from Bio import SeqIO
from Bio import AlignIO
from collections import defaultdict
import pysam
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os


################
# MAFFT
################


def processMAFFT(hmm_list,alignment,output):
  if len(hmm_list) == 0:
    print("Domains were not retrieved correctly, please check if gene exists in Pfam output. Exit.")
    pass
  if len(hmm_list) > 1:
    domainOut = ':'.join(hmm_list)
  else:
    domainOut = hmm_list[0] 
  targetGenes = retrieveGenesPfam(hmm_list)
  if len(targetGenes) == 0:
    print("No other other genes with same domains has been retrieved. Exit.")
    pass
  if os.path.exists("output/"+output+"/"+alignment+"/"+domainOut):
    pass
  else:
    print("Query : "+domainOut+"\n")
    os.makedirs("output/"+output+"/"+alignment+"/"+domainOut)
    geneList,protList = retrieveFastaSeq("input/proteins/metazoanProteinsDb.fasta",targetGenes)
    dtOut = "resolved_pfam_results_metazoan_domtblout.table"
    # This block is for getting only domains of transcripts, running MSA and building tree
    if not os.path.exists("output/"+output+"/"+alignment+"/"+domainOut +"/MSA"):
      os.makedirs("output/"+output+"/"+alignment+"/"+domainOut +"/MSA")  
    if os.path.exists("output/"+output+"/"+alignment+"/"+domainOut +"/MSA/reconstructedDomains.fasta") == False:
      for i in range(len(geneList)):
        listDomain = multiDomain(geneList[i],"output/pfam/"+dtOut)
        concatenateDomain(output,alignment,domainOut,listDomain,geneList[i],protList[i])
    if os.path.exists("output/"+output+"/"+domainOut+"/MSA/MafftMSA.fasta") == False:
      msaInput = mafftMSA(output,domainOut,alignment) 
  return domainOut,msaInput


def retrieveGenesPfam(hmm_list):
  dicGeneDomain = {}
  with open("output/pfam/resolved_pfam_results_metazoan_domtblout.table","r") as f:
    l = f.readline()
    while l != "":
      if l[0] != "#":
        try:
          splitLine = l.split(" ")
          geneName = splitLine[0]
          domainList = []
          domain = splitLine[1]
          eVal = float(splitLine[-1])
          if eVal < 10e-5:
            try:
              dicGeneDomain[geneName].append(domain)
            except KeyError:
              dicGeneDomain[geneName] = domainList
              dicGeneDomain[geneName].append(domain)
          l = f.readline()
        except IndexError:
          l = f.readline()
      else:
        l = f.readline()
  f.close()
  targetGenes = []
  for k,v in dicGeneDomain.items():
    if v == hmm_list:
      targetGenes.append(k) 
  return targetGenes


def retrieveFastaSeq(fasta_file,targetGenes):
  protList = []
  geneList = []
  prot = ""
  with open(fasta_file, 'r') as f:
    line = f.readline()
    while line != "":
      if line[0] == ">":
        spGene = line.split(">")[1] 
        if " " in spGene:
          spGene = spGene.split(" ")[0]
        if "\n" in spGene:
          spGene = spGene.split("\n")[0]  
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


def multiDomain(geneName,dt_out):
  listDomain = []
  if " " in geneName:
    geneCode = geneName.split(' ')[0]
  elif "\n" in geneName:
    geneCode = geneName.split('\n')[0]
  else:
    geneCode = geneName
  with open(dt_out,"r") as f:
    l = f.readline()
    l = l.replace("\n","")
    while l != "":
      if l[0] != "#" and geneCode in l:
        listDomain.append(l)
        l = f.readline()
      else:
        l = f.readline()
  f.close()
  return listDomain


def concatenateDomain(output,alignment,domainOut,listDomain,gene,prot):
  coordTupleList = []
  for i in listDomain:
    splitStr = i.split(" ")
    coord = splitStr[3]
    begin = coord.split("-")[0]
    end = coord.split("-")[1]
    coordTupleList.append((begin,end))
  domainSeq = ""
  gene = ">" + gene
  if '(' in gene:
    gene = gene.replace('(','-')
  if ')' in gene:
    gene = gene.replace(')','-')
  for i in coordTupleList:
    begin = i[0]
    end = i[1]
    seq = prot[int(begin):int(end)]
    seq = seq.replace("\n","")
    domainSeq += seq
  domainSeq = '\n'.join([domainSeq[i:i+80] for i in range(0, len(domainSeq), 80)])
  with open("output/"+output+"/"+alignment+"/"+domainOut +"/MSA/reconstructedDomains.fasta","a") as f:
    if os.stat("output/"+output+"/"+alignment+"/"+domainOut +"/MSA/reconstructedDomains.fasta").st_size != 0:
      f.write("\n")
    f.write(gene)
    f.write("\n")
    f.write(domainSeq)
  f.close()


def mafftMSA(output,domainOut,alignment):
  mafftPath = "/Users/mmeynadier/Documents/PhD/scripts/tools/mafft-mac/mafft.bat"
  msaFile = "output/"+output+"/"+alignment+"/"+domainOut +"/MSA/MafftMSA.fasta" 
  cmd = mafftPath + " output/"+output+"/"+alignment+"/"+domainOut +"/MSA/reconstructedDomains.fasta > " + msaFile
  os.system(cmd)
  return msaFile


################
# HMM
################


def processHMM(hmm_list,output,alignment,mode):
    if len(hmm_list) > 1:
      domainOut = ':'.join(hmm_list)
    else:
      domainOut = hmm_list[0] 
    db_fn = 'input/proteins/metazoanProteinsDb.fasta' 
    if not os.path.exists("output/"+output+"/"+alignment+"/"+domainOut):
        os.makedirs("output/"+output+"/"+alignment+"/"+domainOut)  
    for hmm_name in hmm_list:
      hmm_file = hmm_extract(output,domainOut,hmm_name,'input/pfam/Pfam-A.hmm') 
    target_seqs = []
    fasta = pysam.Fastafile(db_fn) 
    hitDomainDict = retrieveHitList(hmm_list)
    if mode == 'strict':
      hitDomainDict = lengthTrimHitDomainDict(hmm_list,hitDomainDict)
      hitDomainDict = exactTrimHitDomainDict(hmm_list,hitDomainDict)
    hit_list = hitDomainDict.keys()
    for hit_id in hit_list:
      seq = fasta.fetch(hit_id)
      newSeq = domainBoundaries(seq,hit_id,hitDomainDict) # This line reconstruct the sequence with only the domains
      seq_obj = Seq(newSeq)
      seq_rec = SeqRecord(seq_obj,id=hit_id)
      target_seqs.append(seq_rec)
    if not os.path.exists("output/"+output+"/"+alignment+"/"+domainOut+"/MSA"):
        os.makedirs("output/"+output+"/"+alignment+"/"+domainOut+"/MSA") 
    target_db_fn = mk_target_db(output,domainOut,hmm_name,target_seqs,alignment) # targets
    stk_fn = do_hmmalign(output,domainOut,hmm_name,hmm_file,target_db_fn,alignment) # stk
    alimask_fn = do_esl_alimask(output,domainOut,hmm_name,stk_fn,alignment) # alimask
    alimanip_fn = do_esl_alimanip(output,domainOut,hmm_name,alimask_fn,alignment) # alimanip
    mfa_fn = convert_to_mfa(output,domainOut,hmm_name,alimanip_fn,alignment) # convert to mfa
    msaInput = filter_mfa(output,domainOut,hmm_name,mfa_fn,alignment) # filter
    return domainOut,msaInput


def extractQueriesFromFasta(fasta):
  queriesList = []
  with open(fasta,'r') as f:
    l = f.readline()
    while l != "":
      if l[0] == ">":
        gene = l.replace(">","")
        gene = gene.replace("\n","")
        queriesList.append(gene)
      l = f.readline()
  f.close()
  return queriesList


def getQueryDomains(query):
  listDomain = []
  with open("input/pfam/resolved_pfam_results_metazoan_domtblout.table","r") as f:
    l = f.readline()
    while l != "":
      if l[0] != "#":
        geneName = l.split(" ")[0]
        if geneName == query:
          splitList = l.split(" ")
          cVal = float(splitList[6])
          if cVal < 10e-5:
            domain = splitList[1]
            listDomain.append(domain)
      l = f.readline()
  f.close()
  return listDomain


def hmm_extract(output,domainOut,hmm,hmmlib):
   hmm_fn = "output/"+output+"/HMM/"+domainOut +"/"+hmm + '.hmm'
   subprocess.call(['hmmfetch','-o',hmm_fn,hmmlib,hmm])
   return hmm_fn


def retrieveHitList(hmm_list):
   hitDomainDict = {}
   with open('input/pfam/resolved_pfam_results_metazoan_domtblout.table','r') as f:
      l = f.readline()
      while l != "":
         if l[0] == "#":
            l = f.readline()
         splitHit = l.split(" ")
         if splitHit[1] in hmm_list:
            hit = splitHit[0]
            domainBoundary = splitHit[4]
            try:
               hitDomainDict[hit].append(domainBoundary) 
            except KeyError:
               hitDomainDict[hit] = []
               hitDomainDict[hit].append(domainBoundary)
         l = f.readline()
   f.close()
   return hitDomainDict


def lengthTrimHitDomainDict(hmm_list,hitDomainDict):
   domainNumber = len(hmm_list)
   hitDomainDict = {key: value for key, value in hitDomainDict.items() if len(value) == domainNumber}
   return hitDomainDict


def exactTrimHitDomainDict(hmm_list,hitDomainDict):
   with open('input/pfam/resolved_pfam_results_metazoan_domtblout.table','r') as f:
      l = f.readline()
      while l != "":
         if l[0] == "#":
            l = f.readline()
         splitHit = l.split(" ")
         hit = splitHit[0]
         domain = splitHit[1]
         if hit in hitDomainDict:
            if domain not in hmm_list:
               hitDomainDict.pop(hit)
         l = f.readline()
   f.close()
   return hitDomainDict


def domainBoundaries(seq,hit,hit_domain_dict):
   domainList = hit_domain_dict[hit]
   newSeq = ""
   for i in domainList:
      domainBeginning = i.split('-')[0]
      domainEnding = i.split('-')[1]
      newSeq += seq[int(domainBeginning):int(domainEnding)]
   return newSeq


def mk_target_db(output,domainOut,hmm_name,target_seqs,alignment):
   target_db_fn = "output/"+output+"/"+alignment+"/"+domainOut +"/MSA/"+ hmm_name + '_targets.fasta'
   with open(target_db_fn,'w') as fh:
      for seq_rec in target_seqs:
         print(seq_rec.format('fasta'),end='',file=fh)
   return target_db_fn


def do_hmmalign(output,domainOut,hmm_name,hmm_file,target_db_fn,alignment):
   stk_out_fn = "output/"+output+"/"+alignment+"/"+domainOut +"/MSA/"+hmm_name + '.stk'
   subprocess.call(['hmmalign','-o',stk_out_fn,hmm_file,target_db_fn])
   return stk_out_fn


def do_esl_alimask(output,domainOut,hmm_name,stk_fn,alignment):
   alimask_out_fn = "output/"+output+"/"+alignment+"/"+domainOut +"/MSA/"+hmm_name + '.alimask'
   subprocess.call(['esl-alimask','-o',alimask_out_fn,
                    '--pavg','0.5','-p',stk_fn])
   return alimask_out_fn


def do_esl_alimanip(output,domainOut,hmm_name,alimask_fn,alignment):
  alimanip_out_fn = "output/"+output+"/"+alignment+"/"+domainOut +"/MSA/"+hmm_name + '.alimanip'
  print("YOOOOO")
  subprocess.call(['esl-alimanip','-o',alimanip_out_fn,
                    '--minpp','0.3',alimask_fn])

   # need to do this sequentially - first remove bad residues to gap
   # then filter out sequences that are now too short
   # seem to be able to write to same name as input

  subprocess.call(['esl-alimanip','-o',alimanip_out_fn,
                  '--lnfract','0.7',alimanip_out_fn])
  return alimanip_out_fn


def convert_to_mfa(output,domainOut,hmm_name,alimanip_fn,alignment):
   mfa_out_fn = "output/"+output+"/"+alignment+"/"+domainOut +"/MSA/"+hmm_name + '_auto.mfa'
   alignment = AlignIO.read(alimanip_fn,'stockholm')
   AlignIO.write(alignment,mfa_out_fn,'fasta')
   return mfa_out_fn


def filter_mfa(output,domainOut,hmm_name,mfa_fn,alignment):
   filtered_out_fn = "output/"+output+"/"+alignment+"/"+domainOut +"/MSA/"+hmm_name + '_auto_filtered.mfa'
   seq_hash = defaultdict(set)
   reps = []
   for seq_rec in SeqIO.parse(mfa_fn,'fasta'):
      sp_tag = seq_rec.id.split('|')[0]
      seq_str = str(seq_rec.seq)
      if seq_str in seq_hash:
         if sp_tag not in seq_hash[seq_str]:
            seq_hash[seq_str].add(sp_tag)
            reps.append(seq_rec)
      else:
         seq_hash[seq_str].add(sp_tag)
         reps.append(seq_rec)
   with open(filtered_out_fn,'w') as fh:
      for seq_rec in reps:
         print(seq_rec.format('fasta'),end='',file=fh)
   return filtered_out_fn
