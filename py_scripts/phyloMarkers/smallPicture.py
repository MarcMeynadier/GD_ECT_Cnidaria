import argparse
import subprocess,os
import numpy as np
import regex as re
import pandas as pd
import os.path


def subsetTranscriptDomain(domain,dt_out):
    geneList = []
    with open(dt_out,"r") as f:
        l = f.readline()
        while l != "":
            if l[0] != "#" and " "+domain+" " in l:
                gene = l.split(" ")[0]
                geneList.append(gene)
                l = f.readline()
            else:
                l = f.readline()
    f.close()
    geneList2 = []
    for i in geneList:
        if i not in geneList2:
            geneList2.append(i)
    return geneList2


def domainTranscriptSeq(domainTranscript,proteomFasta):
    domainTranscriptSeq = []
    indexlist = []
    with open(proteomFasta,"r") as f:
        l = f.readline()
        while l != "":
            if l[0] == ">":
                g = l.replace(">","") ; g = g.replace("\n","")
                prot = ""
                if g in domainTranscript:
                    indexlist.append(domainTranscript.index(g))
                    l = f.readline()
                    while l[0] != ">":
                        prot += l
                        l = f.readline()
                    domainTranscriptSeq.append(prot)
                else:
                    l = f.readline()
            else:
                l = f.readline()
    f.close()
    domainTranscriptSeq[:] = [domainTranscriptSeq[i] for i in indexlist]
    return domainTranscriptSeq