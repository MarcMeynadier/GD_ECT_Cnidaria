import argparse
import subprocess,os
import numpy as np
import regex as re
import pandas as pd
from copy import deepcopy
from itertools import chain
from nltk import flatten
from domainRetrieve import *



def verifyTopology(orthogroups,output,domain):
  with open("output/"+output+"/MSA/"+domain+"Tree.tree","r") as f:
    newickTree = f.readlines()
  f.close()
  for i in pco:
    for j in i:
            


pco = pcoOrthogroups()
verifyTopology(pco,"inputClytiaCnidocytes","pcoOrthologs")
