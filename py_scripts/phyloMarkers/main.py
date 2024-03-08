#!/usr/bin/env python

import argparse


from MultipleSequenceAlignmentModule import *
from PhylogeneticModule import *
from SingleCellModule import *

###########################################

def main(args):
  alignment = args.alignment ; output = args.output ; query = args.query ; mode = args.mode ; domains = args.domain
  if query != None and domains != None:
    print("It is mandatory to choose between a FASTA input and a list of proteic domains.")
  elif query != None:
    queriesList = extractQueriesFromFasta('input/queries/'+query)
    for query in queriesList:
      hmm_list = getQueryDomains(query)
      if alignment == "HMM":
        domainOut,msaInput = processHMM(hmm_list,output,alignment,mode)
      elif alignment == "MAFFT":
        domainOut,msaInput = processMAFFT(hmm_list,alignment,output)  
      else:
        print("Please choose a valid multiple alignment method\n")
        exit()
      dictTaxo,paralogsFamilies = processPhylogenetic(output,domainOut,msaInput,alignment) 
      scExpr(output,domainOut,dictTaxo,paralogsFamilies,alignment)
  elif domains != None:
    hmm_list = domains.split(',')
    if alignment == "HMM":
      domainOut,msaInput = processHMM(hmm_list,output,alignment,mode)
    elif alignment == "MAFFT":
      domainOut,msaInput = processMAFFT(hmm_list,alignment,output)  
    else:
      print("Please choose a valid multiple alignment method\n")
      exit()
    dictTaxo,paralogsFamilies = processPhylogenetic(output,domainOut,msaInput,alignment) 
    scExpr(output,domainOut,dictTaxo,paralogsFamilies,alignment)

###########################################

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = '')
    parser.add_argument('-a','--alignment',action='store',dest='alignment',help='Alignment method') 
    parser.add_argument('-o','--output',action='store',dest='output',help='Output name')
    parser.add_argument('-q','--query',action='store',dest='query',help='Protein query name')
    parser.add_argument('-d','--domain',action='store',dest='domain',help='Domain name')  
    parser.add_argument('-m','--mode',action='store',dest='mode',help='Mode name')   
    args = parser.parse_args()

    print(args)
    main(args)