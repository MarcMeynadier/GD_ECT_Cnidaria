import matplotlib.pyplot as plt

# Standard library packages
import os

# Import Numpy, Pandas and Seaborn
import numpy as np
import pandas as pd
import seaborn as sns

# Import biopython tools for running local BLASTX
from Bio.Blast.Applications import NcbiblastpCommandline

# Colour scale transformation
from matplotlib.colors import LogNorm



def inputSpecies(sp1,sp2):
    speciesStr = sp1+'_'+sp2
    rPath = os.getcwd()
    s1 = rPath+'/../../../../species/'+sp1+'/raw/longest_isoform_'+sp1+'Proteins.fasta'
    s2 = rPath+'/../../../../species/'+sp2+'/raw/longest_isoform_'+sp2+'Proteins.fasta'
    outdir = rPath+'/../../../../species/crossSpecies/'+speciesStr+'/RBH/'
    os.makedirs(outdir, exist_ok=True)
    # Define output BLAST results
    fwd_out = os.path.join(outdir, '05-fwd-results.tab')
    rev_out = os.path.join(outdir, '05-rev-results.tab')
    return s1,s2,fwd_out,rev_out



def blastCommand(s1,s2,fwd_out,rev_out):
    blastp='/Users/mmeynadier/Documents/PhD/scripts/tools/ncbi-blast/bin/blastp'
    fwd_blastp = NcbiblastpCommandline(cmd=blastp,query=s1, subject=s2, out=fwd_out,
                                    outfmt="6 qseqid sseqid pident qcovs qlen slen length bitscore evalue",
                                    max_target_seqs=1)
    rev_blastp = NcbiblastpCommandline(cmd=blastp,query=s2, subject=s1, out=rev_out,
                                    outfmt="6 qseqid sseqid pident qcovs qlen slen length bitscore evalue",
                                    max_target_seqs=1)
    # Inspect command-lines
    print("FORWARD: %s" % fwd_blastp)
    print("REVERSE: %s" % rev_blastp)
    # Run BLAST searches
    fwd_stdout, fwd_stderr = fwd_blastp()
    rev_stdout, rev_stderr = rev_blastp()
    # Check STDOUT, STDERR
    print("FWD STDOUT: %s" % fwd_stdout)
    print("FWD STDERR: %s" % fwd_stderr)
    print("REV STDOUT: %s" % rev_stdout)
    print("REV STDERR: %s" % rev_stderr)


def discardDuplicate(filePath):
    df = pd.read_csv(filePath,sep='\t',header=None)
    df.columns = ['qseqid','sseqid','pident','qcovs','qlen','slen','length','bitscore','evalue']
    df2 = df.groupby(['qseqid', 'sseqid']).first()
    return df2

def verifyPresence(sp1,sp2):
    path = "../../../../species/crossSpecies/"
    returnFlag = True
    if os.path.exists(path+sp1+"_"+sp2) or os.path.exists(path+sp2+"_"+sp1):
        if os.path.exists(path+sp1+"_"+sp2+'/RBH/') or os.path.exists(path+sp2+"_"+sp1+'/RBH/'): 
            returnFlag = True
        else:
            os.makedirs(path+sp1+"_"+sp2+'/RBH/')
            returnFlag = False
    else:  
        returnFlag = False
        os.makedirs(path+sp1+"_"+sp2)
        os.makedirs(path+sp1+"_"+sp2+'/RBH/')
    return returnFlag


def main():
    rPath = os.getcwd()
    with open('speciesList.txt','r') as f:
        sp = f.readline()
    f.close()
    sp = sp.replace("\n","")
    spList = sp.split(",")
    for i in spList:
        for j in spList:
            if i != j:
                returnFlag = verifyPresence(i,j)
                if returnFlag == False:
                    print("Processing RBH for "+i+" & "+j+" species.\n")
                    s1,s2,fwd_out,rev_out = inputSpecies(i,j)
                    blastCommand(s1,s2,fwd_out,rev_out)
                    try:
                        forwardCured = discardDuplicate(rPath+'/../../../../species/crossSpecies/'+i+'_'+j+'/RBH/05-rev-results.tab') 
                        reverseCured = discardDuplicate(rPath+'/../../../../species/crossSpecies/'+i+'_'+j+'/RBH/05-fwd-results.tab')
                        forwardCured.to_csv(rPath+'/../../../../species/crossSpecies/'+i+'_'+j+'/RBH/'+i+'_'+j+'.txt',sep='\t')
                        reverseCured.to_csv(rPath+'/../../../../species/crossSpecies/'+i+'_'+j+'/RBH/'+j+'_'+i+'.txt',sep='\t')
                    except:
                        forwardCured = discardDuplicate(rPath+'/../../../../species/crossSpecies/'+j+'_'+i+'/RBH/05-rev-results.tab') 
                        reverseCured = discardDuplicate(rPath+'/../../../../species/crossSpecies/'+j+'_'+i+'/RBH/05-fwd-results.tab')
                        forwardCured.to_csv(rPath+'/../../../../species/crossSpecies/'+j+'_'+i+'/RBH/'+i+'_'+j+'.txt',sep='\t')
                        reverseCured.to_csv(rPath+'/../../../../species/crossSpecies/'+j+'_'+i+'/RBH/'+j+'_'+i+'.txt',sep='\t')
    return
    
                

main()
