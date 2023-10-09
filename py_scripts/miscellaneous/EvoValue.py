import numpy as np


def getGeneSeq(proteomPath,species,geneName):
	with open(proteomPath+"/"+species+"Proteins.fasta","r") as f:
		l = f.readline()
		seq = ""
		while l != "":
			if l[0] == ">" and geneName in l:
				l = f.readline()
				while l[0] != ">":
					seq += l
					l = f.readline()
			else:
				l = f.readline()
	f.close()
	seq = seq.replace("\n")
	seq = seq.split()
	return seq


def valueBases():
	surfaceA=40
	surfaceT=32
	surfaceC=14
	surfaceG=54
	dictSurface={}
	dictSurface['A'] = surfaceA
	dictSurface['T'] = surfaceT
	dictSurface['C'] = surfaceC
	dictSurface['G'] = surfaceG
	return dictSurface
	
def geneCompo(dictSurface,geneSeq):
		surfaceGene=0
		for i in geneSeq:
			surfaceGene+=dictSurface[i]
		return surfaceGene

def statistics(surfaceList):
	std = np.std(surfaceList)
	mean = np.mean(surfaceList)
	var = np.var(surfaceList)
	return mean,var,std


def geneGroup(dictSurface,groupGenes):
	surfaceList=[]
	for i in groupGenes:
		surfaceGene = geneCompo(dictSurface,i)
		surfaceList.append(surfaceGene)
	mean,var,std = statistics(surfaceList)
    
def testStat(dictSurface,listGroupGenes):
	mean1,var1,std1=geneGroup(dictSurface,listGroupGenes[0])
	mean2,var2,std2=geneGroup(dictSurface,listGroupGenes[1])
	(t_stat, p) = np.ranksums(mean1,mean2)
	print(p)
	return p
							
						
										
def main():
	geneSeq=['A','T','C']
	dictSurface = valueBases()
	surfaceGene=geneCompo(dictSurface,geneSeq)
	print(surfaceGene)
	
main()