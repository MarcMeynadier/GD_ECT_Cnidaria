def modifyHydractiniaProteom():
    with open("../phyloMarkers/input/proteins/metazoanProteins.fasta","r") as f:
        proteom = f.readlines()
    f.close()
    for i, valeur in enumerate(proteom):
        if valeur[0] == ">" and "Hy|" in valeur:
            nouvelle_valeur = valeur.replace("\n","g\n")
            proteom[i] = nouvelle_valeur
    with open("../phyloMarkers/input/proteins/metazoanProteinsNew.fasta","w") as f:
        for i in proteom:
            f.write(i)
    f.close()


modifyHydractiniaProteom()
    