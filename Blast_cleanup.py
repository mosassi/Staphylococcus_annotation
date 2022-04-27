# -*- coding: utf-8 -*-
"""
Created on Wed Apr 4 2018

@author: Mohamed Sassi, Dimitri Pédron
"""
import os
import pandas as panda
import sys

nArg=0
for Parameter in sys.argv:
    print Parameter
    if "-h" in Parameter:
        print("\nScript d'analyse des ARN régulateurs\n")
        print("In  : fichiers de sortie de Blast")
        print("Out : matrice presence/absence totale format csv")
        print("      fichiers fasta contenant les ARN")
        print("      ")
        print("\nParameter\n")
        print("\n-h : Affiche l'aide")
        print("\n-i : Définit l'input")
        print("\n-o : Définit la sortie")
        print("\n-db : Définit la base de données")
        print("\n-r : Désactive les fichiers de se sortie pour R")
        print("\n-b : (Bavard) Le script affiche la progression en console")
        print("\n-a : (Auromatique) Le script est appelé par le script précédent")
    if "-b" in Parameter:
        modeBavard=True
    if "-r" in Parameter:
        exportR=False
    if "-i" in Parameter:
        input_directory=sys.argv[nArg+1]
        #Where the blast output file is located
    if "-o" in Parameter:
        output_directory=sys.argv[nArg+1]
        #Where the cleaned up blast output and matrix will be output
    if "-db" in Parameter:
        database=sys.argv[nArg+1]
        #Defines the database used to
    nArg+=1

def blast_output_list(input_directory):
    listeFichiers=[]
    for filename in os.walk(input_directory):
        for files in filename[2]:
            if files[0]!=".": #excludes invisible files in Mac
                if "sortieBlast_" in files:
                    files=files.replace("sortieBlast_","")
                    if files not in listeFichiers:
                        listeFichiers.append(files)
    return listeFichiers

def ref_genome_list(listeFichiers):
    listeGenomes=[]
    for noms in listeFichiers:
        print(noms)
        genome=noms.replace(".txt","")
        print(genome)
        if "SRD" in genome:
            genome=genome.replace("SRD_","")
        if "_uid" in genome:
            genome=genome.split("_uid")[0]
        if genome not in listeGenomes:
            listeGenomes.append(genome)
    return listeGenomes

def cleanup_blast(listeFichiers, col):
    """ lit les fichiers de sortie de blast et les nettoie """
    cleanup_directory = output_directory+"Clean/"
    if not os.path.exists(cleanup_directory):
        os.makedirs(cleanup_directory)

    listeGenomesLus=[]
    r=0 #compteur du nombre de fichier npour l'affichage console
    bestFin=0
    for file in listeFichiers:
        dictArnIndices={}
        r+=1
        genome=file.replace(".txt","") #produit le nom du genome

        #montre l'avancement dans la console
        if modeBavard==True:
            print(genome+"\t"+str(r)+"/"+str(len(listeFichiers)))

        #prepare les sorties
        if genome not in listeGenomesLus:
            listeGenomesLus.append(genome)

        #ouvre le fichier produit par Blast
        sortieBlast=open(input_directory+"sortieBlast_"+file, "r")
        tableau=panda.read_csv(sortieBlast , sep="\t", header=None, names=col) #lit le tableau , index_col="NomARN"
        i=tableau["EValue"].count()#compte le nombre de lignes dans le tableau
        tableau = tableau.assign(Coverage=tableau['align_length'] / tableau['query_length'])

        #ouvre le fichier de sortie avec les ARNs netttoyés
        fileClean=open(cleanup_directory+"clean_"+file, "w")

        #lit l'ensemble du fichier et associe dans un dictionnaire chaque ARN
        # avec la liste des indices de lignes correpondant a cet arn
        #permet de gagner en temps d'execution en ne relisant tout le fichier qu'une fois
        for x in range(0,i):
            complete_name=tableau.index[x]
            if database == "ProphageFinder":
                res_type=complete_name.split('_')[0]
            elif database == "PlasmidFinder":
                res_type=complete_name.split('_')[-1]

            if res_type in dictArnIndices:
                dictArnIndices[res_type].append(x)
            else:
                dictArnIndices[res_type]=[]
                dictArnIndices[res_type].append(x)

        #debut de la filtration
        #itere sur chaque ARN sans relire l'ensemble du tableau
        for arn in dictArnIndices:
           listeIndices=dictArnIndices[arn]
           if len(listeIndices)>1:
               arn_table=tableau.iloc[listeIndices]
               arn_table.sort_values(by=['Coverage', 'Identity'])
               fileClean.write(arn_table.iloc[[0]].to_csv(sep="\t", header=False, index=True))
           else:
               arn_index = listeIndices[0]
               fileClean.write(tableau.iloc[[arn_index]].to_csv(sep="\t", header=False, index=True))


def main():
        """ ouvre les principaux fichiers de sortie """
        #si les fichiers ne peuvent être écrits le script plante...
        #fichierSortie2=open(dossier+"matriceIndex.csv", "w")
        #fichierSortie=open(dossier+"matrice.csv", "w")

        #en tetes des colonnes servant  pour lire les fichiers de sortie de blast
        col=["NomARN", "Score", "EValue", "Identity", "DebutRef", "FinRef", "SizeRef", "Debut", "Fin", "query_length", "align_length", "Sequence"]
        #titres des colonnes pour la seconde matrice
        colMatrice2=["Start", "End", "Length", "Sens"]
        print("Step1")
        listeFichiers = blast_output_list(input_directory)
        print("Step2")
        listeGenomes = ref_genome_list(listeFichiers)
        print("Step3")
        cleanup_blast(listeFichiers, col)



if __name__ == "__main__":
    main()
