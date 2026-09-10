# Gestion de projet & Prog avancée - Alignement Multiple CLUSTAL
# MAGROUN Lina (M2 Bioinformatique) - Septembre 2026

# SCRIPT : parser.py
# Rôle : Gère la lecture des fichiers au format FASTA et l'écriture des résultats d'alignement.

def lire_fasta(chemin_fichier):
    """
    Lit un fichier au format FASTA ligne par ligne et retourne une liste contenant 
    les séquences sous forme de dictionnaires avec leur nom et leurs résidus.
    """
    liste_sequences = [] 
    
    with open(chemin_fichier, 'r') as fichier:
        nom_actuel = ""
        sequence_actuelle = ""
        
        for ligne in fichier:
            # On nettoie la ligne (suppression des espaces et sauts de ligne \r/\n)
            ligne = ligne.strip()
            
            # On détecte le début d'une nouvelle séquence (commence par '>')
            if ligne.startswith(">"):
                # On enregistre la séquence précédente si elle existe
                if nom_actuel != "":
                    liste_sequences.append({"nom": nom_actuel, "residus": sequence_actuelle})
                
                # On récupère le nom de la séquence (en sautant le caractère '>')
                nom_actuel = ligne[1:]
                sequence_actuelle = ""
                
            else:
                # Sinon, on accumule les lignes de la séquence (utile si elle est coupée sur plusieurs lignes)
                sequence_actuelle += ligne
                
        # On sauvegarde la toute dernière séquence du fichier à la fin de la boucle
        if nom_actuel != "":
            liste_sequences.append({"nom": nom_actuel, "residus": sequence_actuelle})
            
    return liste_sequences


def ecrire_fasta(liste_sequences, chemin_fichier_sortie):
    """
    Exporte le profil d'alignement multiple final dans un nouveau fichier au format FASTA.
    """
    with open(chemin_fichier_sortie, 'w') as fichier:
        for seq in liste_sequences:
            fichier.write(f">{seq['nom']}\n")
            fichier.write(f"{seq['residus']}\n")
