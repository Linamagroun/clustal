# Gestion de projet & Prog avancée - Alignement Multiple CLUSTAL
# MAGROUN Lina (M2 Bioinformatique) - Septembre 2026

# SCRIPT : main.py
# Rôle : Script principal pour exécuter notre alignement CLUSTAL.
# On identifie le type de séquence (ADN ou protéine) pour adapter le système de score.
# On enchaîne ici les différentes étapes : lecture FASTA, détection du type, calcul 
# de la matrice, construction de l'arbre guide et alignement progressif.

import os
from parser import lire_fasta, ecrire_fasta
from alignment import (
    calculer_matrice_similarite,
    determiner_ordre_alignement_avec_scores,
    lancer_alignement_multiple,
    construire_arbre_newick,
    generer_arbre_ascii,
    dessiner_dendrogramme,
    calculer_statistiques_alignement
)


def determiner_type_sequences(liste_sequences):
    """
    Analyse les résidus des séquences pour identifier s'il s'agit d'ADN 
    ou de protéines, afin d'adapter automatiquement le système de score.
    """
    nucleotides = {'A', 'T', 'C', 'G', 'N', 'U', '-'}
    for seq in liste_sequences:
        residus = seq["residus"].upper()
        # Dès qu'un seul résidu ne fait pas partie de l'alphabet ADN,
        # on arrête la vérification : c'est une protéine.
        if not all(char in nucleotides for char in residus):
            return "proteine"
    return "adn"


def main():
    # Nom du fichier d'entrée contenant nos séquences à aligner
    chemin_fasta = "sequences.fasta"

    # On génère un jeu de test par défaut si le fichier "sequences.fasta" n'existe pas
    if not os.path.exists(chemin_fasta):
        with open(chemin_fasta, "w") as f:
            f.write(">Humain\nMALWMRLLPLLALLALWGPDPAAA\n>Souris\nMALLVHFLPLLALLALWEPKPTQA\n>Vache\nMALWTRLRPLLALLALWPPPPARA\n")

    print("--- DÉBUT DE L'ALIGNEMENT CLUSTAL (NEEDLEMAN-WUNSCH + GAP AFFINE DE GOTOH) ---")

    # ÉTAPE 1 : lecture et chargement des séquences depuis le fichier FASTA
    print("\n1. Lecture du fichier FASTA d'entrée...")
    mes_sequences = lire_fasta(chemin_fasta)
    print(f"-> {len(mes_sequences)} séquences chargées en mémoire.")
    for s in mes_sequences:
        print(f"   - [{s['nom']}] : {s['residus']}")

    # ÉTAPE 2 : On détecte automatiquement si on travaille sur de l'ADN ou
    # des protéines, pour utiliser le bon système de score (voir get_score
    # dans alignment.py) et les bonnes pénalités de gap par défaut
    type_seq = determiner_type_sequences(mes_sequences)
    print(f"-> Détection automatique du type biologique : {type_seq.upper()}")

    # ÉTAPE 3 : On calcule la matrice de similarité pairwise (2 à 2) qui servira 
    # à définir l'arbre guide. On ne précise pas gap_open/gap_extend ici : l'algorithme 
    # applique les valeurs adaptées selon le type de séquence (voir alignment.py)
    print(f"\n2. Calcul de la matrice de similarité (scoring : {type_seq.upper()})...")
    matrice = calculer_matrice_similarite(mes_sequences, type_seq=type_seq)
    print("Matrice obtenue (identités de résidus après alignement, moins pénalité de gap) :")
    for ligne in matrice:
        print(f"   {ligne}")

    # ÉTAPE 4 : On utilise cette matrice pour définir l'arbre guide et l'ordre 
    # dans lequel les séquences vont être intégrées (les paires les plus proches en premier)
    print("\n3. Calcul de l'ordre d'embranchement séquentiel...")
    ordre, scores_fusion = determiner_ordre_alignement_avec_scores(matrice)
    print(f"-> Ordre d'intégration calculé (par indices) : {ordre}")
    noms_ordonnes = [mes_sequences[idx]["nom"] for idx in ordre]
    print(f"-> Ordre d'alignement : {' -> '.join(noms_ordonnes)}")

    # On affiche l'arbre obtenu sous deux formats différents : le format standard Newick et le schéma ASCII
    print("\n[RÉSULTAT SUPPLÉMENTAIRE] Arbre d'embranchement construit :")
    print(f"  Format Newick : {construire_arbre_newick(noms_ordonnes)}")
    print("  Schéma ASCII de l'arbre guide :")
    print(generer_arbre_ascii(noms_ordonnes))

    # On génère et on sauvegarde le dendrogramme au format graphique PNG.
    # Il représente visuellement l'ordre d'embranchement et les scores
    # de similarité utilisés lors des différentes intégrations.
    chemin_dendrogramme = "dendrogramme.png"
    dessiner_dendrogramme(ordre, scores_fusion, mes_sequences, chemin_dendrogramme)
    print(f"  Dendrogramme sauvegardé dans '{chemin_dendrogramme}'.")

    # ÉTAPE 5 : On lance l'alignement multiple progressif selon l'ordre d'embranchement
    # On applique la programmation dynamique avec gaps affines (méthode de Gotoh)
    # pour aligner pas à pas les séquences et les profils
    print(f"\n4. Exécution de l'alignement multiple progressif (scoring : {type_seq.upper()})...")
    alignement_final = lancer_alignement_multiple(mes_sequences, ordre, type_seq=type_seq)
    print("Séquences alignées obtenues :")
    for s in alignement_final:
        print(f"   - [{s['nom']}] : {s['residus']}")

    # On calcule les statistiques de l'alignement final (résidus, gaps et proportions)
    print("\n[RÉSULTAT SUPPLÉMENTAIRE] Statistiques de l'alignement final :")
    stats = calculer_statistiques_alignement(alignement_final)
    print(f"   Longueur totale de l'alignement : {stats[0]['longueur_alignement']} colonnes")
    print(f"   {'Séquence':20s} {'Résidus':>8s} {'Gaps':>6s} {'% gaps':>8s}")
    for s in stats:
        print(f"   {s['nom']:20s} {s['nb_residus']:8d} {s['nb_gaps']:6d} {s['taux_gaps']:7.1f}%")

    
    # ÉTAPE 6 : On exporte l'alignement final dans un fichier au format FASTA
    print("\n5. Création du fichier de sortie...")
    chemin_sortie = "resultat_alignement.fasta"
    ecrire_fasta(alignement_final, chemin_sortie)
    print(f"-> Fichier de sortie '{chemin_sortie}' généré avec succès.")

    print("\n--- FIN DU PROGRAMME CLUSTAL ---")


if __name__ == "__main__":
    main()