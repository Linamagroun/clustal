# Gestion de projet & Prog avancée - Alignement Multiple CLUSTAL
# MAGROUN Lina (M2 Bioinformatique) - Septembre 2026

# SCRIPT : benchmark.py
## Rôle : Tester la vitesse du programme en modifiant la longueur et le nombre de séquences.
# Permet de tracer des graphiques pour voir comment le temps de calcul augmente.
#
# Ce script permet de mesurer les performances de notre code. 
# Le but est de vérifier qu'on retrouve bien les mêmes logiques que dans l'article 
# de Higgins & Sharp (1989) : 
# - Le temps augmente avec le carré de la longueur à cause de la programmation dynamique.
# - Le calcul de la matrice de similarité devient l'étape limitante quand le nombre de séquences augmente.

import time
import random
from alignment import calculer_matrice_similarite, determiner_ordre_alignement, lancer_alignement_multiple


def generer_sequence_fictive(longueur, type_seq="proteine"):
    """
    Génère une séquence aléatoire de la longueur voulue pour faire nos mesures 
    sans avoir besoin de charger de vrais fichiers à chaque test.
    """
    if type_seq == "adn":
        alphabet = ['A', 'C', 'T', 'G']
    else:
        # Les 20 acides aminés dans l'ordre de notre matrice PAM250
        alphabet = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    return "".join(random.choice(alphabet) for _ in range(longueur))


def test_longueur_sequences():
    """
    On mesure le temps total d'exécution du pipeline complet en faisant varier 
    la longueur (L) des séquences, avec un nombre de séquences fixé à N = 5.
    
    On cherche à observer l'évolution du temps de calcul lorsque la longueur
    des séquences augmente. Vu que la programmation dynamique 
    remplit une grille de taille n*m pour chaque paire, le temps de calcul 
    augmente plus vite que la taille des séquences (il est proportionnel au carré).
    """
    print("\n--- TEST 1 : Temps d'alignement global en fonction de la LONGUEUR des séquences (N = 5) ---")
    print("Longueur (L)\tTemps de calcul total (secondes)")
    
    # On teste différentes tailles pour observer l'augmentation du temps en fonction de la longueur
    longueurs = [50, 100, 200, 300]
    nb_sequences = 5

    for L in longueurs:
        # On génère un jeu de séquences pour la longueur en cours
        set_test = []
        for idx in range(nb_sequences):
            set_test.append({
                "nom": f"seq_{idx}",
                "residus": generer_sequence_fictive(L, "proteine")
            })

        # On lance le chronomètre juste avant de démarrer les calculs
        debut = time.time()

        # on fait tourner le pipeline : matrice, ordre, alignement final
        matrice = calculer_matrice_similarite(set_test, type_seq="proteine")
        ordre = determiner_ordre_alignement(matrice)
        _ = lancer_alignement_multiple(set_test, ordre, type_seq="proteine")

        # On affiche le temps total pour cette longueur
        duree = time.time() - debut
        print(f"{L} aa/nt\t\t{duree:.4f} s")


def test_nombre_sequences():
    """
    On mesure séparément le temps de chaque étape du pipeline en faisant varier 
    le nombre de séquences (N), avec une longueur fixée à L = 50.
    
    On cherche à observer l'évolution du temps de calcul lorsque le nombre
    de séquences augmente. Quand le nombre de séquences grandit, 
    c'est le calcul de la matrice de similarité qui prend le plus de temps 
    parce qu'il y a beaucoup de paires à comparer.
    """
    print("\n--- TEST 2 : Décomposition du temps de calcul selon le NOMBRE de séquences (L = 50) ---")
    print("Nb Séquences (N)\tMatrice (s)\tArbre (s)\tAlignement (s)\tTotal (s)")

    tailles_set = [3, 5, 8, 12]
    longueur_fixe = 50

    for N in tailles_set:
        # On génère le jeu de séquences pour la taille N en cours
        set_test = []
        for idx in range(N):
            set_test.append({
                "nom": f"seq_{idx}",
                "residus": generer_sequence_fictive(longueur_fixe, "proteine")
            })

        # Étape 1 : On mesure le temps de calcul de la matrice de similarité (O(N^2 * L^2) au total)
        t_debut_matrice = time.time()
        matrice = calculer_matrice_similarite(set_test, type_seq="proteine")
        t_matrice = time.time() - t_debut_matrice

        # Étape 2 : On mesure le temps de construction de l'ordre d'embranchement
        t_debut_arbre = time.time()
        ordre = determiner_ordre_alignement(matrice)
        t_arbre = time.time() - t_debut_arbre

        # Étape 3 : On mesure le temps de l'alignement progressif final
        t_debut_align = time.time()
        _ = lancer_alignement_multiple(set_test, ordre, type_seq="proteine")
        t_align = time.time() - t_debut_align

        # On additionne tout pour avoir le total et on affiche la ligne de résultats
        total = t_matrice + t_arbre + t_align
        print(f"{N} séquences\t\t{t_matrice:.4f}\t\t{t_arbre:.4f}\t\t{t_align:.4f}\t\t{total:.4f}")


if __name__ == "__main__":
    print("=" * 80)
    print("DEBUT DU BENCHMARK COMPARATIF (ESPRIT FIGURES 1 & 2 DE HIGGINS & SHARP 1989)")
    print("=" * 80)

    test_longueur_sequences()
    test_nombre_sequences()

    print("\n" + "=" * 80)
    print("FIN DU BENCHMARK - DONNÉES PRÊTES POUR LES FIGURES DU RAPPORT")
    print("=" * 80)
