# Gestion de projet & Prog avancée - Alignement Multiple CLUSTAL
# MAGROUN Lina (M2 Bioinformatique) - Septembre 2026

# SCRIPT : alignment.py
# Rôle : Fonctions de calcul de l'alignement (matrice, arbre guide et alignement progressif).
#
# Deux adaptations sont faites conformément au sujet :
#   1. On remplace l'arbre hiérarchique UPGMA par un embranchement séquentiel glouton.
#   2. On remplace l'heuristique de Wilbur & Lipman par une programmation
#      dynamique exacte (Needleman-Wunsch avec gaps affines de Gotoh).
#
# Les autres choix de score et de pénalités correspondent aux paramètres
# retenus pour cette implémentation.

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Matrice Dayhoff PAM 250 brute, tirée de l'article.
# Les valeurs d'origine (log-odds) vont de -8 à +17.
# L'article indique qu'il faut ramener ces valeurs sur une échelle positive de 0 à 25.
# Pour cela, on ajoutera +8 à chaque valeur lors du parsing.
PAM250_RAW = """
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
A  2 -2  0  0 -2  0  0  1 -1 -1 -2 -1 -1 -3  1  1  1 -6 -3  0
R -2  6  0 -1 -4  1 -1 -3  2 -2 -3  3  0 -4  0  0 -1  2 -1 -2
N  0  0  2  2 -4  1  1  0  2 -2 -3  1 -2 -4  0  1  0 -4 -2 -2
D  0 -1  2  4 -5  2  3  0  1 -2 -4  0 -3 -6 -1  0  0 -7 -4 -2
C -2 -4 -4 -5 12 -5 -5 -3 -3 -2 -6 -5 -5 -4 -3  0 -2 -8  0 -2
Q  0  1  1  2 -5  4  2 -1  3 -2 -2  1 -1 -5  0 -1 -1 -5 -4 -2
E  0 -1  1  3 -5  2  4 -1  1 -2 -3  0 -2 -5 -1  0  0 -7 -4 -2
G  1 -3  0  0 -3 -1 -1  5 -2 -3 -4 -2 -3 -5  0  1  0 -7 -5 -1
H -1  2  2  1 -3  3  1 -2  6 -2 -2  0 -2 -2  0 -1 -1 -3  0 -2
I -1 -2 -2 -2 -2 -2 -2 -3 -2  5  2 -2  2  1 -2 -1  0 -5 -1  4
L -2 -3 -3 -4 -6 -2 -3 -4 -2  2  6 -3  4  2 -3 -3 -2 -2  1  2
K -1  3  1  0 -5  1  0 -2  0 -2 -3  5  0 -5 -1  0  0 -3 -4 -2
M -1  0 -2 -3 -5 -1 -2 -3 -2  2  4  0  6  0 -2 -1 -1 -4 -2  2
F -3 -4 -4 -6 -4 -5 -5 -5 -2  1  2 -5  0  9 -5 -3 -3  0  7 -1
P  1  0  0 -1 -3  0 -1  0  0 -2 -3 -1 -2 -5  6  1  0 -6 -5 -1
S  1  0  1  0  0 -1  0  1 -1 -1 -3  0 -1 -3  1  2  1 -2 -3 -1
T  1 -1  0  0 -2 -1  0  0 -1  0 -2  0 -1 -3  0  1  3 -5 -3  0
W -6  2 -4 -7 -8 -5 -7 -7 -3 -5 -2 -3 -4  0 -6 -2 -5 17  0 -6
Y -3 -1 -2 -4  0 -4 -4 -5  0 -1  1 -4 -2  7 -5 -3 -3  0 10 -2
V  0 -2 -2 -2 -2 -2 -2 -1 -2  4  2 -2  2 -1 -1 -1  0 -6 -2  4
"""


def parse_pam250():
    """
    On convertit la matrice brute PAM250 en un dictionnaire Python
    (clé = tuple d'acides aminés, valeur = score ajusté).

    On lance ce parsing une seule fois lors de l'importation du script
    pour éviter de ralentir les calculs suivants (voir la variable globale).
    """
    lines = [l.strip() for l in PAM250_RAW.strip().split('\n') if l.strip()]
    headers = lines[0].split()

    matrix = {}
    for line in lines[1:]:
        parts = line.split()
        row_aa = parts[0]
        scores = [int(x) for x in parts[1:]]
        # On parcourt chaque score de la ligne pour l'associer à sa colonne correspondante
        for col_aa, score in zip(headers, scores):
            # On ajoute +8 pour n'avoir que des valeurs positives 
            matrix[(row_aa, col_aa)] = score + 8
    return matrix


# On charge la matrice en mémoire dès l'importation du script
PAM250 = parse_pam250()


def parametres_gap_par_defaut(type_seq):
    """
    On définit les pénalités d'ouverture et d'extension de gap par défaut.

    Les valeurs (10/2 pour l'ADN, 25/10 pour les protéines) sont adaptées 
    aux échelles de score respectives (0-10 pour l'ADN, 0-25 pour PAM250) 
    pour éviter une insertion excessive de gaps artificiels.
    """
    if type_seq == "adn":
        return 10, 2
    return 25, 10


def get_score(a, b, type_seq='proteine'):
    """
     On évalue la distance (ou coût de substitution) entre deux résidus 'a' et 'b'.
    Comme notre algorithme cherche à minimiser la distance globale:

    - ADN : On utilise le modèle à 3 distances (Identité = 0, Transition = 5, 
      Transversion = 10).
    - Protéine : On convertit la matrice de similarité en matrice de distance 
      par la formule "25 - score" pour ramener les meilleurs scores près de 0.
    """
    # Si l'un des caractères est un gap, on retourne 0 car les pénalités 
    # de gap (ouverture et extension) sont gérées à part dans l'algorithme de Gotoh (matrices Ix et Iy)
    if a == '-' or b == '-':
        return 0

    if type_seq == 'adn':
        if a == b:
            return 0
        purines = {'A', 'G'}
        pyrimidines = {'C', 'T'}
        if (a in purines and b in purines) or (a in pyrimidines and b in pyrimidines):
            return 5
        return 10
    else:
        return 25 - PAM250.get((a, b), 0)


def alignement_dynamique(seq1, seq2, type_seq='proteine', gap_open=None, gap_extend=None):
    """
    On réalise l'alignement global optimal de deux séquences (Needleman-Wunsch) en 
    intégrant le modèle de pénalités affines de Gotoh.
    
    Le principe des gaps affines est d'avoir un coût d'ouverture de trou (gap_open) 
    qui pénalise la création d'une nouvelle délétion/insertion, et un coût d'extension 
    (gap_extend) plus faible pour la prolonger, ce qui est biologiquement plus réaliste.

    Pour cela on maintient trois matrices distinctes :
      - M  : On y gère les alignements où les deux séquences avancent ensemble (match ou mismatch).
      - Ix : On y stocke les scores lorsqu'un gap est ouvert ou prolongé dans la séquence 2.
      - Iy : On y stocke les scores lorsqu'un gap est ouvert ou prolongé dans la séquence 1.
    """
    if gap_open is None or gap_extend is None:
        defaut_open, defaut_extend = parametres_gap_par_defaut(type_seq)
        gap_open = defaut_open if gap_open is None else gap_open
        gap_extend = defaut_extend if gap_extend is None else gap_extend

    s1 = seq1["residus"]
    s2 = seq2["residus"]
    n = len(s1)
    m = len(s2)
    
    # On utilise l'infini positif car on cherche à minimiser la distance globale
    INF = float("inf")

    # On prépare nos trois matrices de programmation dynamique
    M = [[INF] * (m + 1) for _ in range(n + 1)]
    Ix = [[INF] * (m + 1) for _ in range(n + 1)]
    Iy = [[INF] * (m + 1) for _ in range(n + 1)]

    # On fixe notre point de départ à un coût nul (aucune lettre consommée)
    M[0][0] = 0

    # On initialise la première ligne et la première colonne : cela correspond 
    # à un gap ouvert dès le début de l'alignement puis prolongé tout le long.
    for i in range(1, n + 1):
        Ix[i][0] = gap_open + (i - 1) * gap_extend
    for j in range(1, m + 1):
        Iy[0][j] = gap_open + (j - 1) * gap_extend

    # On remplit les trois matrices de programmation dynamique (phase Forward)
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            # On récupère le coût de substitution (distance) entre nos deux résidus
            d = get_score(s1[i - 1], s2[j - 1], type_seq)

            # M[i][j] correspond à la diagonale : on cherche la transition minimale parmi les 3 états précédents
            M[i][j] = min(M[i - 1][j - 1], Ix[i - 1][j - 1], Iy[i - 1][j - 1]) + d
            
            # On détermine s'il est plus rentable d'ouvrir un gap depuis M ou de prolonger un gap existant
            Ix[i][j] = min(M[i - 1][j] + gap_open, Ix[i - 1][j] + gap_extend)
            Iy[i][j] = min(M[i][j - 1] + gap_open, Iy[i][j - 1] + gap_extend)

    # Le score final optimal correspond à la valeur minimale trouvée dans la cellule en bas à droite de nos trois grilles
    score_final = min(M[n][m], Ix[n][m], Iy[n][m])

    # On détermine la matrice de départ en bas à droite pour lancer la remontée (backtracking)
    if M[n][m] <= Ix[n][m] and M[n][m] <= Iy[n][m]:
        etat = "M"
    elif Ix[n][m] <= Iy[n][m]:
        etat = "Ix"
    else:
        etat = "Iy"

    align1 = ""
    align2 = ""
    i, j = n, m

    # On remonte pas à pas le chemin optimal pour reconstruire les séquences alignées
    while i > 0 or j > 0:
        if etat == "M":
            d = get_score(s1[i - 1], s2[j - 1], type_seq)
            align1 = s1[i - 1] + align1
            align2 = s2[j - 1] + align2

            # On retranche le coût de substitution pour savoir de quelle matrice on provenait
            precedent = M[i][j] - d
            i -= 1
            j -= 1
            if precedent == M[i][j]:
                etat = "M"
            elif precedent == Ix[i][j]:
                etat = "Ix"
            else:
                etat = "Iy"

        elif etat == "Ix":
            align1 = s1[i - 1] + align1
            align2 = "-" + align2
            
            # Si le score correspond à une ouverture de gap, on repasse dans l'état de substitution M
            # Sinon, c'est que le gap a été prolongé, on reste donc dans l'état Ix
            etat = "M" if Ix[i][j] == M[i - 1][j] + gap_open else "Ix"
            i -= 1

        else:  # etat == "Iy"
            align1 = "-" + align1
            align2 = s2[j - 1] + align2
            
            # Même logique : retour à l'état M s'il s'agit d'une ouverture, sinon prolongation dans Iy
            etat = "M" if Iy[i][j] == M[i][j - 1] + gap_open else "Iy"
            j -= 1

    return align1, align2, score_final


def calculer_matrice_similarite(liste_sequences, type_seq='proteine', gap_open=None, gap_extend=None, penalite_gap=1):
    """
     On calcule la matrice de similarité par paires (pairwise) qui servira de base pour l'arbre guide.
    
    Note : l'alignement de Gotoh calcule une distance que l'on cherche à minimiser.
    Pour déterminer l'ordre d'intégration, on utilise ensuite un score de similarité
    calculé à partir du nombre d'identités et des colonnes contenant des gaps :

    Score = Identités - (nb_gaps * pénalité)
    
    """
    n = len(liste_sequences)
    matrice_sim = [[0 for _ in range(n)] for _ in range(n)]

    # La matrice est symétrique, on ne calcule donc que la moitié supérieure
    for i in range(n):
        for j in range(i + 1, n):
            align1, align2, _ = alignement_dynamique(
                liste_sequences[i], liste_sequences[j], type_seq, gap_open, gap_extend
            )

            # On compte les correspondances exactes de résidus (hors gaps)
            identites = sum(1 for c1, c2 in zip(align1, align2) if c1 == c2 and c1 != "-")
            # On compte le nombre de colonnes contenant au moins un gap
            nb_gaps = sum(1 for c1, c2 in zip(align1, align2) if c1 == "-" or c2 == "-")

            # On applique le score combiné 
            score_sim = identites - nb_gaps * penalite_gap

            matrice_sim[i][j] = score_sim
            matrice_sim[j][i] = score_sim

    return matrice_sim


def determiner_ordre_alignement(matrice_sim):
    """
    Conçoit l'ordre d'intégration progressif (notre alternative gloutonne à UPGMA).
    
    La démarche consiste à identifier d'abord la paire ayant la plus forte similarité 
    pour initialiser l'alignement, puis à intégrer séquentiellement la séquence 
    restante la plus similaire à l'un des membres du groupe déjà formé.
    
    Note : Notre alignement Gotoh calcule une distance (plus c'est bas, plus elles sont proches),
    mais Clustal a besoin de scores de similarité pour l'arbre (plus le score est haut, mieux c'est).
    On applique la formule : Score = Identités - (gaps * pénalité).
    """
    n = len(matrice_sim)
    meilleur_score = float("-inf")
    seq_a, seq_b = 0, 1

    # Recherche du couple de départ optimal
    for i in range(n):
        for j in range(i + 1, n):
            if matrice_sim[i][j] > meilleur_score:
                meilleur_score = matrice_sim[i][j]
                seq_a, seq_b = i, j

    ordre = [seq_a, seq_b]
    visite = {seq_a, seq_b}

    # Intégration itérative des autres séquences
    while len(ordre) < n:
        meilleure_suivante = None
        meilleur_score_suivant = float("-inf")

        for k in range(n):
            if k not in visite:
                # On évalue la séquence par rapport à toutes celles déjà intégrées
                for deja_choisie in ordre:
                    if matrice_sim[k][deja_choisie] > meilleur_score_suivant:
                        meilleur_score_suivant = matrice_sim[k][deja_choisie]
                        meilleure_suivante = k

        ordre.append(meilleure_suivante)
        visite.add(meilleure_suivante)

    return ordre


def alignement_profil(seq, profil, type_seq='proteine', gap_open=None, gap_extend=None):
    """
    On aligne une séquence individuelle contre un profil (un groupe de séquences déjà alignées).
    
    Pour comparer la séquence au profil, on multiplie les pénalités de gap par le nombre 
    de séquences 'K' du profil. Cela nous évite d'utiliser des moyennes à virgule (flottants) 
    qui créent des erreurs d'arrondi lors du backtracking en Python. Les proportions restent identiques.
    
    Les gaps déjà présents dans le profil restent figés entre eux lors de la remontée.
    """
    if gap_open is None or gap_extend is None:
        defaut_open, defaut_extend = parametres_gap_par_defaut(type_seq)
        gap_open = defaut_open if gap_open is None else gap_open
        gap_extend = defaut_extend if gap_extend is None else gap_extend

    s = seq["residus"]
    n = len(s)
    m = len(profil[0]["residus"])
    k_sequences = len(profil)

    # On ajuste les pénalités de gap en les multipliant par la taille du profil
    p_gap_open = gap_open * k_sequences
    p_gap_extend = gap_extend * k_sequences

    INF = float("inf")
    M = [[INF] * (m + 1) for _ in range(n + 1)]
    Ix = [[INF] * (m + 1) for _ in range(n + 1)]
    Iy = [[INF] * (m + 1) for _ in range(n + 1)]

    M[0][0] = 0

    for i in range(1, n + 1):
        Ix[i][0] = p_gap_open + (i - 1) * p_gap_extend
    for j in range(1, m + 1):
        Iy[0][j] = p_gap_open + (j - 1) * p_gap_extend

    # On isole les colonnes du profil pour simplifier le calcul des distances
    colonnes_profil = [[p_seq["residus"][j] for p_seq in profil] for j in range(m)]

    def score_somme_colonne(residu, colonne):
        """
        On additionne les distances entre notre résidu et chaque lettre de la colonne.
        On ignore les gaps internes du profil car ils ont déjà été comptés.
        """
        somme = 0
        for char_profil in colonne:
            # On ne calcule pas la distance avec un gap pour ne pas le pénaliser deux fois
            if char_profil != "-":
                somme += get_score(residu, char_profil, type_seq)
        return somme

    # On calcule les scores dans nos trois matrices de programmation dynamique
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            # On utilise notre fonction pour sommer les distances avec les lettres de la colonne
            d = score_somme_colonne(s[i - 1], colonnes_profil[j - 1])

            # On applique les transitions de Gotoh (on prend le minimum car on cherche à minimiser un coût)
            M[i][j] = min(M[i - 1][j - 1], Ix[i - 1][j - 1], Iy[i - 1][j - 1]) + d
            Ix[i][j] = min(M[i - 1][j] + p_gap_open, Ix[i - 1][j] + p_gap_extend)
            Iy[i][j] = min(M[i][j - 1] + p_gap_open, Iy[i][j - 1] + p_gap_extend)

    # Le meilleur score final correspond au minimum tout en bas à droite des trois matrices
    score_final = min(M[n][m], Ix[n][m], Iy[n][m])

    # On cherche l'état de départ en bas à droite pour commencer la remontée
    if M[n][m] <= Ix[n][m] and M[n][m] <= Iy[n][m]:
        etat = "M"
    elif Ix[n][m] <= Iy[n][m]:
        etat = "Ix"
    else:
        etat = "Iy"
    
    # On prépare nos chaînes vides pour reconstruire l'alignement à l'envers
    align_seq = ""
    align_profil = ["" for _ in profil]
    i, j = n, m

    # On remonte le chemin pour reconstruire l'alignement entre la séquence et le profil
    while i > 0 or j > 0:
        if etat == "M":
            d = score_somme_colonne(s[i - 1], colonnes_profil[j - 1])
            align_seq = s[i - 1] + align_seq
            
            # On ajoute la colonne entière du profil existant à notre alignement
            for k in range(k_sequences):
                align_profil[k] = colonnes_profil[j - 1][k] + align_profil[k]

            # On retrouve la cellule précédente en soustrayant le coût de substitution
            precedent = M[i][j] - d
            i -= 1
            j -= 1
            # On décale sur la diagonale (i et j décrémentés) et on vérifie de quelle matrice on venait
            if precedent == M[i][j]:
                etat = "M"
            elif precedent == Ix[i][j]:
                etat = "Ix"
            else:
                etat = "Iy"

        elif etat == "Ix":
            align_seq = s[i - 1] + align_seq
            # On ajoute un gap à chaque séquence du profil
            for k in range(k_sequences):
                align_profil[k] = "-" + align_profil[k]
            # On teste si ce gap vient d'être ouvert (on repasse dans M) ou s'il s'agit d'une extension (on reste dans Ix)
            etat = "M" if Ix[i][j] == M[i - 1][j] + p_gap_open else "Ix"
            i -= 1

        else:  # etat == "Iy"
            align_seq = "-" + align_seq
            # On ajoute un gap dans la nouvelle séquence
            for k in range(k_sequences):
                align_profil[k] = colonnes_profil[j - 1][k] + align_profil[k]
            # Même logique pour savoir si on ferme ou si on continue le gap dans la séquence
            etat = "M" if Iy[i][j] == M[i][j - 1] + p_gap_open else "Iy"
            j -= 1

    # On crée notre liste finale avec les séquences du profil et la nouvelle séquence alignées
    nouveau_profil = []
    for k in range(k_sequences):
        nouveau_profil.append({"nom": profil[k]["nom"], "residus": align_profil[k]})
    nouveau_profil.append({"nom": seq["nom"], "residus": align_seq})

    return nouveau_profil


def lancer_alignement_multiple(liste_sequences, ordre, type_seq='proteine', gap_open=None, gap_extend=None):
    """
    On réalise l'alignement progressif : on commence par aligner les deux premières 
    séquences, puis on ajoute les autres au profil au fur et à mesure.
    """

    # On récupère les deux séquences les plus proches selon notre arbre guide
    seq_a = liste_sequences[ordre[0]]
    seq_b = liste_sequences[ordre[1]]

    # On fait le premier alignement par paires avec notre algorithme de Gotoh
    align1, align2, _ = alignement_dynamique(seq_a, seq_b, type_seq, gap_open, gap_extend)

    # On crée notre premier profil avec cet alignement initial de deux séquences
    profil = [
        {"nom": seq_a["nom"], "residus": align1},
        {"nom": seq_b["nom"], "residus": align2},
    ]

    # On ajoute les autres séquences une par une au profil existant (à partir du 3e élément)
    for index in ordre[2:]:
        nouvelle_seq = liste_sequences[index]
        profil = alignement_profil(nouvelle_seq, profil, type_seq, gap_open, gap_extend)

    return profil


# FONCTIONS D'EXPORT ET DE VISUALISATION 

def construire_arbre_newick(noms_ordonnes):
    """
    On écrit l'arbre au format Newick standard d'après notre ordre d'alignement.
    Comme on ajoute les séquences une par une de manière gloutonne, l'arbre a 
    une structure en peigne (chaque nœud se greffe sur le bloc précédent).
    """
    if len(noms_ordonnes) < 2:
        return ""

    # On commence par le couple de séquences les plus proches d'après notre matrice
    arbre = f"({noms_ordonnes[0]},{noms_ordonnes[1]})"

    # On intègre les autres séquences une par une autour de ce premier noyau
    for nom in noms_ordonnes[2:]:
        arbre = f"({arbre},{nom})"
    # On ferme l'arbre avec le point-virgule requis par le format Newick
    return arbre + ";"


def generer_arbre_ascii(noms_ordonnes):
    """
    On génère un dessin de l'arbre en caractères ASCII pour pouvoir visualiser 
    la structure en peigne directement dans la console.
    """
    n = len(noms_ordonnes)
    if n < 2:
        return ""

    # On dessine la première séparation verticale reliant nos deux premières séquences
    lines = [
        "      +--- " + noms_ordonnes[0],
        "   +--|",
        "   |  +--- " + noms_ordonnes[1]
    ]

    # Pour chaque séquence suivante, on allonge le schéma vers la droite et le bas
    for idx in range(2, n):
        # On décale le dessin existant vers la droite pour faire de la place aux nouvelles branches
        for i in range(len(lines)):
            lines[i] = "   |  " + lines[i]

        # On dessine le nouvel embranchement vertical
        lines.append("   +--|")
        # On y rattache la séquence que l'on vient d'ajouter au profil
        lines.append("      +--- " + noms_ordonnes[idx])

    return "\n".join(lines)


def determiner_ordre_alignement_avec_scores(matrice_sim):
    """
    On cherche l'ordre d'alignement des séquences par embranchement glouton, tout en 
    conservant les scores de similarité de chaque fusion pour les afficher sur notre dendrogramme.
    """
    n = len(matrice_sim)

    # On utilise l'infini négatif pour trouver un maximum (on cherche le score le plus haut)
    meilleur_score = float("-inf")
    seq_a, seq_b = 0, 1

    # On cherche d'abord les deux séquences les plus proches (qui se ressemblent le plus)
    for i in range(n):
        for j in range(i + 1, n):
            if matrice_sim[i][j] > meilleur_score:
                meilleur_score = matrice_sim[i][j]
                seq_a, seq_b = i, j

    # On initialise notre ordre d'alignement avec ces deux séquences de départ
    ordre = [seq_a, seq_b]

    # On garde le score de leur ressemblance pour notre futur graphique
    scores = [meilleur_score]

    # On note que ces séquences ont déjà été traitées
    visite = {seq_a, seq_b}

    # On cherche à connecter toutes les autres séquences restantes, une par une
    while len(ordre) < n:
        meilleure_suivante = None
        meilleur_score_suivant = float("-inf")

        # On cherche la séquence non visitée qui ressemble le plus à une des séquences déjà choisies
        for k in range(n):
            if k not in visite:
                for deja_choisie in ordre:
                    if matrice_sim[k][deja_choisie] > meilleur_score_suivant:
                        meilleur_score_suivant = matrice_sim[k][deja_choisie]
                        meilleure_suivante = k

        # On ajoute la séquence sélectionnée et son score à nos listes de suivi
        ordre.append(meilleure_suivante)
        scores.append(meilleur_score_suivant)
        visite.add(meilleure_suivante)

    return ordre, scores


def calculer_statistiques_alignement(profil):
    """
    On calcule les statistiques de notre alignement final (nombre de gaps, taux de gaps) 
    pour mesurer précisément l'impact de nos paramètres de pénalités.
    """
    # Toutes les séquences étant alignées, elles ont exactement la même longueur
    longueur_alignement = len(profil[0]["residus"])
    stats = []

    for seq in profil:
        # On compte les tirets (gaps) introduits dans cette séquence
        nb_gaps = seq["residus"].count("-")

        # On en déduit la longueur d'origine sans les trous
        nb_residus = longueur_alignement - nb_gaps

        # On calcule le pourcentage de gaps par rapport à la longueur totale
        taux_gaps = (nb_gaps / longueur_alignement) * 100

        # On regroupe les résultats de la séquence dans notre dictionnaire
        stats.append({
            "nom": seq["nom"],
            "longueur_alignement": longueur_alignement,
            "nb_gaps": nb_gaps,
            "nb_residus": nb_residus,
            "taux_gaps": taux_gaps,
        })

    return stats


def dessiner_dendrogramme(ordre, scores, liste_sequences, chemin_sortie="dendrogramme.png"):
    """
    On génère le graphique de notre arbre guide (dendrogramme) avec Matplotlib,
    en affichant les scores de fusions successives pour justifier notre ordre d'intégration.
    """
    n = len(ordre)
    noms = [liste_sequences[i]["nom"] for i in ordre]
    
    # On fixe une hauteur Y unique pour chaque séquence sur notre graphique
    y_positions = {ordre[i]: i for i in range(n)}

    # On adapte la hauteur de l'image selon le nombre de séquences pour éviter les chevauchements
    fig, ax = plt.subplots(figsize=(9, 0.5 * n + 2))

    # Premier nœud de fusion ( noyau) : on se place au milieu vertical des deux premières séquences
    y_cluster = (y_positions[ordre[0]] + y_positions[ordre[1]]) / 2

    # On dessine vers la gauche (X négatifs) pour que les feuilles soient à 0 (à droite) et la racine à gauche
    ax.plot([0, -1], [y_positions[ordre[0]], y_positions[ordre[0]]], color="black")
    ax.plot([0, -1], [y_positions[ordre[1]], y_positions[ordre[1]]], color="black")

    # On trace la barre verticale de connexion
    ax.plot([-1, -1], [y_positions[ordre[0]], y_positions[ordre[1]]], color="black")

    # On affiche notre premier score de similarité juste à côté de l'embranchement
    ax.text(-1.05, y_cluster, f"{scores[0]}", fontsize=7, va="center", ha="right", color="gray")

    x_precedent = -1
    
    # On ajoute les séquences restantes une par une vers la gauche
    for k in range(2, n):
        nouvelle_feuille = ordre[k]
        x_nouveau = -k
        y_nouvelle = y_positions[nouvelle_feuille]

        # On trace les branches horizontales et la connexion verticale du nouveau groupe
        ax.plot([0, x_nouveau], [y_nouvelle, y_nouvelle], color="black")
        ax.plot([x_precedent, x_nouveau], [y_cluster, y_cluster], color="black")
        ax.plot([x_nouveau, x_nouveau], [y_cluster, y_nouvelle], color="black")

        # On écrit le score de cette nouvelle fusion
        ax.text(x_nouveau - 0.05, (y_cluster + y_nouvelle) / 2, f"{scores[k - 1]}",
                fontsize=7, va="center", ha="right", color="gray")

        # On recalcule la hauteur moyenne pour le prochain embranchement
        y_cluster = (y_cluster * k + y_nouvelle) / (k + 1)
        x_precedent = x_nouveau

    # On ajuste la mise en page (on masque les axes inutiles pour ne garder que les noms à gauche)
    ax.set_yticks(range(n))
    ax.set_yticklabels(noms, fontsize=9)
    ax.set_xlabel("Ordre d'embranchement séquentiel (en gris : score de similarité à chaque fusion)")
    ax.set_title("Arbre par embranchement séquentiel")

    # On retire le cadre pour l'esthétique 
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.set_xticks([])

    plt.tight_layout()
    plt.savefig(chemin_sortie, dpi=150)
    plt.close()