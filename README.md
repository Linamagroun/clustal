# Projet CLUSTAL - Alignement Multiple Progressif

**Matière :** Gestion de projet & Programmation avancée  
**Auteur :** MAGROUN Lina (M2 Bioinformatique)  
**Date :** Septembre 2026  

## Présentation du projet

Ce projet consiste à réaliser un alignement multiple de séquences en s'inspirant de la méthode **CLUSTAL** décrite par Higgins et Sharp en 1989. L'objectif est de construire progressivement un alignement à partir d'une matrice de similarité.

Deux adaptations ont été faites par rapport à la méthode originale pour répondre au sujet :
* L'arbre UPGMA est remplacé par un **embranchement séquentiel glouton**.
* L'heuristique de Wilbur et Lipman est remplacée par une **programmation dynamique exacte avec gaps affines (Gotoh)**.

## Organisation du projet

Le code est réparti dans les fichiers suivants :
* `main.py` : Lance l'ensemble du programme.
* `parser.py` : Lecture et écriture des fichiers au format FASTA.
* `alignment.py` : Contient les fonctions de calcul (scores, matrice de similarité, arbre guide et alignement progressif).
* `benchmark.py` : Tests des temps d'exécution.
* `sequences.fasta` : Le jeu de séquences utilisé pour les tests.

## Choix techniques et gestion des scores

Pour que l'algorithme fonctionne correctement, plusieurs choix ont été faits concernant les calculs :

* **Matrice PAM250 :** Les scores bruts de Dayhoff (qui vont de -8 à +17) sont décalés de +8 au moment de la lecture pour ne garder que des entiers positifs (de 0 à 25). Vu que le programme cherche le chemin le plus court (minimisation d'une distance), on applique la formule `25 - score`. Deux acides aminés identiques ont donc un coût de substitution proche de 0.
* **Système de scores ADN :** On applique le barème à trois niveaux de l'article : 0 pour une identité, 5 pour une transition et 10 pour une transversion.
* **Pénalités de gaps :** Les distances ADN vont de 0 à 10, alors que les protéines vont jusqu'à 25. Pour éviter que le programme ne remplisse les protéines de trous artificiels, l'échelle des pénalités d'ouverture et d'extension a été adaptée au type biologique (10/2 pour l'ADN, 25/10 pour les protéines). Ces valeurs sont un choix d'implémentation.
* **Algorithme de Gotoh :** Il gère l'ouverture et l'extension des gaps séparément via trois matrices de programmation dynamique (M, I_x, I_y). Contrairement à l'article qui évoque l'optimisation mémoire linéaire de Myers & Miller, ce script utilise des matrices complètes.
* **Comparaison Séquence-Profil :** Insérer un gap face à un profil de $K$ séquences multiplie simplement la pénalité par K. Cela permet de rester sur des calculs d'entiers et d'éviter les bugs d'arrondis liés aux `float` lors de la phase de remontée. Les gaps du profil sont conservés tels quels.
* **Matrice de similarité :** L'ordre d'alignement est déterminé par un score calculé selon la formule de l'article : `Identités - (nb_gaps * pénalité)`. La pénalité est ici fixée à 1.

## Fonctionnement

Le programme suit ces étapes :
1. Lecture des séquences FASTA.
2. Détermination automatique du type de séquences (ADN ou protéines).
3. Calcul de la matrice de similarité entre toutes les paires.
4. Détermination de l'ordre d'embranchement (glouton : on commence par la meilleure paire, puis on ajoute la séquence la plus proche de celles déjà sélectionnées).
5. Alignement progressif (séquence contre séquence, puis séquence contre profil).
6. Calcul des statistiques finales (pourcentage de gaps).
7. Sauvegarde du résultat dans `resultat_alignement.fasta`.
8. Génération du graphique de l'arbre guide dans `dendrogramme.png`.

## Résultats obtenus

Testé sur les 10 séquences de la protéine p53 fournies, l'alignement final contient **411 colonnes**. 

L'ordre d'intégration obtenu est :
`p53_Humain → p53_Macaque → p53_Lapin → p53_Rat → p53_Souris → p53_Vache → p53_Chien → p53_Poulet → p53_Poisson_Zebre → p53_Cheval`

*Note :* La séquence du cheval, récupérée telle quelle sur UniProt, est naturellement plus courte (280 résidus). Cela a permis de vérifier empiriquement que l'algorithme gère très bien les grands gaps terminaux et la place correctement en groupe externe. Après vérification par suppression des gaps, les séquences de l'alignement final correspondent parfaitement aux originaux.

## Benchmark

Le fichier `benchmark.py` permet d'étudier l'évolution du temps de calcul :
* en fonction de la longueur des séquences.
* en fonction du nombre de séquences.

Les résultats montrent une augmentation rapide du temps de calcul (la complexité explose avec la taille des matrices), ce qui est cohérent avec l'étude expérimentale de l'article d'origine.

## Exécution

### Prérequis

Le programme utilise la bibliothèque graphique `matplotlib` pour générer l'arbre. Avant de lancer le script, assurez-vous de l'installer avec :
```bash
pip3 install matplotlib
```

Pour lancer le programme principal :
```bash
python3 main.py
```

Pour lancer les tests de performances temporelles :
```bash
python3 benchmark.py
```
