# Projet MPRO
Projet: optimisation robuste du problème de tournée de véhicule

Binôme: Tadeo Delapalme et Dimitri de Saint Guilhem

rappel deadlines:
 - 19/11 création du dépôt github
 - 15/12 modélisation papier
 - 15/02 rapport
 - 14/02 soutenance.

## Dossiers
- data: instances
- fig: figures des sous-tours pour des instances générées.
- illustrations: figures pour la présentation et le rapport
- rapports_soutenance: contient la modélisation papier, le rapport ainsi que la soutenance
- results: contient les résultats des méthodes exactes sur 600s
- results_heuristiques: contient les résultats des heuristiques.

## Notebook
- illustrations.ipynb: extraction et visualisation des resultats des heuristiques
- output.ipynb: extraction et visualisation des performances des différentes méthodes.

## Code Julia
- dual.jl: implémente l'algorithme de dualisation.
- export_results.jl
- heurstiques.jl
- lecture.jl
- main_exactes.jl: exécute les méthodes exactes de résolution avec le même temps d'exécution pour chaque.
- optim.jl:
- plans_coupants.jl: implémente l'algorithme des plans_coupants et de branch-and-cut.
- static.jl: résout le problème statique.
-test.jl: