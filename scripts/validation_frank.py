#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Projet GM4 -- Validation du code splic1 sur une fonction test (tâche T6, BF8)
============================================================================

On échantillonne la fonction de Franke en N points aléatoires de [0,1]^2,
on lance le code Fortran, puis on compare la surface obtenue à la vraie
fonction sur la grille d'évaluation.

Trois usages :

1) Préparer un cas (fichier de données + fichier de réponses) :
       python3 validation_franke.py preparer --n 500 --maillage 10 --eps 1e-6
   puis lancer le code :
       ./splic1 < entrees_test.txt

2) Calculer l'erreur d'un résultat déjà obtenu :
       python3 validation_franke.py erreur resultat_test.dat

3) Faire toute l'étude de convergence automatiquement
   (l'exécutable ./splic1, s1don et c1ref doivent être dans le dossier) :
       python3 validation_franke.py etude

Dépendances : numpy, matplotlib (pour la figure de l'étude).
"""

import argparse
import os
import subprocess
import sys
import time

import numpy as np


# ---------------------------------------------------------------------
#  Fonction test
# ---------------------------------------------------------------------
def franke(x, y):
    """Fonction de Franke, fonction test classique en approximation 2D."""
    t1 = 0.75 * np.exp(-((9 * x - 2) ** 2) / 4 - ((9 * y - 2) ** 2) / 4)
    t2 = 0.75 * np.exp(-((9 * x + 1) ** 2) / 49 - (9 * y + 1) / 10)
    t3 = 0.50 * np.exp(-((9 * x - 7) ** 2) / 4 - ((9 * y - 3) ** 2) / 4)
    t4 = -0.20 * np.exp(-((9 * x - 4) ** 2) - (9 * y - 7) ** 2)
    return t1 + t2 + t3 + t4


# ---------------------------------------------------------------------
#  Fichiers d'entrée du code
# ---------------------------------------------------------------------
def ecrire_donnees(nom, x, y, z):
    """Format attendu par splic1 : N sur la 1re ligne, puis x y z."""
    with open(nom, "w") as f:
        f.write(f"{len(x)}\n")
        for a, b, c in zip(x, y, z):
            f.write(f"{a:.10f} {b:.10f} {c:.10f}\n")


def ecrire_reponses(nom, fichier_donnees, fichier_resultat, maillage, eps,
                    npev=50, domaine=(0.0, 0.0, 1.0, 1.0)):
    """Réponses aux questions de splic1, dans l'ordre où il les pose."""
    xmin, ymin, xmax, ymax = domaine
    lignes = [
        "o",                    # ajustement de données ponctuelles
        xmin, ymin, xmax, ymax,  # domaine d'étude
        maillage, maillage,     # nombre d'intervalles en x et en y
        fichier_donnees,        # fichier des points
        "o", f"{eps:g}",        # valeur particulière du paramètre de lissage
        "o",                    # estimation sur une grille régulière
        fichier_resultat,       # fichier résultat (format Tecplot)
        npev, npev,             # nombre de points d'évaluation en x et en y
    ]
    with open(nom, "w") as f:
        f.write("\n".join(str(l) for l in lignes) + "\n")


# ---------------------------------------------------------------------
#  Lecture du résultat (format Tecplot écrit par EVGRI)
# ---------------------------------------------------------------------
def lire_resultat(nom):
    xs, ys, zs = [], [], []
    with open(nom) as f:
        for ligne in f:
            ligne = ligne.strip().replace("D", "E")
            if not ligne or ligne.upper().startswith(("VARIABLES", "ZONE")):
                continue
            morceaux = ligne.split()
            if len(morceaux) < 3:
                continue
            xs.append(float(morceaux[0]))
            ys.append(float(morceaux[1]))
            zs.append(float(morceaux[2]))
    return np.array(xs), np.array(ys), np.array(zs)


def erreurs(nom_resultat):
    """Erreur quadratique moyenne (RMSE) et erreur max par rapport à Franke."""
    x, y, s = lire_resultat(nom_resultat)
    if len(x) == 0:
        raise ValueError(f"Aucune valeur lue dans {nom_resultat}")
    e = s - franke(x, y)
    return np.sqrt(np.mean(e ** 2)), np.max(np.abs(e)), len(x)


# ---------------------------------------------------------------------
#  Commandes
# ---------------------------------------------------------------------
def cmd_preparer(args):
    rng = np.random.default_rng(args.graine)
    x, y = rng.random(args.n), rng.random(args.n)
    z = franke(x, y) + args.bruit * rng.standard_normal(args.n)
    donnees = f"test_N{args.n}.dat"
    ecrire_donnees(donnees, x, y, z)
    ecrire_reponses("entrees_test.txt", donnees, "resultat_test.dat",
                    args.maillage, args.eps, args.npev)
    print(f"Écrit : {donnees} ({args.n} points) et entrees_test.txt")
    print("Lance ensuite :  ./splic1 < entrees_test.txt")
    print("Puis :           python3 validation_franke.py erreur resultat_test.dat")


def cmd_erreur(args):
    rmse, emax, k = erreurs(args.fichier)
    print(f"{k} points d'évaluation")
    print(f"Erreur quadratique moyenne (RMSE) : {rmse:.3e}")
    print(f"Erreur maximale                   : {emax:.3e}")


def lancer_cas(exe, n, maillage, eps, npev, graine, bruit):
    rng = np.random.default_rng(graine)
    x, y = rng.random(n), rng.random(n)
    z = franke(x, y) + bruit * rng.standard_normal(n)
    donnees, resultat, reponses = "etude_donnees.dat", "etude_resultat.dat", "etude_entrees.txt"
    ecrire_donnees(donnees, x, y, z)
    ecrire_reponses(reponses, donnees, resultat, maillage, eps, npev)
    if os.path.exists(resultat):
        os.remove(resultat)
    t0 = time.time()
    with open(reponses) as entree:
        subprocess.run([exe], stdin=entree, stdout=subprocess.DEVNULL, check=True)
    duree = time.time() - t0
    rmse, emax, _ = erreurs(resultat)
    return rmse, emax, duree


def cmd_etude(args):
    exe = os.path.abspath(args.exe)
    for requis in (exe, "s1don", "c1ref"):
        if not os.path.exists(requis):
            sys.exit(f"Fichier introuvable : {requis} "
                     "(lance le script dans le dossier du code compilé).")

    lignes = ["etude,N,maillage,h,eps,rmse,emax,duree_s"]

    print("\n1) Convergence en h (N fixé)")
    print(f"{'maillage':>8} {'h':>8} {'RMSE':>11} {'Emax':>11} {'temps (s)':>10}")
    for m in args.maillages:
        rmse, emax, d = lancer_cas(exe, args.n, m, args.eps, args.npev, args.graine, args.bruit)
        h = 1.0 / m
        print(f"{m:>8d} {h:>8.4f} {rmse:>11.3e} {emax:>11.3e} {d:>10.1f}")
        lignes.append(f"h,{args.n},{m},{h},{args.eps},{rmse},{emax},{d}")

    print("\n2) Influence du nombre de points (maillage fixé)")
    print(f"{'N':>8} {'RMSE':>11} {'Emax':>11} {'temps (s)':>10}")
    for n in args.liste_n:
        rmse, emax, d = lancer_cas(exe, n, args.maillage_fixe, args.eps, args.npev, args.graine, args.bruit)
        print(f"{n:>8d} {rmse:>11.3e} {emax:>11.3e} {d:>10.1f}")
        lignes.append(f"N,{n},{args.maillage_fixe},{1/args.maillage_fixe},{args.eps},{rmse},{emax},{d}")

    with open("etude_validation.csv", "w") as f:
        f.write("\n".join(lignes) + "\n")
    print("\nRésultats enregistrés dans etude_validation.csv")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        data = [l.split(",") for l in lignes[1:]]
        hs = [(float(d[3]), float(d[5])) for d in data if d[0] == "h"]
        ns = [(int(d[1]), float(d[5])) for d in data if d[0] == "N"]
        fig, ax = plt.subplots(1, 2, figsize=(10, 4))
        ax[0].loglog([a for a, _ in hs], [b for _, b in hs], "o-")
        ax[0].set_xlabel("pas du maillage h")
        ax[0].set_ylabel("RMSE")
        ax[0].set_title(f"Convergence en h (N = {args.n})")
        ax[0].grid(True, which="both", alpha=0.3)
        ax[1].loglog([a for a, _ in ns], [b for _, b in ns], "s-")
        ax[1].set_xlabel("nombre de points N")
        ax[1].set_ylabel("RMSE")
        ax[1].set_title(f"Influence de N (maillage {args.maillage_fixe}x{args.maillage_fixe})")
        ax[1].grid(True, which="both", alpha=0.3)
        fig.tight_layout()
        fig.savefig("etude_validation.png", dpi=150)
        print("Figure enregistrée dans etude_validation.png")
    except ImportError:
        print("matplotlib absent : pas de figure.")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = p.add_subparsers(dest="commande", required=True)

    a = sp.add_parser("preparer", help="écrire un jeu de test et les réponses du code")
    a.add_argument("--n", type=int, default=500, help="nombre de points (défaut 500)")
    a.add_argument("--maillage", type=int, default=10, help="intervalles par direction, <= 30")
    a.add_argument("--eps", type=float, default=1e-6, help="paramètre de lissage")
    a.add_argument("--npev", type=int, default=50, help="points d'évaluation par direction")
    a.add_argument("--bruit", type=float, default=0.0, help="écart-type d'un bruit ajouté")
    a.add_argument("--graine", type=int, default=0)
    a.set_defaults(func=cmd_preparer)

    b = sp.add_parser("erreur", help="erreur d'un fichier résultat")
    b.add_argument("fichier")
    b.set_defaults(func=cmd_erreur)

    c = sp.add_parser("etude", help="étude complète (lance ./splic1)")
    c.add_argument("--exe", default="./splic1")
    c.add_argument("--n", type=int, default=2000, help="N pour l'étude en h")
    c.add_argument("--maillages", type=int, nargs="+", default=[5, 10, 15, 20, 25, 30])
    c.add_argument("--liste-n", type=int, nargs="+", default=[100, 250, 500, 1000, 2000, 4000])
    c.add_argument("--maillage-fixe", type=int, default=20)
    c.add_argument("--eps", type=float, default=1e-6)
    c.add_argument("--npev", type=int, default=50)
    c.add_argument("--bruit", type=float, default=0.0)
    c.add_argument("--graine", type=int, default=0)
    c.set_defaults(func=cmd_etude)

    args = p.parse_args()
    if getattr(args, "maillage", 1) > 30 or any(m > 30 for m in getattr(args, "maillages", [])):
        sys.exit("Le code accepte au plus 30 intervalles par direction.")
    args.func(args)


if __name__ == "__main__":
    main()