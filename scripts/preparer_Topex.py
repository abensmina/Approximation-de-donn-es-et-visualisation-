#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Projet GM4 -- Prétraitement des données Topex (tâche T7, BF1 à BF3)
==================================================================

Entrée : le fichier texte téléchargé sur https://topex.ucsd.edu/cgi-bin/get_data.cgi
         (trois colonnes : longitude, latitude, profondeur en mètres).

Ce que fait le script :
  1. lit le fichier en ignorant les lignes non numériques ;
  2. ramène les longitudes dans [-180, 180] ;
  3. nettoie : valeurs manquantes, doublons, valeurs hors de [-11 500, 9 000] m ;
  4. sépare les points contraints par un sondage de navire (profondeur impaire)
     des points prédits par la gravimétrie (profondeur paire) ;
  5. vérifie la limite du code (41 000 points) et sous-échantillonne si besoin ;
  6. écrit les fichiers au format du code, les réponses du code et un aperçu.

Exemple :
    python3 preparer_topex.py topex_brut.txt
    python3 preparer_topex.py topex_brut.txt --km --sous-ech 2 --maillage 20 --eps 1e-6

Fichiers produits :
    donnees.dat          toutes les données nettoyées (format splic1)
    sondages.dat         uniquement les points contraints par des sondages
    entrees_donnees.txt  réponses du code pour donnees.dat
    entrees_sondages.txt réponses du code pour sondages.dat
    apercu_zone.png      carte des données et des sondages
    resume_zone.txt      statistiques de la zone (à recopier dans le rapport)
"""

import argparse
import sys

import numpy as np

MAX_POINTS = 41000          # limite MXNP du code Fortran
Z_MIN, Z_MAX = -11500.0, 9000.0


def lire_topex(nom):
    lignes = []
    with open(nom) as f:
        for ligne in f:
            morceaux = ligne.replace(",", " ").split()
            if len(morceaux) < 3:
                continue
            try:
                lignes.append([float(m) for m in morceaux[:3]])
            except ValueError:
                continue                     # en-tête ou ligne de texte
    if not lignes:
        sys.exit(f"Aucune ligne numérique trouvée dans {nom}.")
    return np.array(lignes)


def nettoyer(d):
    rapport = {"lus": len(d)}
    d = d[np.all(np.isfinite(d), axis=1)]
    rapport["valeurs manquantes"] = rapport["lus"] - len(d)

    lon = d[:, 0].copy()
    lon[lon > 180] -= 360
    d[:, 0] = lon

    _, idx = np.unique(np.round(d[:, :2], 6), axis=0, return_index=True)
    avant = len(d)
    d = d[np.sort(idx)]
    rapport["doublons"] = avant - len(d)

    ok = (d[:, 2] >= Z_MIN) & (d[:, 2] <= Z_MAX)
    rapport["aberrantes"] = int((~ok).sum())
    d = d[ok]
    rapport["conservés"] = len(d)
    return d, rapport


def sous_echantillonner(d, k):
    """Garde un point sur k dans chaque direction de la grille."""
    if k <= 1:
        return d
    lons, lats = np.unique(d[:, 0]), np.unique(d[:, 1])
    keep_lon = set(np.round(lons[::k], 6))
    keep_lat = set(np.round(lats[::k], 6))
    m = np.array([(round(a, 6) in keep_lon) and (round(b, 6) in keep_lat)
                  for a, b in d[:, :2]])
    return d[m]


def en_km(d, lon0, lat0):
    """Projection locale équirectangulaire : x, y en km autour de (lon0, lat0)."""
    x = (d[:, 0] - lon0) * 111.32 * np.cos(np.radians(lat0))
    y = (d[:, 1] - lat0) * 110.57
    return np.column_stack([x, y, d[:, 2]])


def ecrire_donnees(nom, d):
    with open(nom, "w") as f:
        f.write(f"{len(d)}\n")
        for x, y, z in d:
            f.write(f"{x:.6f} {y:.6f} {z:.1f}\n")


def ecrire_reponses(nom, donnees, resultat, domaine, maillage, eps, npev):
    xmin, ymin, xmax, ymax = domaine
    lignes = ["o", xmin, ymin, xmax, ymax, maillage, maillage, donnees,
              "o", f"{eps:g}", "o", resultat, npev, npev]
    with open(nom, "w") as f:
        f.write("\n".join(str(l) for l in lignes) + "\n")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("fichier", help="fichier brut téléchargé sur Topex")
    p.add_argument("--sous-ech", type=int, default=1,
                   help="garder 1 point sur k par direction (défaut : automatique si > 41 000 points)")
    p.add_argument("--km", action="store_true",
                   help="convertir longitude/latitude en kilomètres (mise à l'échelle)")
    p.add_argument("--maillage", type=int, default=20, help="intervalles par direction (<= 30)")
    p.add_argument("--eps", type=float, default=1e-6, help="paramètre de lissage")
    p.add_argument("--npev", type=int, default=150, help="points d'évaluation par direction")
    args = p.parse_args()

    d, rapport = nettoyer(lire_topex(args.fichier))

    # Sous-échantillonnage pour respecter la limite du code
    k = max(1, args.sous_ech)
    while len(sous_echantillonner(d, k)) > MAX_POINTS:
        k += 1
    if k > 1:
        print(f"Sous-échantillonnage : 1 point sur {k} par direction.")
    d = sous_echantillonner(d, k)

    lon_min, lon_max = d[:, 0].min(), d[:, 0].max()
    lat_min, lat_max = d[:, 1].min(), d[:, 1].max()
    lon0, lat0 = (lon_min + lon_max) / 2, (lat_min + lat_max) / 2

    impair = (np.round(d[:, 2]).astype(int) % 2) != 0
    s = d[impair]

    if args.km:
        d_out, s_out = en_km(d, lon0, lat0), en_km(s, lon0, lat0)
    else:
        d_out, s_out = d, s
    # Domaine englobant, très légèrement élargi : le code ignore sans message
    # les points situés hors du domaine.
    marge_x = 1e-4 * (d_out[:, 0].max() - d_out[:, 0].min())
    marge_y = 1e-4 * (d_out[:, 1].max() - d_out[:, 1].min())
    domaine = (round(d_out[:, 0].min() - marge_x, 6), round(d_out[:, 1].min() - marge_y, 6),
               round(d_out[:, 0].max() + marge_x, 6), round(d_out[:, 1].max() + marge_y, 6))

    ecrire_donnees("donnees.dat", d_out)
    ecrire_donnees("sondages.dat", s_out)
    ecrire_reponses("entrees_donnees.txt", "donnees.dat", "resultat_donnees.dat",
                    domaine, args.maillage, args.eps, args.npev)
    ecrire_reponses("entrees_sondages.txt", "sondages.dat", "resultat_sondages.dat",
                    domaine, args.maillage, args.eps, args.npev)

    unite = "km" if args.km else "degrés"
    resume = [
        "Résumé de la zone d'étude",
        "=========================",
        f"Longitudes : {lon_min:.4f} à {lon_max:.4f}",
        f"Latitudes  : {lat_min:.4f} à {lat_max:.4f}",
        f"Profondeur : {d[:, 2].min():.0f} m à {d[:, 2].max():.0f} m",
        f"Points lus : {rapport['lus']}",
        f"  valeurs manquantes retirées : {rapport['valeurs manquantes']}",
        f"  doublons retirés            : {rapport['doublons']}",
        f"  valeurs aberrantes retirées : {rapport['aberrantes']}",
        f"Points conservés (après sous-échantillonnage 1/{k}) : {len(d)}",
        f"Points contraints par des sondages (impairs) : {len(s)} "
        f"({100 * len(s) / max(len(d), 1):.1f} %)",
        f"Domaine pour le code ({unite}) : xmin={domaine[0]}, ymin={domaine[1]}, "
        f"xmax={domaine[2]}, ymax={domaine[3]}",
    ]
    if args.km:
        resume.append(f"Projection locale centrée en lon={lon0:.4f}, lat={lat0:.4f}")
    if len(s) == 0 or len(s) == len(d):
        resume.append("ATTENTION : la règle pair/impair ne semble pas respectée "
                      "dans ce fichier ; vérifier la documentation de la grille.")
    with open("resume_zone.txt", "w") as f:
        f.write("\n".join(resume) + "\n")
    print("\n".join(resume))

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(12, 5))
        sc = ax[0].scatter(d_out[:, 0], d_out[:, 1], c=d_out[:, 2], s=2, cmap="viridis")
        fig.colorbar(sc, ax=ax[0], label="profondeur (m)")
        ax[0].set_title(f"Toutes les données ({len(d)} points)")
        ax[1].scatter(d_out[:, 0], d_out[:, 1], c="0.85", s=1)
        ax[1].scatter(s_out[:, 0], s_out[:, 1], c="k", s=2)
        ax[1].set_title(f"Sondages de navires ({len(s)} points)")
        for a in ax:
            a.set_xlabel("x (km)" if args.km else "longitude")
            a.set_ylabel("y (km)" if args.km else "latitude")
            a.set_aspect("equal" if args.km else 1 / np.cos(np.radians(lat0)))
        fig.tight_layout()
        fig.savefig("apercu_zone.png", dpi=130)
        print("Aperçu enregistré dans apercu_zone.png")
    except ImportError:
        print("matplotlib absent : pas d'aperçu.")

    print("\nÉtape suivante :  ./splic1 < entrees_donnees.txt")


if __name__ == "__main__":
    main()