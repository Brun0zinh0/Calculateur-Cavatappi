# -*- coding: utf-8 -*-
"""Inspection brute des fichiers expérimentaux de Sacha."""
import pandas as pd
import numpy as np

BASE = r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents\Stage muscle artificièle\modèle\modèle chinois\Espace de travail\expérimentale\Sacha"

files_xlsx = [
    "Maintient sous pression D 20 mn.xlsx",
    "Mise sous pression muscle B 200 g.xlsx",
    "Variation rappide de position (D).xlsx",
]
files_csv = [
    "Mise sous pression muscle I 200 g.csv",
    "Mise sous pression muscle J 200 g.csv",
]

for f in files_xlsx:
    path = BASE + "\\" + f
    print("=" * 80)
    print("FICHIER:", f)
    xl = pd.ExcelFile(path)
    print("  feuilles:", xl.sheet_names)
    for sh in xl.sheet_names:
        df = xl.parse(sh, header=None)
        print(f"  -- feuille '{sh}': shape={df.shape}")
        print(df.head(12).to_string())
        print("  ...")
        print(df.tail(5).to_string())

for f in files_csv:
    path = BASE + "\\" + f
    print("=" * 80)
    print("FICHIER:", f)
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        lines = fh.readlines()
    print("  nb lignes:", len(lines))
    for l in lines[:12]:
        print("  |", l.rstrip())
    print("  ...")
    for l in lines[-5:]:
        print("  |", l.rstrip())
