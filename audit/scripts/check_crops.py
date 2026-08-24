"""Verifie la plausibilite de la numerisation : taille de l'image figure 7 et crops."""
import sys
from pathlib import Path

from PIL import Image, ImageDraw
from pypdf import PdfReader
from io import BytesIO

WORKSPACE = Path(r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents\Stage muscle artificièle\modèle\modèle chinois\Espace de travail")
SCRATCH = Path(r"C:\Users\B3DCB~1.PER\AppData\Local\Temp\claude\C--Users-b-pereiraazevedo-OneDrive---House-Of-HR-NV-Documents-Stage-muscle-artifici-le-mod-le-mod-le-chinois-Espace-de-travail\c5852e10-6d38-4703-b218-a1649a5b12c9\scratchpad")
PDF = WORKSPACE / "Article support" / "Blocked actuation of twisted and coiled tube-based polymer actuators driven by hydraulic pressure.pdf"

reader = PdfReader(str(PDF))
page = reader.pages[11]
print(f"page 12 du PDF: {len(page.images)} image(s) embarquee(s)")
img = Image.open(BytesIO(page.images[0].data)).convert("RGB")
print(f"taille image figure 7: {img.size}")

PANELS = {
    "a_force": (120, 51, 704, 350),
    "a_torque": (120, 351, 704, 647),
    "a_pressure": (120, 648, 704, 946),
    "b_force": (899, 51, 1484, 350),
    "b_torque": (899, 351, 1484, 647),
    "b_pressure": (899, 648, 1484, 946),
}
overlay = img.copy()
d = ImageDraw.Draw(overlay)
for name, (x0, y0, x1, y1) in PANELS.items():
    d.rectangle([x0, y0, x1, y1], outline=(255, 0, 0), width=3)
out = SCRATCH / "fig7_crop_overlay.png"
overlay.save(out)
print(f"overlay sauvegarde: {out}")
