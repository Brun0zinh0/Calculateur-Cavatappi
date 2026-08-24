from io import BytesIO
from pathlib import Path

from PIL import Image
from pypdf import PdfReader

PDF = Path(
    r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents"
    r"\Stage muscle artificièle\modèle\modèle chinois\Espace de travail"
    r"\Article support"
    r"\Blocked actuation of twisted and coiled tube-based polymer actuators driven by hydraulic pressure.pdf"
)
OUT = Path(
    r"C:\Users\B3DCB~1.PER\AppData\Local\Temp\claude"
    r"\C--Users-b-pereiraazevedo-OneDrive---House-Of-HR-NV-Documents-Stage-muscle-artifici-le-mod-le-mod-le-chinois-Espace-de-travail"
    r"\c5852e10-6d38-4703-b218-a1649a5b12c9\scratchpad"
)

reader = PdfReader(str(PDF))
page = reader.pages[11]
print("n images:", len(page.images))
img = Image.open(BytesIO(page.images[0].data)).convert("RGB")
print("size:", img.size)
img.save(OUT / "figure7_extracted.png")
print("saved")
text = page.extract_text()
print("---- page text (first 3500 chars) ----")
print(text[:3500].encode("ascii", "replace").decode())
