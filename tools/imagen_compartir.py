#!/usr/bin/env python
"""Genera img/compartir.png (1200 x 630), la imagen que muestran los enlaces a la página.

Usa las mismas equipotenciales del fondo (img/equipotenciales.svg) y el título de la página.

    python tools/imagen_compartir.py
"""
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

RAIZ = Path(__file__).resolve().parent.parent
svg = (RAIZ / "img" / "equipotenciales.svg").read_text()
fig = plt.figure(figsize=(12, 6.3), dpi=100, facecolor="#f6f3ec")
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 1000); ax.set_ylim(525, 0); ax.axis("off")
for d, ancho in re.findall(r'<path d="M([^"]+)" stroke-width="([0-9.]+)"', svg):
    pts = [tuple(map(float, p.split())) for p in d.split(" L")]
    ax.plot([p[0] for p in pts], [p[1] for p in pts], color="#1d4f8f", lw=float(ancho) * 1.3, alpha=0.2, solid_capstyle="round")
fig.text(0.07, 0.60, "¿PEST se detuvo\nantes de tiempo?", fontsize=62, fontweight="bold", color="#152232",
         family=["Archivo", "Helvetica Neue", "Arial"], linespacing=1.05, va="center")
fig.text(0.07, 0.24, "Jacobian Helper: diagnóstico de una calibración estancada\na partir del Jacobiano ya calculado",
         fontsize=23, color="#4a5868", family=["Public Sans", "Helvetica Neue", "Arial"], linespacing=1.4, va="center")
fig.savefig(RAIZ / "img" / "compartir.png", dpi=100, facecolor=fig.get_facecolor())
print("escrito", RAIZ / "img" / "compartir.png")
