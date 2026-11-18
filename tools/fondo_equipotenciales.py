#!/usr/bin/env python
"""Genera img/equipotenciales.svg, el fondo de la cabecera de la página.

Son las curvas equipotenciales de la capa superior del modelo de Freyberg (MODFLOW 6) usado en el
estudio, en el último paso de tiempo. El modelo se dibuja apaisado (el norte queda a la izquierda).

    python tools/fondo_equipotenciales.py ruta/al/freyberg6_freyberg.hds ruta/al/activo.npy
"""
import sys
from pathlib import Path

import flopy
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import distance_transform_edt, gaussian_filter, zoom

hds, activo = sys.argv[1], sys.argv[2]
SALIDA = Path(__file__).resolve().parent.parent / "img" / "equipotenciales.svg"
ANCHO, ALTO, FACTOR, NIVELES = 1000.0, 500.0, 8, 22

h = flopy.utils.HeadFile(hds).get_data()[0].astype(float)          # capa superior, último tiempo
act = np.load(activo)
h[~act] = np.nan
# las celdas inactivas se rellenan con el valor activo más cercano para que las curvas no se corten
idx = distance_transform_edt(np.isnan(h), return_distances=False, return_indices=True)
h = gaussian_filter(zoom(h[tuple(idx)], FACTOR, order=3), FACTOR * 0.6)
h = h.T[::-1]                                                        # apaisado
ny, nx = h.shape
cs = plt.contour(np.linspace(0, ANCHO, nx), np.linspace(0, ALTO, ny), h,
                 levels=np.linspace(np.nanpercentile(h, 1), np.nanpercentile(h, 99), NIVELES))
trazos = []
for i, nivel in enumerate(cs.allsegs):
    for seg in nivel:
        if len(seg) < 8:
            continue
        seg = seg[::3]
        d = "M" + " L".join(f"{x:.1f} {y:.1f}" for x, y in seg)
        trazos.append(f'<path d="{d}" stroke-width="{1.6 if i % 5 == 0 else 0.8}"/>')
SALIDA.write_text(
    f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {ANCHO:.0f} {ALTO:.0f}" preserveAspectRatio="xMidYMid slice" '
    f'fill="none" stroke="#1d4f8f" stroke-linecap="round" stroke-linejoin="round">\n' + "\n".join(trazos) + "\n</svg>\n")
print("escrito", SALIDA, f"{SALIDA.stat().st_size / 1024:.0f} kB,", len(trazos), "trazos")
