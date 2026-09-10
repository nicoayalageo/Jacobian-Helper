"""Lectores mínimos de archivos PEST / PEST++ sin pyemu, para Jacobian-Helper (corre en Pyodide).

    leer_pst(texto, externos={})  -> dict(control=..., par=[...], obs=[...])   (.pst v1 y v2, con CSV externos)
    leer_jco(datos_bytes)         -> (J ndarray, filas, columnas)               (formato binario PEST, port de pyemu.Matrix.read_binary)
    leer_res(texto)               -> dict nombre -> dict(residual, weight, measured, modelled)   (.res de PEST, .rei de PEST++)
    leer_par(texto)               -> dict nombre -> valor
    leer_rec(texto)               -> dict(razon, nopt, noptmax, topo_limite)     (PEST y PEST++)
    preparar(pst, jco, res, par)  -> lo que necesita discriminador.discriminar_grupo
    diagnosticar(archivos)        -> resultado completo (dict) a partir de {nombre: bytes|str}

Todo se valida contra pyemu en test_pestio.py sobre las corridas de este repo.
"""
import io, re, struct
import numpy as np

# ----------------------------------------------------------------------------- .pst
_CTRL_V1 = {  # (línea dentro de '* control data', índice del token): clave  — PEST man. I §4.2
    (4, 0): "relparmax", (4, 1): "facparmax", (4, 2): "facorig",
    (6, 0): "noptmax", (6, 1): "phiredstp", (6, 2): "nphistp", (6, 3): "nphinored", (6, 4): "relparstp", (6, 5): "nrelpar",
}
_PAR_COLS = ["parnme", "partrans", "parchglim", "parval1", "parlbnd", "parubnd", "pargp", "scale", "offset", "dercom"]
_OBS_COLS = ["obsnme", "obsval", "weight", "obgnme"]
_NUM = {"parval1", "parlbnd", "parubnd", "scale", "offset", "obsval", "weight"}


def _secciones(lineas):
    """Devuelve lista de (titulo, [lineas]) para cada bloque '* ...'."""
    out, tit, buf = [], None, []
    for l in lineas:
        s = l.strip()
        if not s: continue
        if s.startswith("*"):
            if tit is not None: out.append((tit, buf))
            tit, buf = s[1:].strip().lower(), []
        elif tit is not None:
            buf.append(l.rstrip("\n"))
    if tit is not None: out.append((tit, buf))
    return out


def _csv(texto):
    """CSV mínimo (encabezado + filas) -> lista de dict, claves en minúscula."""
    filas = [l for l in texto.splitlines() if l.strip()]
    hdr = [h.strip().lower() for h in filas[0].split(",")]
    out = []
    for l in filas[1:]:
        v = [x.strip() for x in l.split(",")]
        out.append(dict(zip(hdr, v)))
    return out


def _tabla(lineas, cols, externos, externa):
    """Lee una tabla de parámetros u observaciones, en línea (v1 o v2 con encabezado) o externa (CSV)."""
    filas = []
    if externa:
        for l in lineas:
            nombre = l.split()[0]
            clave = next((k for k in externos if k.lower() == nombre.lower() or k.lower().endswith("/" + nombre.lower())), None)
            if clave is None:
                raise FileNotFoundError(f"el .pst pide el archivo externo '{nombre}' y no se entregó")
            filas += _csv(externos[clave])
    else:
        hdr = None
        for l in lineas:
            t = l.split()
            if not t: continue
            if hdr is None and t[0].lower() == cols[0]:      # v2 en línea: primera fila es el encabezado
                hdr = [x.lower() for x in t]; continue
            if hdr is not None:
                filas.append(dict(zip(hdr, t)))
            elif len(t) >= len(cols):                        # v1 (las líneas de 'tied' tienen 2 tokens y se omiten)
                filas.append(dict(zip(cols, t)))
    out = []
    for f in filas:
        d = {}
        for c in cols:
            v = f.get(c)
            if v is None: continue
            d[c] = float(v) if c in _NUM else v.strip().lower()
        out.append(d)
    return out


def leer_pst(texto, externos=None):
    externos = externos or {}
    lineas = texto.splitlines()
    ctrl, par, obs = {}, [], []
    v2 = "version=2" in lineas[0].lower()
    for tit, buf in _secciones(lineas):
        if tit.startswith("control data keyword"):
            for l in buf:
                t = l.split()
                if len(t) >= 2:
                    k = t[0].lower()
                    try: ctrl[k] = float(t[1])
                    except ValueError: ctrl[k] = t[1].lower()
        elif tit.startswith("control data"):
            for (i, j), k in _CTRL_V1.items():
                if i < len(buf):
                    t = buf[i].split()
                    if j < len(t):
                        try: ctrl[k] = float(t[j])
                        except ValueError: pass
            if len(buf) > 4:                                  # ABSPARMAX(N)=valor en la línea 5
                for m in re.finditer(r"absparmax\((\d+)\)\s*=\s*([0-9.eE+-]+)", buf[4], re.I):
                    ctrl.setdefault("absparmax", {})[int(m.group(1))] = float(m.group(2))
        elif tit.startswith("parameter data"):
            par = _tabla(buf, _PAR_COLS, externos, "external" in tit)
        elif tit.startswith("observation data"):
            obs = _tabla(buf, _OBS_COLS, externos, "external" in tit)
    for k in ("relparmax", "facparmax", "facorig", "noptmax", "phiredstp", "nphistp"):
        ctrl.setdefault(k, {"relparmax": 3.0, "facparmax": 3.0, "facorig": 0.001, "noptmax": 30, "phiredstp": 0.005, "nphistp": 4}[k])
    return dict(control=ctrl, par=par, obs=obs, version=2 if v2 else 1)


# ----------------------------------------------------------------------------- .jco / .jcb
def leer_jco(datos):
    """Port de pyemu.Matrix.read_binary (formatos nuevo y antiguo de PEST)."""
    f = io.BytesIO(datos)
    itemp1, itemp2, icount = np.frombuffer(f.read(12), dtype=[("a", "<i4"), ("b", "<i4"), ("c", "<i4")])[0]
    if itemp1 == 0 and itemp2 == icount:
        raise ValueError("formato 'dense' de PEST++ no soportado todavía")
    ncol, nrow = abs(int(itemp1)), abs(int(itemp2))
    x = np.zeros((nrow, ncol))
    if itemp1 >= 0:                                           # formato nuevo: (i, j, valor) en base 0
        rec = np.frombuffer(f.read(16 * int(icount)), dtype=[("i", "<i4"), ("j", "<i4"), ("d", "<f8")])
        x[rec["i"], rec["j"]] = rec["d"]; lp, lo = 200, 200
    else:                                                     # formato antiguo: índice lineal en base 1
        rec = np.frombuffer(f.read(12 * int(icount)), dtype=[("j", "<i4"), ("d", "<f8")])
        ic = ((rec["j"] - 1) // nrow) + 1; ir = rec["j"] - ((ic - 1) * nrow)
        x[ir - 1, ic - 1] = rec["d"]; lp, lo = 12, 20
    cols = [f.read(lp).decode(errors="ignore").strip().lower() for _ in range(ncol)]
    filas = [f.read(lo).decode(errors="ignore").strip().lower() for _ in range(nrow)]
    return x, filas, cols


# ----------------------------------------------------------------------------- .res / .rei / .par / .rec
def leer_res(texto):
    out, hdr = {}, None
    for l in texto.splitlines():
        t = l.split()
        if not t: continue
        if hdr is None:
            if t[0].lower() == "name": hdr = [x.lower() for x in t]
            continue
        d = dict(zip(hdr, t))
        try:
            out[d["name"].lower()] = dict(residual=float(d["residual"]), weight=float(d["weight"]),
                                          measured=float(d["measured"]), modelled=float(d["modelled"]))
        except (KeyError, ValueError):
            continue
    return out


def leer_par(texto):
    out = {}
    for l in texto.splitlines()[1:]:
        t = l.split()
        if len(t) >= 2:
            try: out[t[0].lower()] = float(t[1])
            except ValueError: pass
    return out


def leer_rec(t):
    info = {}
    m = re.search(r"Reason for terminating[^\n]*:\s*(.+)", t)
    if m: info["razon"] = m.group(1).strip(); info["codigo"] = "PEST++"
    m = re.search(r"NOPTMAX\s*=\s*(\d+)\s*;\s*NOPT at termination\s*=\s*(\d+)", t)
    if m: info["noptmax"], info["nopt"] = int(m.group(1)), int(m.group(2))
    if "razon" not in info:
        m = re.search(r"(Optimisation complete:[^\n]*(?:\n[^\n]*)?)", t)
        if m: info["razon"] = " ".join(m.group(1).split()); info["codigo"] = "PEST"
        m = re.search(r"Maximum number of optimisation iterations\s*:\s*(\d+)", t)
        its = re.findall(r"OPTIMISATION ITERATION NO\.\s*:\s*(\d+)", t)
        if m and its: info["noptmax"], info["nopt"] = int(m.group(1)), int(its[-1])
    phis = re.findall(r"Lowest phi this iteration:\s*([0-9.eE+-]+)", t) or re.findall(r"^\s*(\d+)\s+[0-9.eE+-]+\s+([0-9.eE+-]+)", t, re.M)
    if phis: info["phis"] = [float(p if isinstance(p, str) else p[-1]) for p in phis]
    lim = re.findall(r"Maximum\s+(relative|factor|absolute)\s+change[^\n]*?\[\"?([^\]\"]+)\"?\]", t)
    if lim: info["topo_limite"] = lim[-1][1].strip()
    return info


# ----------------------------------------------------------------------------- armar el problema
def preparar(pst, jco, res, par=None, relparmax=None, facparmax=None, absparmax=None):
    J, filas, cols = jco
    pdat = {p["parnme"]: p for p in pst["par"]}
    odat = {o["obsnme"]: o for o in pst["obs"]}
    aj = [c for c in cols if c in pdat and pdat[c].get("partrans", "none") not in ("fixed", "tied")]
    on = [f for f in filas if f in res and f in odat]
    if not aj: raise ValueError("ninguna columna del Jacobiano corresponde a un parámetro ajustable del .pst")
    if not on: raise ValueError("ninguna fila del Jacobiano tiene residual en el .res/.rei")
    fi = {f: k for k, f in enumerate(filas)}; ci = {c: k for k, c in enumerate(cols)}
    Jm = J[np.ix_([fi[f] for f in on], [ci[c] for c in aj])]
    r = np.array([res[f]["residual"] for f in on]); w = np.array([odat[f]["weight"] for f in on], float)
    vals = np.array([(par or {}).get(c, pdat[c]["parval1"]) for c in aj], float)
    tipos = ["absolute" if str(pdat[c].get("parchglim", "relative")).startswith("abs") else str(pdat[c].get("parchglim", "relative")) for c in aj]
    ctrl = pst["control"]
    rel = relparmax if relparmax else float(ctrl.get("relparmax", 3.0))
    fac = facparmax if facparmax else float(ctrl.get("facparmax", 3.0))
    ab = absparmax if absparmax else float(ctrl.get("absparmax", {1: 0.3}).get(1, 0.3)) if isinstance(ctrl.get("absparmax"), dict) else (absparmax or 0.3)
    return dict(J=Jm, r=r, w=w, par_names=aj, par_values=vals, tipos=tipos, rel=rel, fac=fac, abs_limit=ab,
                n_obs=len(on), n_par=len(aj), n_par_total=len(cols), usa_par_final=par is not None)


def diagnosticar(archivos, relparmax=None, facparmax=None):
    """archivos: {nombre: bytes|str}. Detecta por extensión: .pst, .jco/.jcb, .res/.rei, .par, .rec, .csv (externos)."""
    from discriminador import discriminar_grupo
    txt = {k: (v.decode("utf-8", errors="ignore") if isinstance(v, (bytes, bytearray)) else v) for k, v in archivos.items()}
    def uno(exts):
        c = [k for k in archivos if k.lower().rsplit(".", 1)[-1] in exts]
        return c[0] if c else None
    k_pst, k_jco, k_res, k_par, k_rec = uno(("pst",)), uno(("jco", "jcb")), None, uno(("par",)), uno(("rec",))
    res_c = sorted([k for k in archivos if k.lower().rsplit(".", 1)[-1] in ("res", "rei")], key=lambda k: (k.lower().endswith(".res"), k))
    k_res = res_c[-1] if res_c else None                      # si hay .res y .rei, prefiere el .res (final)
    faltan = [n for n, k in (("pst", k_pst), ("jco/jcb", k_jco), ("res/rei", k_res)) if k is None]
    if faltan: raise FileNotFoundError("faltan archivos: " + ", ".join(faltan))
    externos = {k: txt[k] for k in txt if k.lower().endswith(".csv")}
    pst = leer_pst(txt[k_pst], externos)
    jco = leer_jco(archivos[k_jco] if isinstance(archivos[k_jco], (bytes, bytearray)) else archivos[k_jco].encode("latin-1"))
    res = leer_res(txt[k_res]); par = leer_par(txt[k_par]) if k_par else None
    P = preparar(pst, jco, res, par, relparmax, facparmax)
    d = discriminar_grupo(jco=P["J"], residuals=P["r"], weights=P["w"], par_names=P["par_names"], par_values=P["par_values"],
                          tipos=P["tipos"], rel_limit=P["rel"], factor_limit=P["fac"], abs_limit=P["abs_limit"])
    d = {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in d.items()}
    d.update(rec=leer_rec(txt[k_rec]) if k_rec else None, n_obs=P["n_obs"], n_par=P["n_par"], n_par_total=P["n_par_total"],
             relparmax=P["rel"], facparmax=P["fac"], usa_par_final=P["usa_par_final"],
             valores={n: float(v) for n, v in zip(P["par_names"], P["par_values"])}, tipos=dict(zip(P["par_names"], P["tipos"])),
             control=pst["control"], archivos=dict(pst=k_pst, jco=k_jco, res=k_res, par=k_par, rec=k_rec))
    return d
