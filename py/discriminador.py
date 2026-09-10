"""Discriminador de convergencia: ¿PEST paró por falta de informacion o por el limitador?

La logica del triage original pregunta "¿hay un K cuyo paso LIMITADO de buena
ganancia?". Con un limite relativo diminuto beta es minusculo para todo K, asi
que ningun K califica y devuelve NO_INFORMATION aunque la mejora exista.

La pregunta correcta es cuanto se PIERDE por culpa del limitador:

    g_libre    = mejor ganancia sobre K, SIN limitador (beta = 1)
    g_limitado = mejor ganancia sobre K, CON el limitador real
    perdida    = g_libre - g_limitado

    g_libre < min_gain                    -> NO_INFORMATION  (no hay nada que ganar)
    perdida  < tol                        -> HEALTHY         (el limitador no estorba)
    en otro caso                          -> STRANGLED       (el limitador es el culpable)

Se reporta ademas el parametro que fija beta a rango completo (el "rehen").
"""
import numpy as np


def _limites(valores, tipos, rel_limit=3.0, factor_limit=5.0, abs_limit=0.3):
    lim = np.empty(len(tipos))
    for i, m in enumerate(tipos):
        if m == "factor":
            lim[i] = np.log10(factor_limit)
        elif m == "absolute":
            lim[i] = abs_limit
        elif m == "relative":
            v = abs(valores[i])
            lim[i] = rel_limit * v if v > 0.0 else np.inf
        else:
            raise ValueError(m)
    return lim


def _beta(d, lim):
    m = np.abs(d) > 0
    return 1.0 if not m.any() else float(min(1.0, np.min(lim[m] / np.abs(d[m]))))


def discriminar(jco, residuals, weights, par_names, par_values, tipos,
                rel_limit=3.0, factor_limit=5.0, abs_limit=0.3,
                min_gain=0.05, tol_perdida=0.05, rank_tol=1e-10):
    WJ = weights[:, None] * np.asarray(jco, float)
    wr = weights * residuals
    phi0 = float(wr @ wr)
    U, S, Vt = np.linalg.svd(WJ, full_matrices=False)
    rank = int(np.sum(S > S[0] * rank_tol))
    lim = _limites(par_values, tipos, rel_limit, factor_limit, abs_limit)

    def phi(d):
        e = wr - WJ @ d
        return float(e @ e)

    # g_lim   : mejor ganancia con TODOS los limites reales
    # g_libre : mejor ganancia liberando SOLO al parametro que fija beta.
    #           No se usa "sin ningun limite": eso siempre da 100% amplificando
    #           ruido al dividir entre valores singulares diminutos (ver
    #           test_sin_informacion_cuando_el_jacobiano_es_casi_nulo).
    g_libre, k_libre, g_lim, k_lim = -np.inf, 0, -np.inf, 0
    for k in range(1, rank + 1):
        d = Vt[:k].T @ ((U[:, :k].T @ wr) / S[:k])
        gc = 1 - phi(_beta(d, lim) * d) / phi0
        with np.errstate(divide="ignore", invalid="ignore"):
            rr = np.where(np.abs(d) > 0, lim / np.abs(d), np.inf)
        lim_libre = lim.copy()
        lim_libre[int(np.argmin(rr))] = np.inf      # se libera al que fija beta
        gl = 1 - phi(_beta(d, lim_libre) * d) / phi0
        if gl > g_libre: g_libre, k_libre = gl, k
        if gc > g_lim:   g_lim, k_lim = gc, k

    d_full = Vt[:rank].T @ ((U[:, :rank].T @ wr) / S[:rank])
    b_full = _beta(d_full, lim)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratios = np.where(np.abs(d_full) > 0, lim / np.abs(d_full), np.inf)
    rehen = par_names[int(np.argmin(ratios))] if len(par_names) else None

    perdida = g_libre - g_lim
    if g_libre < min_gain:
        veredicto = "NO_INFORMATION"
    elif perdida < tol_perdida:
        veredicto = "HEALTHY"
    else:
        veredicto = "STRANGLED"

    return dict(veredicto=veredicto, phi0=phi0, rango=rank,
                g_libre=g_libre, k_libre=k_libre, g_lim=g_lim, k_lim=k_lim,
                perdida=perdida, beta_full=b_full, rehen=rehen)


def discriminar_grupo(jco, residuals, weights, par_names, par_values, tipos,
                      rel_limit=3.0, factor_limit=5.0, abs_limit=0.3,
                      min_gain=0.10, tol_perdida=0.05, rank_tol=1e-10, max_sostener=None,
                      solo_relativos=True):
    """Igual que discriminar(), pero identifica el GRUPO de rehenes, no solo el primero.

    Motivo (medido el 9-sep-2026, caso `pest_caso_bordes5`): cuando varios parametros
    lineales nacen cerca de cero, sostener uno solo no libera el paso porque otro pasa a
    fijar beta. Es lo que advierte PEST man. I §6.2.2: "once a particular troublesome
    parameter has been identified and held... another insensitive parameter then dominates".

    Procedimiento: se calcula la ganancia con todos los limites activos; se identifica quien
    fija beta; se le sostiene (se saca del problema); se recalcula. Se repite hasta que la
    ganancia deja de mejorar. Devuelve la secuencia de rehenes y cuanto aporta cada uno.

    Devuelve, ademas de las claves de discriminar():
      rehenes    : lista de nombres en el orden en que conviene sostenerlos
      ganancias  : ganancia acumulada tras sostener 0, 1, 2, ... de ellos
      n_sostener : cuantos conviene sostener (donde la ganancia deja de subir de forma util)
      g_grupo      : ganancia sosteniendo esos n
      g_sin_limite : ganancia si el limitador no existiera (mejor sobre K sin restricciones).
                     Es la medida correcta de lo que el limitador esta bloqueando: en el caso
                     estandar con 44 puntos piloto predijo 32.7% y el rescate real con
                     OFFSET+log dio 28.3% en la primera iteracion, mientras que la version
                     que solo liberaba un parametro predecia 3.8% (falso negativo).

    min_gain=0.10 esta calibrado sobre 101 corridas con verdad conocida (`fam_validar.py`): una
    ganancia predicha por debajo del 10% desde un solo paso linealizado cae dentro del ruido de la
    linealizacion (con 44 puntos piloto el metodo predecia ~6% estando ya en el optimo). Los casos
    que importan muestran entre 25% y 38%. Con este umbral: 99 aciertos de 101, 0 falsos positivos.

    Con solo_relativos=True (por defecto) unicamente se consideran candidatos a sostener los
    parametros de limite relativo o absoluto (los que nacen cerca de cero); sostener una K
    log-transformada en su valor inicial es una intervencion mucho mas fuerte y se evita.
    """
    WJ = weights[:, None] * np.asarray(jco, float)
    wr = weights * residuals
    phi0 = float(wr @ wr)
    lim = _limites(par_values, tipos, rel_limit, factor_limit, abs_limit)
    n = len(par_names)
    if max_sostener is None:
        max_sostener = max(1, n - 3)          # nunca dejar menos de 3 parametros libres

    def sin_limite(libres):
        """Mejor ganancia sobre K SIN limites. La truncacion regulariza, asi que no se dispara a 100%."""
        if not libres: return 0.0
        A = WJ[:, libres]
        U, S, Vt = np.linalg.svd(A, full_matrices=False)
        rank = int(np.sum(S > S[0] * rank_tol)) if S.size else 0
        g = -np.inf
        for k in range(1, rank + 1):
            d = Vt[:k].T @ ((U[:, :k].T @ wr) / S[:k])
            e = wr - A @ d
            g = max(g, 1.0 - float(e @ e) / phi0)
        return g

    def mejor(libres):
        """(ganancia limitada, quien fija beta) usando solo las columnas `libres`."""
        if len(libres) == 0:
            return 0.0, None
        A = WJ[:, libres]
        U, S, Vt = np.linalg.svd(A, full_matrices=False)
        rank = int(np.sum(S > S[0] * rank_tol)) if S.size else 0
        l = lim[libres]
        g_best, rehen = -np.inf, None
        for k in range(1, rank + 1):
            d = Vt[:k].T @ ((U[:, :k].T @ wr) / S[:k])
            b = _beta(d, l)
            e = wr - A @ (b * d)
            g = 1.0 - float(e @ e) / phi0
            if g > g_best:
                g_best = g
        d = Vt[:rank].T @ ((U[:, :rank].T @ wr) / S[:rank])
        with np.errstate(divide="ignore", invalid="ignore"):
            rr = np.where(np.abs(d) > 0, l / np.abs(d), np.inf)
        rehen = par_names[libres[int(np.argmin(rr))]]
        return g_best, rehen

    candidatos = set(range(n)) if not solo_relativos else {i for i, t in enumerate(tipos) if t != "factor"}
    libres = list(range(n))
    rehenes, ganancias = [], []
    g, rehen = mejor(libres)
    ganancias.append(g)
    for _ in range(max_sostener):
        if rehen is None:
            break
        i = par_names.index(rehen)
        if i not in libres:
            break
        if i not in candidatos:
            # el que fija beta no es candidato a sostener (p.ej. una K log-transformada):
            # se sostiene al mejor candidato restante en vez de abandonar la busqueda
            resto = [q for q in libres if q in candidatos]
            if not resto:
                break
            mejores = [(mejor([x for x in libres if x != q])[0], q) for q in resto]
            g_q, i = max(mejores)
            rehen = par_names[i]
        libres = [j for j in libres if j != i]
        rehenes.append(rehen)
        g, rehen = mejor(libres)
        ganancias.append(g)

    # cuantos conviene sostener: hasta que sostener uno mas aporte menos de 1 punto de phi
    n_sost = 0
    for k in range(1, len(ganancias)):
        if ganancias[k] - ganancias[n_sost] > 0.01:
            n_sost = k
    g_grupo = ganancias[n_sost]
    g_sin_lim = sin_limite(list(range(n)))       # lo que se alcanzaria si el limitador no existiera
    perdida = max(g_grupo, g_sin_lim) - ganancias[0]

    base = discriminar(jco, residuals, weights, par_names, par_values, tipos,
                       rel_limit, factor_limit, abs_limit, min_gain, tol_perdida, rank_tol)
    if max(g_grupo, g_sin_lim) < min_gain:
        veredicto = "NO_INFORMATION"
    elif perdida < tol_perdida:
        veredicto = "HEALTHY"
    else:
        veredicto = "STRANGLED"
    base.update(veredicto=veredicto, rehenes=rehenes[:max(n_sost, 1)], ganancias=ganancias,
                n_sostener=n_sost, g_grupo=g_grupo, g_lim=ganancias[0], perdida=perdida,
                g_sin_limite=g_sin_lim)
    return base
