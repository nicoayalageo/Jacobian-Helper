# Jacobian Helper

**¿PEST se detuvo antes de tiempo?** Diagnóstico de una calibración estancada a partir del Jacobiano ya calculado, sin correr el modelo. Corre entero en el navegador (Pyodide): los archivos no salen de tu computador.

**Página:** https://nicoayalageo.github.io/Jacobian-Helper/

## Qué necesita

| archivo | obligatorio | de dónde |
|---|---|---|
| `caso.pst` | sí | el archivo de control (v1 de PEST/PEST_HP o v2 de PEST++; si usa tablas externas, arrastra también los CSV) |
| `caso.jco` / `caso.jcb` | sí | el Jacobiano final que escribe PEST / PEST++ |
| `caso.res` / `caso.rei` | sí | residuales de la última iteración |
| `caso.par` | recomendado | valores finales de los parámetros (si falta, se usan los iniciales del .pst y se avisa) |
| `caso.rec` | opcional | para leer cómo terminó y quién topó el límite |

## Qué responde

1. **Cómo terminó** (del `.rec`): razón, iteraciones usadas de NOPTMAX, qué parámetro topó el límite.
2. **Qué dice el Jacobiano**: el factor de acortamiento β = mín(límite_i/|Δ_i|) del paso de Gauss-Newton y quién lo fija.
3. **Conviene intervenir**: ganancia de phi tal como está, sosteniendo al grupo de parámetros que fijan β, y sin el limitador. Veredicto `STRANGLED` / `HEALTHY` / `NO_INFORMATION`.
4. **Qué hacer**: reinicio sin recalcular el Jacobiano (`pest caso /i`, `++base_jacobian`), cambio del tipo de límite (OFFSET+log, `absolute(N)`), o revisar el parámetro.

## Método

`py/discriminador.py` es el código del estudio, sin modificar (hash en `py/discriminador.sha256`). `py/pestio.py` lee los archivos de PEST sin dependencias; se validó contra pyemu en 188 corridas (`test_pestio.py` en el repositorio del estudio). La ganancia es una predicción lineal: valida siempre con una corrida.

Ayala, N. (Arcadis) & López, D. (GWLab), 2026. *Antes de recalibrar: cómo saber si una corrida PEST estancada se puede rescatar desde el Jacobiano ya calculado.* Referencias: Doherty (2018) PEST Manual Part I; White et al. (2020) PEST++ v5, USGS T&M 7-C26.

## Ejemplos incluidos

- `ejemplos/congelado/`: Freyberg MODFLOW 6 + puntos piloto + decaimiento de K con la profundidad por unidad; PEST 18.25 declara "Optimisation complete" en la iteración 2 con phi 1338 (óptimo 629). Veredicto `STRANGLED`.
- `ejemplos/sano/`: el mismo modelo con 44 puntos piloto, convergido a phi 597. Veredicto `NO_INFORMATION`.

Licencia MIT.
