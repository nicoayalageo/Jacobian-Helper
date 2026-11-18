# Jacobian Helper

Diagnóstico de una calibración de PEST o PEST++ estancada a partir del Jacobiano ya calculado, sin correr el modelo. Funciona en el navegador con Pyodide y no envía los archivos a ningún servidor.

Página: https://nicoayalageo.github.io/Jacobian-Helper/

## Archivos que lee

| Archivo | Obligatorio | Contenido |
|---|---|---|
| `caso.pst` | Sí | Archivo de control, versión 1 de PEST o versión 2 de PEST++. Si usa tablas externas hay que agregar también los CSV |
| `caso.jco` o `caso.jcb` | Sí | Jacobiano que escribe PEST o PEST++ |
| `caso.res` o `caso.rei` | Sí | Residuales de la última iteración |
| `caso.par` | Recomendado | Valores finales de los parámetros. Si falta se usan los iniciales del archivo de control |
| `caso.rec` | Opcional | Registro de la corrida, para informar cómo terminó |

## Qué entrega

1. Cómo terminó la corrida, según el registro: motivo, iteraciones usadas y último parámetro que alcanzó su límite.
2. El factor con que se acorta el paso de Gauss-Newton, β = mín(límite_i / |Δ_i|), y el parámetro que lo fija.
3. La mejora de la función objetivo tal como está, sosteniendo los parámetros que fijan el paso y sin el límite de cambio, con un veredicto: estrangulada, sana o nada que recuperar.
4. Qué hacer: reinicio desde el mismo Jacobiano (`pest caso /i`, `++base_jacobian`) o cambio del tipo de límite.

La mejora es una predicción lineal y debe comprobarse con una iteración antes de adoptarla.

## Código

| Ruta | Contenido |
|---|---|
| `index.html` | Página completa: estilos, textos en español e inglés y lógica |
| `py/discriminador.py` | Cálculo del diagnóstico. Es el mismo archivo del repositorio del estudio; su SHA-256 está en `py/discriminador.sha256` |
| `py/pestio.py` | Lectura de los archivos de PEST sin dependencias externas |
| `ejemplos/` | Cuatro corridas de ejemplo |
| `tools/` | Scripts que generan el fondo de la cabecera y la imagen para compartir |
| `img/` | Imágenes de la página |

## Ejemplos

- `ejemplos/detenida`: Freyberg en MODFLOW 6 con puntos piloto y decaimiento de K con la profundidad por unidad. PEST 18.25 declara la optimización completa en la iteración 2 con phi 1338. Veredicto: estrangulada.
- `ejemplos/pestpp`: el mismo caso, con el Jacobiano de la primera iteración de PEST++GLM. Veredicto: estrangulada.
- `ejemplos/sana`: Freyberg con una unidad de roca bajo el acuífero y exponentes iniciados en −1, tras la primera iteración de PEST++GLM. Veredicto: sana.
- `ejemplos/convergida`: una variante del mismo modelo que llegó a phi 596. Veredicto: nada que recuperar.

## Validación

El veredicto se comparó con lo ocurrido en 183 corridas de PEST 18.25 y PEST++ 5.2.16 terminadas por su propio criterio, de 14 variantes del modelo de Freyberg, y coincidió en 182. Las variantes comparten el modelo base, por lo que el resultado no reemplaza la comprobación en un modelo propio.

## Referencia

Ayala, N. y López, D. (2026). Antes de recalibrar, cómo saber si una corrida PEST estancada se puede rescatar desde el Jacobiano ya calculado. XVII Congreso Geológico Chileno, Valdivia, sesión AT2.6.

Doherty, J. (2018). PEST, Model-Independent Parameter Estimation. User Manual Part I. Watermark Numerical Computing.

White, J. T., Hunt, R. J., Fienen, M. N. y Doherty, J. E. (2020). Approaches to highly parameterized inversion: PEST++ Version 5. U.S. Geological Survey Techniques and Methods 7-C26.

## Herramientas

La revisión de los manuales de PEST y PEST++ se apoyó en [MODFLOW AI](https://agent.modflow.ai) (GWLab), un servidor de consulta de la documentación de MODFLOW, PEST y PEST++. No interviene en el cálculo del diagnóstico.

Licencia MIT.
