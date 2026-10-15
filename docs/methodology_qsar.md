# Metodología QSAR implementada

## Objetivo

Predecir `log10(EC50_mg_L)` de sustancias individuales para *Daphnia magna*, agua dulce, 48 h e inmovilización. La categoría de riesgo no se usa como etiqueta de entrenamiento.

## Origen y depuración de datos

La fuente es la descarga ASCII de ECOTOX del 15 de septiembre de 2026. Se reconstruyeron las tablas de pruebas, resultados, especies y sustancias mediante `test_id` y `test_cas`.

Se conservaron únicamente los registros con:

1. Especie *Daphnia magna*.
2. Medio `FW` (agua dulce).
3. Exposición de 48 h.
4. Endpoint `EC50`.
5. Medición `IMBL` o `IMBL/` (inmovilización).
6. Valor de concentración numérico, positivo y sin operador de desigualdad.
7. Unidad convertible a mg/L.

Las unidades mg/L, µg/L, ng/L, ppm, ppb y equivalentes de ingrediente activo se normalizaron a mg/L. Los valores se transformaron con `log10` y se consolidaron por CASRN mediante la mediana. Las unidades molares, los intervalos y los límites como `> x` se excluyeron.

## Estructuras y descriptores

Los CASRN se consultaron en PubChem para obtener SMILES. RDKit validó cada estructura; las estructuras no interpretables se excluyeron.

Para cada SMILES se calcularon diez descriptores 2D: peso molecular, logP, TPSA, donantes y aceptores de enlaces de hidrógeno, enlaces rotables, átomos pesados, anillos, anillos aromáticos y fracción sp3. También se generó un fingerprint Morgan de radio 2 y 2.048 bits para representar subestructuras químicas.

## Entrenamiento y evaluación

Se compararon Random Forest, Extra Trees con descriptores 2D, Extra Trees con fingerprints Morgan y Extra Trees híbrido. La evaluación usa validación cruzada de cinco particiones agrupada por andamio de Murcko; así, compuestos de una misma familia estructural no se reparten entre entrenamiento y evaluación de una misma partición.

La selección se realiza por el menor MAE promedio. El modelo Extra Trees híbrido obtuvo MAE de 0.963 log10(mg/L), RMSE de 1.237 y R² de 0.423, por lo que se ajustó con todas las sustancias elegibles.

## Dominio de aplicabilidad

Para cada estructura nueva se calcula la similitud de Tanimoto de su fingerprint Morgan (radio 2, 2.048 bits) frente a las sustancias de entrenamiento. Se reporta la mayor similitud y se advierte cuando está por debajo del percentil 5 de las similitudes entre vecinos más cercanos del conjunto de entrenamiento. Esta advertencia se complementa con la revisión de los rangos de los descriptores.

## Interpretación

Las importancias de variables del modelo describen cuánto utiliza cada descriptor para reducir la variabilidad en sus árboles. No establecen causalidad toxicológica. La aplicabilidad del modelo queda limitada a sustancias químicamente similares a las estructuras incluidas en el entrenamiento.
