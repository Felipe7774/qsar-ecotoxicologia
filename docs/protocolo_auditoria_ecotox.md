# Protocolo de auditoría de calidad ECOTOX

## Objetivo

Verificar que una muestra representativa del conjunto QSAR conserve correctamente los criterios de inclusión y la trazabilidad desde ECOTOX.

## Muestra

`src/create_ecotox_audit_sample.py` selecciona 30 sustancias de forma reproducible, con diez sustancias de cada estrato de toxicidad basada en `log10(EC50 mg/L)`.

## Procedimiento por sustancia

1. Buscar el CASRN o el `result_id` en ECOTOX y confirmar la identidad química.
2. Confirmar que la especie sea *Daphnia magna*.
3. Confirmar agua dulce y exposición de 48 h.
4. Confirmar endpoint EC50 y medición de inmovilización.
5. Comparar valor y unidad de concentración con el registro exportado.
6. Registrar `Sí`, `No` o `No verificable` en cada campo de revisión.
7. En `decision`, registrar `Aceptar`, `Excluir` o `Revisar` y justificar la decisión en `notas`.

## Regla de decisión

Un registro se acepta cuando todos los criterios de inclusión se verifican. Un error de especie, duración, endpoint, medición o unidad implica excluir ese registro y regenerar el conjunto antes de reentrenar. Si la fuente no permite verificarlo, se marca como `Revisar` y no se usa para justificar la calidad del conjunto.

## Evidencia

Conservar el CSV de la auditoría diligenciado, la fecha, el revisor y una copia o enlace estable de la consulta ECOTOX cuando sea posible.

## Auditoría automática implementada

`src/run_ecotox_audit.py` consulta directamente las tablas `results.txt` y `tests.txt` de la descarga ASCII de ECOTOX. Para cada una de las 30 sustancias verifica CASRN, especie, agua dulce, 48 h, EC50, inmovilización, operador, valor positivo y unidad permitida. Esta verificación comprueba la trazabilidad del procesamiento, pero no sustituye la consulta de la publicación científica primaria.
