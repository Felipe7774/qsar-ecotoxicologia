# Descarga reproducible de ECOTOX

El conjunto de la tesis necesita registros comparables para el endpoint QSAR, no todos los registros de contaminantes emergentes.

## Filtros que se deben aplicar en ECOTOX

1. Hábitat: acuático.
2. Especie: *Daphnia magna*.
3. Endpoint: `EC50`.
4. Efecto: inmovilización.
5. Duración: `48 h`.
6. Medio: agua dulce.
7. Resultado: concentración numérica con unidad.

## Campos que se deben exportar

Activa la salida completa, no la versión condensada. Conserva como mínimo: nombre químico, CASRN, SMILES o identificador químico, especie, endpoint, efecto, duración de exposición, medio, concentración, unidad, referencia y año de publicación.

Guarda el resultado descargado en `data/raw/` como `ecotox_export.csv` o `ecotox_export.xlsx` y ejecuta:

```powershell
.\.venv\Scripts\python.exe .\src\prepare_ecotox.py
```

La selección y exportación se hacen desde las funciones Search o Explore de ECOTOX. La documentación de EPA explica que esos filtros pueden combinarse y que el resultado se puede exportar como CSV o Excel.
