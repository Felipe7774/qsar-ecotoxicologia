# Modelo QSAR de toxicidad acuática

Proyecto de tesis para entrenar y validar un modelo supervisado que estime un endpoint ecotoxicológico medido experimentalmente.

## Estado actual

- Endpoint QSAR definido: EC50 de inmovilización a 48 h en *Daphnia magna* de agua dulce.
- `src/filter_training_data.py` audita el consolidado y genera los candidatos reproducibles para ese endpoint.
- El entrenamiento sólo se habilita cuando haya al menos 300 sustancias únicas con el endpoint completo y estructura química verificable.

## Preparar datos desde el consolidado

En la terminal de VS Code, desde esta carpeta, ejecuta:

```powershell
.\.venv\Scripts\python.exe .\src\filter_training_data.py
```

El script no completa ni supone campos faltantes. Guarda los ensayos depurados en `data/processed/` y el conteo de elegibilidad en `results/filtering_report.json`.

## Resolver identidades químicas de los candidatos

Después del filtrado, ejecuta:

```powershell
.\.venv\Scripts\python.exe .\src\resolve_candidate_structures.py
```

El proceso consulta PubChem con CAS o nombre y crea `data/processed/qsar_manual_review.csv`. La estructura obtenida se marca como pendiente hasta que se confirme contra ECOTOX; allí se registra también la duración y la fuente que la respalda.

## Entrenar el QSAR seleccionado

Cuando `qsar_training_ready.csv` tenga al menos 300 sustancias únicas, EC50 a 48 h verificable y una columna `smiles`, ejecuta:

```powershell
.\.venv\Scripts\python.exe .\src\train_hybrid_qsar.py
```

El modelo seleccionado combina descriptores 2D y fingerprints Morgan en Extra Trees, con validación cruzada agrupada por familia estructural. Guarda el modelo, los resultados por partición y las métricas promedio en `models/` y `results/`.

Para predecir una sustancia individual después de entrenar:

```powershell
.\.venv\Scripts\python.exe .\src\predict_qsar.py "CC(=O)NC1=CC=C(O)C=C1"
```

La salida incluye EC50 estimada en mg/L, su valor logarítmico y una advertencia si alguno de los descriptores cae fuera del intervalo observado durante el entrenamiento.

## Priorizar sustancias para el contexto ambiental

La prioridad combina la EC50 predicha con logKow, concentración observada y remoción cuando esos datos existen. La metodología y los umbrales están en `docs/prioritization_methodology.md`.

```powershell
.\.venv\Scripts\python.exe .\src\prioritize_environmental_risk.py --input .\data\processed\priority_input.csv
```

Para construir esa entrada con los candidatos ya presentes en el consolidado:

```powershell
.\.venv\Scripts\python.exe .\src\build_priority_input.py
```

## Usar la aplicación local

Instala las dependencias y ejecuta:

```powershell
.\.venv\Scripts\python.exe -m pip install -r .\requirements.txt
.\.venv\Scripts\streamlit.exe run .\app.py
```

La aplicación abre una pestaña local donde puedes ingresar SMILES y, de forma opcional, logKow, concentración observada y remoción. También acepta un CSV de hasta 1.000 sustancias y permite descargar el resultado.

## Publicar la aplicación

La configuración reproducible para alojarla con Python y RDKit está en
`Dockerfile` y `render.yaml`. Antes de publicar ejecuta:

```powershell
.\.venv\Scripts\python.exe .\src\preflight_public_release.py
docker build -t qsar-ecotoxicologia .
```

La guía completa, incluido qué datos no deben incluirse en el repositorio, está
en `DEPLOYMENT.md`.

## Generar figuras para la tesis

Después de entrenar, genera las figuras de validación y explicabilidad con:

```powershell
.\.venv\Scripts\python.exe .\src\create_qsar_figures.py
```

## Bitácora

El estado del proyecto, evidencia disponible y tareas pendientes se mantienen en `docs/bitacora_desarrollo.md`.

El informe técnico y la guía de uso se encuentran en `docs/informe_tecnico_qsar.md` y `docs/guia_usuario_local.md`.

## Auditoría de calidad ECOTOX

Para crear la muestra de revisión de 30 sustancias:

```powershell
.\.venv\Scripts\python.exe .\src\create_ecotox_audit_sample.py
```

El protocolo de revisión está en `docs/protocolo_auditoria_ecotox.md`.

Para ejecutar la auditoría automática contra la descarga original:

```powershell
.\.venv\Scripts\python.exe .\src\run_ecotox_audit.py
```

## Comparar fingerprints moleculares

Para comparar descriptores 2D, fingerprints Morgan y un modelo híbrido usando la misma validación por andamio:

```powershell
.\.venv\Scripts\python.exe .\src\benchmark_fingerprint_models.py
```

## Extraer el endpoint desde ECOTOX completo

Con el archivo ASCII oficial en `data/raw/ecotox_ascii_09_15_2026.zip`, ejecuta:

```powershell
.\.venv\Scripts\python.exe .\src\extract_ecotox_endpoint.py
```

El script reconstruye el vínculo entre pruebas, resultados y catálogo químico. Primero genera un archivo de auditoría con todos los EC50 de *Daphnia magna*, 48 h y agua dulce; el código de efecto y las unidades se revisan antes de convertirlo en el conjunto de entrenamiento.

Para aplicar ese control de calidad y consolidar la mediana de los ensayos por CASRN:

```powershell
.\.venv\Scripts\python.exe .\src\build_qsar_training_table.py
```

El resultado `qsar_training_before_structures.csv` contiene una fila por sustancia. Sólo faltará obtener y validar los SMILES para entrenar.

Para resolver los SMILES por CASRN y validarlos con RDKit:

```powershell
.\.venv\Scripts\python.exe .\src\resolve_training_structures.py
```

Este paso guarda `qsar_training_with_structures.csv`; las sustancias sin SMILES interpretable se documentan y no entran al modelo.

Para publicar el conjunto definitivo antes del entrenamiento:

```powershell
.\.venv\Scripts\python.exe .\src\finalize_training_dataset.py
```

## Estructura

- `data/raw/`: descargas originales sin modificar.
- `data/processed/`: conjuntos depurados y preparados para el modelo.
- `notebooks/`: exploración y visualizaciones reproducibles.
- `src/`: scripts de limpieza, creación de descriptores, entrenamiento y evaluación.
- `models/`: modelos entrenados.
- `results/`: métricas, tablas y figuras finales.
- `docs/`: metodología, diccionario de datos y fuentes.
