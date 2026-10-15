# Paso 4: adquisición y filtrado de ECOTOX

La fuente primaria será la base ECOTOX de la EPA. ECOTOX permite exportar los resultados de una búsqueda en CSV o Excel y seleccionar los campos de especie, efecto, endpoint, duración, concentración y unidad. Para este proyecto se conservarán únicamente registros de *Daphnia magna*, EC50, 48 horas e inmovilización.

## Preparar el archivo

1. En ECOTOX realiza una búsqueda con especie `Daphnia magna`, endpoint `EC50` y efecto `immobilization`.
2. Exporta todos los campos disponibles en CSV o Excel.
3. Guarda el archivo como `data/raw/ecotox_export.csv` dentro de este proyecto.
4. Ejecuta:

```powershell
.\.venv\Scripts\python.exe src\prepare_ecotox.py
```

El programa informará las columnas detectadas, aplicará las conversiones de concentración y guardará los ensayos válidos en `data/processed/ecotox_daphnia_48h_raw.csv`.

Fuente: EPA ECOTOX Knowledgebase Resource Hub: https://www.epa.gov/comptox-tools/ecotoxicology-ecotox-knowledgebase-resource-hub
