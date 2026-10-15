# Guía de uso local

## Abrir la herramienta

En PowerShell, ubicado en la carpeta `qsar_ecotoxicologia`, ejecutar:

```powershell
.\.venv\Scripts\streamlit.exe run .\app.py
```

Abrir `http://localhost:8501` en el navegador.

## Ingresar una sustancia

1. Escribir nombre y CAS si se conocen.
2. Pegar un SMILES válido. Este campo es obligatorio.
3. Agregar logKow medido cuando esté disponible.
4. Agregar concentración observada en µg/L y remoción en porcentaje cuando existan fuentes documentadas.
5. Seleccionar **Predecir y priorizar**.

## Interpretar la respuesta

- **EC50 estimada:** toxicidad aguda predicha para el endpoint definido.
- **Prioridad:** Baja, Media o Alta según los componentes disponibles.
- **Completitud ambiental:** porcentaje de los cuatro componentes incluidos.
- **Advertencia de aplicabilidad:** aparece si la sustancia se aleja del espacio químico de entrenamiento.

Una clasificación preliminar requiere datos de ocurrencia y remoción antes de usarse para comparar sustancias en un contexto ambiental.
