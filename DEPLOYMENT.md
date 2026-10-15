# Publicar la aplicación QSAR

La aplicación pública debe ejecutar Python, RDKit y el archivo
`models/qsar_daphnia_ec50.joblib`. Por esa razón se publica como aplicación
Streamlit empaquetada en Docker, no como una página estática.

## Comprobación previa

Desde la carpeta del proyecto:

```powershell
.\.venv\Scripts\python.exe .\src\preflight_public_release.py
docker build -t qsar-ecotoxicologia .
docker run --rm -p 8501:8501 qsar-ecotoxicologia
```

Abre `http://localhost:8501` y verifica una predicción individual y una carga
CSV. El informe de la comprobación queda en `results/public_release_check.json`.

## Render

1. Crea un repositorio GitHub para esta carpeta e incluye el modelo
   `models/qsar_daphnia_ec50.joblib`.
2. En Render, crea un servicio desde el repositorio y selecciona el Blueprint
   `render.yaml`.
3. El servicio toma el puerto definido por Render y publica una URL HTTPS.
4. Revisa la URL generada con una predicción individual y una carga CSV antes
   de compartirla.

## Streamlit Community Cloud

Como alternativa, conecta el mismo repositorio y selecciona `app.py` como
punto de entrada. La plataforma instala `requirements.txt`. Mantén Python 3.13
para que las versiones del modelo serializado coincidan con las usadas en el
entrenamiento.

## Qué debe permanecer privado

No subas `data/raw/`, `data/processed/`, el consolidado con información de
terceros ni credenciales. El `.gitignore` ya excluye esos datos. El único
artefacto de modelo permitido en el repositorio público es
`models/qsar_daphnia_ec50.joblib`, necesario para que la aplicación funcione.

## Uso responsable

La interfaz muestra el endpoint, el dominio de aplicabilidad y las
limitaciones. Toda publicación debe conservar esos mensajes: el resultado es
apoyo a priorización y no sustituye datos experimentales ni una decisión
regulatoria.
