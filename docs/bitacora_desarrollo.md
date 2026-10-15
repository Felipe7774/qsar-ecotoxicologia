# Bitácora de desarrollo — QSAR y priorización ecotoxicológica

Última actualización: 15 de octubre de 2026.

## Propósito del proyecto

Desarrollar una herramienta local que estime `log10(EC50 mg/L)` para inmovilización de *Daphnia magna* en agua dulce a 48 h. La predicción se integra después con una priorización ambiental explicable.

## Estado general

| Área | Estado | Evidencia principal |
|---|---|---|
| Entorno local | Completado | Entorno virtual, dependencias y estructura del proyecto |
| Datos ECOTOX | Completado | 1.449 ensayos controlados y 675 sustancias antes de estructuras |
| Auditoría de trazabilidad | Completado | 30 de 30 registros aceptados automáticamente contra tablas ECOTOX originales |
| Estructuras químicas | Completado | 408 SMILES validados con RDKit |
| QSAR inicial | Completado | Extra Trees entrenado y guardado |
| Validación | Completado | Validación cruzada agrupada por andamio estructural |
| Dominio de aplicabilidad | Completado | Rangos de descriptores y similitud de Tanimoto |
| Priorización ambiental | Completado, preliminar | Reglas explicables; faltan datos ambientales en la mayoría de casos |
| Aplicación local y por lote | Completado | Streamlit disponible en `http://localhost:8501`, con carga y descarga CSV |
| Paquete de publicación | Completado | `Dockerfile`, `render.yaml`, comprobación previa y guía de despliegue |
| Acceso público | En preparación de alojamiento | Falta vincular un host con cuenta y verificar la URL HTTPS |

## Datos y trazabilidad

1. Se descargó la versión ASCII de ECOTOX de septiembre de 2026.
2. Se seleccionaron pruebas con *Daphnia magna*, agua dulce, exposición de 48 h, endpoint EC50 y medición de inmovilización.
3. Se excluyeron valores con operadores de desigualdad, unidades no convertibles y concentraciones no numéricas.
4. Se consolidaron resultados repetidos por CASRN con la mediana de `log10(EC50 mg/L)`.
5. PubChem resolvió estructuras químicas por CASRN; RDKit aceptó 408 de las 675 sustancias.

## Modelo actual

| Elemento | Resultado |
|---|---|
| Algoritmo seleccionado | Extra Trees híbrido: descriptores 2D + fingerprints Morgan |
| Sustancias de entrenamiento | 408 |
| Ensayos incluidos para esas sustancias | 801 |
| Variables estructurales | 10 descriptores moleculares 2D + fingerprint Morgan de 2.048 bits |
| Validación | GroupKFold por andamio de Murcko, 5 particiones |
| MAE promedio | 0.963 log10(mg/L) |
| RMSE promedio | 1.237 log10(mg/L) |
| R² promedio | 0.423 |
| Aporte de representación | 20.4 % descriptores 2D; 79.6 % fingerprints Morgan |

### Interpretación registrada

El modelo puede apoyar una primera priorización de sustancias comparables a las usadas en el entrenamiento. Sus métricas son moderadas; no debe presentarse como una sustitución de un ensayo ecotoxicológico experimental.

## Entregables disponibles

| Entregable | Ubicación |
|---|---|
| Modelo entrenado | `models/qsar_daphnia_ec50.joblib` |
| Conjunto listo para entrenar | `data/processed/qsar_training_ready.csv` |
| Predicciones fuera de muestra | `results/out_of_fold_predictions.csv` |
| Comparación de modelos | `results/model_comparison.csv` |
| Importancia de variables | `results/feature_importance.csv` |
| Figuras de tesis | `results/observed_vs_predicted.png`, `results/residuals_vs_predicted.png`, `results/feature_importance.png` |
| Priorización preliminar | `results/environmental_priority.csv` |
| Interfaz local | `app.py` |

## Pendientes prioritarios

### 1. Auditoría de calidad de datos

- [x] Revisar automáticamente una muestra estratificada de 30 CASRN frente a las tablas originales de ECOTOX.
- [x] Confirmar especie, duración, endpoint, medición, unidad y valor de EC50 de cada registro de la muestra.
- [ ] Revisar los casos con alta dispersión entre estudios (`log_ec50_sd`) y decidir si deben excluirse o documentarse.
- [ ] Conservar fecha, versión y criterios de consulta de ECOTOX en el manuscrito.

### 2. Mejorar y comparar el modelo

- [x] Comparar fingerprints moleculares y un modelo híbrido de fingerprints más descriptores 2D. El modelo híbrido fue seleccionado por su menor MAE.
- [ ] Evaluar ajuste de hiperparámetros dentro de la validación cruzada, sin usar el conjunto de evaluación para escoger parámetros.
- [ ] Comparar el rendimiento con una línea base simple, por ejemplo la media del entrenamiento.
- [ ] Si se consigue una fuente independiente, realizar validación externa por sustancias no usadas durante el desarrollo.

### 3. Completar la priorización ambiental

- [ ] Reunir concentraciones ambientales medidas en Colombia u otra zona de estudio.
- [ ] Verificar porcentajes de remoción por tecnología de tratamiento.
- [ ] Actualizar `priority_input.csv` con esos valores documentados.
- [ ] Convertir la prioridad preliminar en una prioridad completa sólo cuando existan los cuatro componentes: toxicidad, logKow, ocurrencia y remoción.

### 4. Mejorar la herramienta y publicar

- [x] Añadir carga por lote de un CSV con varias sustancias.
- [x] Permitir descargar el resultado desde la interfaz.
- [ ] Añadir validación de CASRN y búsqueda asistida de SMILES.
- [x] Incorporar una pantalla de alcance y limitaciones para usuarios no técnicos.
- [x] Preparar contenedor Docker, configuración de Render y comprobación previa.
- [ ] Vincular un alojamiento público, verificar la URL HTTPS y registrar la versión publicada.

### 5. Preparar la tesis

- [x] Redactar métodos, resultados, alcance y reproducibilidad en `docs/informe_tecnico_qsar.md`.
- [x] Elaborar guía de uso local en `docs/guia_usuario_local.md`.
- [ ] Reportar explícitamente la limitación del R² moderado y la necesidad de validación externa.
- [ ] Añadir discusión sobre el uso de IA como apoyo a la priorización y reducción de ensayos, sin reemplazar evidencia experimental.

### 6. Publicación y acceso externo

- [x] Preparar el despliegue público después de completar auditoría, pruebas y documentación.
- [x] Definir acceso público para cualquier persona con la URL.
- [x] Añadir mensajes de alcance, dominio de aplicabilidad y registro de versión del modelo.
- [ ] Conectar el repositorio a un host público y comprobar el enlace final.

## Próximo hito recomendado

Ejecutar la comprobación previa, vincular la carpeta a un repositorio de GitHub y desplegarla en Render o Streamlit Community Cloud. Antes de compartir la URL se debe probar una predicción individual y un CSV desde el enlace HTTPS.

## Cómo actualizar esta bitácora

Después de cada cambio relevante, registrar fecha, scripts ejecutados, datos añadidos, resultado, decisión tomada y archivos de evidencia. Mantener los resultados reproducibles dentro de `data/processed/`, `models/` y `results/`.
