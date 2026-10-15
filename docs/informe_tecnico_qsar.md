# Informe técnico — modelo QSAR de priorización ecotoxicológica

Fecha de cierre técnico: 15 de octubre de 2026.

## 1. Finalidad

La herramienta estima la toxicidad acuática aguda de una sustancia química expresada como `log10(EC50 mg/L)` para inmovilización de *Daphnia magna* en agua dulce, a 48 horas. El resultado se incorpora después a una priorización ambiental; no se entrenó el modelo con categorías de riesgo derivadas.

## 2. Datos y control de calidad

La fuente de entrenamiento fue la descarga ASCII de ECOTOX publicada en septiembre de 2026. Se vincularon pruebas, resultados y sustancias mediante los identificadores de ECOTOX. Se conservaron registros con *Daphnia magna*, agua dulce, 48 h, EC50, medición de inmovilización, concentración positiva y unidades normalizables a mg/L.

Se excluyeron intervalos, operadores de desigualdad y unidades no convertibles. Los resultados repetidos se consolidaron por CASRN mediante la mediana de `log10(EC50 mg/L)`.

| Etapa | Resultado |
|---|---:|
| Resultados EC50 previos a control de calidad | 1.860 |
| Registros válidos de inmovilización | 1.449 |
| Sustancias antes de estructura | 675 |
| Sustancias con SMILES RDKit válido | 408 |
| Ensayos asociados a esas sustancias | 801 |

La auditoría automática verificó una muestra estratificada de 30 sustancias contra las tablas originales `tests.txt` y `results.txt` de ECOTOX. Las 30 cumplieron especie, medio, duración, endpoint, medición, valor, operador y unidad. Esta comprobación valida la trazabilidad del procesamiento; no sustituye la lectura individual de las publicaciones científicas originales.

## 3. Representación química y modelo

Cada estructura se representa mediante diez descriptores 2D de RDKit y un fingerprint Morgan de radio 2 con 2.048 bits. Se compararon Extra Trees con descriptores, Extra Trees con fingerprints y el modelo híbrido.

El modelo seleccionado es **Extra Trees híbrido**, porque alcanzó el menor error de validación. El modelo usa aproximadamente 20.4 % de su importancia total en los descriptores 2D y 79.6 % en los fingerprints Morgan.

## 4. Validación

La evaluación utiliza cinco particiones `GroupKFold` agrupadas por andamio de Murcko. Así se evita que sustancias de una misma familia estructural entren a la vez al ajuste y a la evaluación de una misma partición.

| Métrica fuera de muestra | Resultado |
|---|---:|
| MAE | 0.963 log10(mg/L) |
| RMSE | 1.237 log10(mg/L) |
| R² | 0.423 |

El rendimiento es moderado. La herramienta es apropiada para apoyar una primera priorización y para orientar qué sustancias requieren atención, pero no reemplaza un ensayo experimental ni una evaluación regulatoria.

## 5. Dominio de aplicabilidad

La aplicación calcula la similitud de Tanimoto entre la estructura consultada y las estructuras de entrenamiento. Se genera una advertencia cuando la mayor similitud es menor a 0.172, correspondiente al percentil 5 de los vecinos estructurales más cercanos en el entrenamiento. También se revisan los rangos de los descriptores moleculares.

Una predicción con advertencia debe tratarse como exploratoria y no utilizarse como evidencia suficiente para tomar una decisión ambiental.

## 6. Priorización ambiental

La prioridad combina, cuando están disponibles, cuatro componentes: EC50 predicha, logKow, concentración observada y porcentaje de remoción. Cada componente se transforma en una puntuación de 1 a 3; la categoría final se expresa como Baja, Media o Alta.

Si faltan concentración observada o remoción, la aplicación conserva la salida pero la etiqueta como **preliminar** y muestra su porcentaje de completitud. No se reemplazan datos faltantes por cero.

## 7. Uso local

Desde la carpeta del proyecto:

```powershell
.\.venv\Scripts\streamlit.exe run .\app.py
```

La interfaz se abre en `http://localhost:8501`. Se requiere un SMILES válido; los demás campos son opcionales, pero completarlos mejora la priorización ambiental.

## 8. Reproducibilidad

Los scripts principales se encuentran en `src/`. La secuencia completa es:

1. `extract_ecotox_endpoint.py`
2. `build_qsar_training_table.py`
3. `resolve_training_structures.py`
4. `finalize_training_dataset.py`
5. `train_hybrid_qsar.py`
6. `create_qsar_figures.py`
7. `run_ecotox_audit.py`

El modelo entrenado, las tablas procesadas, las métricas y las figuras se guardan en `models/`, `data/processed/` y `results/`.

## 9. Condiciones antes de publicación pública

Antes de habilitar acceso externo se debe mantener este informe junto con la metodología, incluir el alcance y las limitaciones en la interfaz, probar carga por lote y descarga de resultados, y decidir el tipo de acceso. La publicación no debe presentar la herramienta como un sustituto de evidencia experimental.
