# Priorización ambiental de contaminantes emergentes domésticos en Colombia

## Propósito

Esta etapa recibe la toxicidad estimada por el QSAR y variables ambientales del contexto colombiano. No reclasifica ni reentrena el modelo de toxicidad.

## Separación necesaria: peligro y riesgo local

La descarga de ECOTOX se utiliza únicamente como evidencia de **peligro ecotoxicológico intrínseco**: una concentración de efecto medida bajo condiciones experimentales definidas. Que los estudios provengan de otros países no invalida ese componente si el organismo, endpoint, duración, medio y unidad se controlan; una EC50 o LC50 no es una propiedad del país donde se consultó la base.

Sin embargo, ECOTOX no permite afirmar el **riesgo ambiental en Colombia**. La presencia real de contaminantes domésticos, sus concentraciones en afluentes, efluentes y cuerpos receptores, las prácticas de consumo y la eficiencia de remoción de cada PTAR sí dependen del territorio y del proceso. Por ello, la aplicación solo debe describir una prioridad colombiana completa cuando `occurrence_ug_l` y `removal_pct` proceden de datos locales, con fuente, lugar, fecha y matriz documentados. En ausencia de esos datos, el resultado se denomina **prioridad preliminar basada en peligro**, no riesgo local.

## Componentes de la puntuación

Cada componente aporta entre 1 y 3 puntos, donde 3 representa mayor preocupación ambiental.

| Componente | 1 punto | 2 puntos | 3 puntos |
|---|---:|---:|---:|
| EC50 predicha (mg/L) | > 10 | > 1 a 10 | ≤ 1 |
| logKow | < 2 | 2 a < 4 | ≥ 4 |
| Concentración observada en Colombia (µg/L) | < 0.1 | 0.1 a < 1 | ≥ 1 |
| Remoción en la PTAR/proceso colombiano de interés (%) | > 80 | > 50 a 80 | ≤ 50 |

La puntuación final es la media de los componentes con dato válido. La clasificación es Baja si la media es menor de 1.67, Media desde 1.67 hasta menor de 2.34 y Alta desde 2.34.

## Datos faltantes y alcance

La toxicidad y logP se calculan para cada SMILES válido. Si falta concentración observada en Colombia o remoción del proceso evaluado, el resultado se marca como preliminar y muestra el porcentaje de completitud. Las reglas no sustituyen una evaluación regulatoria ni establecen causalidad.

Junto con cada valor local se debe conservar, fuera o dentro del CSV, la referencia, municipio/cuenca, matriz (afluente, efluente o agua superficial), fecha de muestreo, método analítico y unidad original. No se deben reutilizar como si fueran colombianos valores medidos en otro territorio.

## Cociente de riesgo ecotoxicológico de cribado

Cuando existe una concentración ambiental local compatible con la matriz y el sitio evaluados, la aplicación calcula por separado para cada grupo:

`RQ = concentración ambiental (mg/L) / concentración de efecto predicha (mg/L)`

La salida no mezcla EC50 de crustáceos, LC50 de peces ni EC50 de algas. Muestra un RQ independiente por grupo, junto con el endpoint, el desempeño de validación y el dominio de aplicabilidad. Los intervalos de lectura inicial son `RQ < 0,1` bajo, `0,1 ≤ RQ < 1` intermedio y `RQ ≥ 1` alto, todos con la etiqueta **de cribado**. Estos intervalos no constituyen límites regulatorios colombianos. La interpretación final requiere confirmar el endpoint, la matriz, la representatividad temporal y el criterio de protección aplicable.

Si no se suministra exposición local, la aplicación muestra peligro QSAR y no calcula riesgo. Si una predicción está fuera del dominio químico o su grupo no superó la validación, el RQ de ese grupo se marca como no confiable.

## Uso

Crear un CSV con estas columnas:

```text
compound_name,chemical_key,smiles,logkow,occurrence_ug_l,removal_pct
```

Luego ejecutar:

```powershell
.\.venv\Scripts\python.exe .\src\prioritize_environmental_risk.py --input .\data\processed\priority_input.csv
```
