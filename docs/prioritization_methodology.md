# Priorización ambiental posterior al QSAR

## Propósito

Esta etapa recibe la toxicidad estimada por el QSAR y variables ambientales. No reclasifica ni reentrena el modelo de toxicidad.

## Componentes de la puntuación

Cada componente aporta entre 1 y 3 puntos, donde 3 representa mayor preocupación ambiental.

| Componente | 1 punto | 2 puntos | 3 puntos |
|---|---:|---:|---:|
| EC50 predicha (mg/L) | > 10 | > 1 a 10 | ≤ 1 |
| logKow | < 2 | 2 a < 4 | ≥ 4 |
| Concentración observada (µg/L) | < 0.1 | 0.1 a < 1 | ≥ 1 |
| Remoción en tratamiento (%) | > 80 | > 50 a 80 | ≤ 50 |

La puntuación final es la media de los componentes con dato válido. La clasificación es Baja si la media es menor de 1.67, Media desde 1.67 hasta menor de 2.34 y Alta desde 2.34.

## Datos faltantes y alcance

La toxicidad y logP se calculan para cada SMILES válido. Si falta concentración observada o remoción, el resultado se marca como preliminar y muestra el porcentaje de completitud. Las reglas no sustituyen una evaluación regulatoria ni establecen causalidad.

## Uso

Crear un CSV con estas columnas:

```text
compound_name,chemical_key,smiles,logkow,occurrence_ug_l,removal_pct
```

Luego ejecutar:

```powershell
.\.venv\Scripts\python.exe .\src\prioritize_environmental_risk.py --input .\data\processed\priority_input.csv
```
