# Endpoints del modelo QSAR multi-organismo

## Pregunta predictiva

Estimar concentraciones de efecto agudo para contaminantes emergentes en grupos de organismos acuáticos y priorizar la sustancia usando el grupo más sensible con validación y dominio de aplicabilidad aceptables.

## Endpoints comparables

| Grupo | Ensayo incluido | Variable objetivo |
|---|---|---|
| Crustáceos | EC50 de inmovilización, 48 h, agua dulce | `log10(EC50_mg_L)` |
| Peces | LC50 de mortalidad, 96 h, agua dulce | `log10(LC50_mg_L)` |
| Algas | EC50 de respuesta poblacional o crecimiento, 72 h, agua dulce | `log10(EC50_mg_L)` |

Cada grupo se entrena en un modelo independiente. No se combinan EC50 y LC50 ni respuestas de organismos distintos en una sola variable de entrenamiento.

## Criterios de inclusión

1. CASRN y SMILES químicamente válido por RDKit.
2. Duración, medio acuático, endpoint y tipo de efecto explícitos para cada grupo.
3. Concentración positiva y unidad convertible de forma documentada a mg/L.
4. Resultado numérico sin operadores de censura (`>`, `<`, `~`) y con trazabilidad de ECOTOX.

## Integración ambiental

La aplicación muestra el resultado de cada grupo. La prioridad general usa la menor concentración de efecto entre los grupos cuyo desempeño de validación es aceptable y cuya molécula está dentro del dominio de aplicabilidad. Un grupo con desempeño insuficiente se conserva como evidencia complementaria, pero no determina la prioridad.

## Alcance

La salida es una estimación QSAR de apoyo a decisión para organismos acuáticos. No sustituye ensayos experimentales ni permite afirmar toxicidad en todas las especies acuáticas.
