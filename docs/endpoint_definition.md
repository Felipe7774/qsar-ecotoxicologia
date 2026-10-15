# Endpoint del modelo QSAR

## Pregunta predictiva

Estimar la toxicidad acuática aguda de un compuesto químico como `log10(EC50_mg_L)` para *Daphnia magna*.

## Endpoint seleccionado

- Organismo: *Daphnia magna*.
- Medio: agua dulce.
- Duración de exposición: 48 horas.
- Efecto: inmovilización.
- Medida: EC50.
- Unidad original permitida: mg/L o una unidad convertible de forma documentada a mg/L.
- Variable objetivo: `log_ec50_mg_l = log10(ec50_mg_l)`.

El logaritmo reduce el efecto de concentraciones extremas y permite entrenar una regresión más estable.

## Criterios de inclusión

1. Sustancia única, identificada con CASRN y estructura química o SMILES.
2. Ensayo de *Daphnia magna* con exposición de 48 h.
3. Resultado EC50 numérico, positivo y expresado en concentración.
4. Unidad convertible a mg/L con el factor de conversión registrado.
5. Referencia o identificador de fuente conservado.

## Criterios de exclusión

1. Mezclas, formulaciones comerciales o sustancias sin estructura química definida.
2. Ensayos con otra especie, duración o efecto biológico.
3. Valores cualitativos, intervalos sin valor central o resultados sin unidad.
4. Registros duplicados del mismo ensayo.

## Consolidación por compuesto

Un mismo compuesto puede tener varios ensayos válidos. Se conservarán los registros crudos y, para el entrenamiento inicial, se calculará la mediana de `log_ec50_mg_l` por estructura química canónica. La cantidad de ensayos y la variabilidad entre estudios se conservarán como campos de trazabilidad.

## Alcance

El modelo estimará toxicidad acuática aguda para el endpoint definido. La posterior priorización ambiental en seis categorías será una capa de interpretación separada y no la variable objetivo con la que se entrena el modelo.
