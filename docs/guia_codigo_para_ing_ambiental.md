# Guía sencilla del código y de la aplicación

Esta gl quuía explica el proyecto para una persona de ingeniería ambientae no programa. No es necesario memorizar Python: cada archivo representa una etapa de una evaluación ecotoxicológica.

## La idea antes del código

El sistema no crea datos experimentales nuevos. Usa ensayos ecotoxicológicos ya existentes para aprender una relación entre la **estructura química** de una sustancia y su concentración de efecto.

```mermaid
flowchart LR
    A[Datos ECOTOX y base consolidada] --> B[Depuración y selección de ensayos comparables]
    B --> C[CAS y estructura SMILES verificable]
    C --> D[Descriptores y fingerprints]
    D --> E[Modelo por grupo acuático]
    E --> F[Predicción de una sustancia nueva]
    F --> G[Prioridad ambiental]
```

El modelo se entrenó una sola vez. Para consultar una sustancia nueva no hace falta cargar miles de datos: basta con su SMILES; los datos de ocurrencia, logKow y remoción mejoran la etapa de prioridad ambiental.

## Glosario mínimo

| Término | Explicación sencilla |
|---|---|
| Python | Lenguaje con el que se escribió el proyecto. |
| Script | Archivo que ejecuta una tarea concreta, por ejemplo limpiar datos o entrenar. |
| CSV | Tabla de datos, similar a Excel, separada por comas. |
| CAS | Código identificador de una sustancia química. |
| SMILES | Texto que representa la estructura química de una molécula. |
| QSAR | Modelo que relaciona estructura química y actividad o toxicidad. |
| Descriptor | Número calculado desde el SMILES, por ejemplo masa molecular o polaridad. |
| Fingerprint | Cadena numérica que resume fragmentos estructurales de una molécula. |
| Entrenamiento | Etapa en que el algoritmo aprende usando sustancias con resultados conocidos. |
| Predicción | Estimación para una sustancia sin ensayo experimental dentro de la aplicación. |

## 1. El archivo de la aplicación: `app.py`

Este es el archivo que crea la pantalla que ve la persona usuaria. No entrena el modelo cada vez: carga el modelo ya entrenado y lo utiliza.

### Líneas 1 a 22: preparación

- Línea 1: describe que el archivo corresponde a la interfaz multi-organismo.
- Líneas 5 a 10: importan librerías. `pandas` trabaja con tablas; `numpy` hace cálculos numéricos; `streamlit` crea la página web; `joblib` abre el modelo guardado.
- Línea 12: importa las funciones que predicen para crustáceos, peces y algas.
- Línea 13: importa las reglas usadas para convertir toxicidad, logKow, ocurrencia y remoción en prioridad ambiental.
- Líneas 16 y 17: localizan el archivo `qsar_multiorganismo.joblib`. Este archivo es el modelo aprendido; no es una nueva base de datos.
- Líneas 20 a 22: `load_model` abre el modelo una vez y lo mantiene en memoria para no cargarlo nuevamente en cada clic.

### Líneas 25 a 32: lectura segura de números
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
La función `as_number` transforma textos como `0,5` o `0.5` en el número 0.5. Si una casilla se deja vacía, devuelve “sin dato” en lugar de inventar un cero.

### Líneas 35 a 53: cálculo de prioridad ambiental

`environmental_priority` recibe una concentración de efecto predicha y datos ambientales opcionales.

- Líneas 36 a 38: obtiene logKow, ocurrencia y remoción del formulario.
- Línea 39: si no se proporcionó logKow medido, usa el logP calculado desde la estructura. Se indica después que fue estimado.
- Línea 40: inicia dos puntajes: efecto ecotoxicológico y logKow.
- Líneas 41 a 44: solo añade ocurrencia y remoción cuando fueron informadas y son valores coherentes.
- Línea 45: calcula el promedio de los puntajes disponibles.
- Líneas 46 a 53: devuelve el nivel de prioridad, el puntaje, la completitud y el origen del logKow.

La idea ambiental es clara: una sustancia puede ser tóxica, pero su prioridad real en Colombia aumenta si aparece en agua colombiana y si la PTAR o proceso local la remueve con dificultad.

### Líneas 56 a 80: procesamiento de una fila de CSV

Esta función permite procesar varias sustancias a la vez.

- Línea 57: conserva primero los datos originales de la fila.
- Líneas 58 a 62: intenta leer el SMILES. Si es inválido, informa el error y no inventa una predicción.
- Línea 63: reúne las predicciones de todos los grupos acuáticos.
- Líneas 64 a 71: agrega a la salida el grupo más sensible y una columna por cada grupo.
- Líneas 72 a 79: solo calcula prioridad general si existe una estimación general confiable.

### Líneas 83 a 99: construcción del formulario

`st` significa Streamlit. Estas líneas dibujan la página:

- Líneas 83 a 85: título y explicación corta.
- Líneas 87 a 89: detienen la aplicación si el archivo del modelo no está disponible.
- Líneas 91 a 99: crean los campos visibles. Solo el SMILES tiene asterisco porque representa la estructura necesaria para el QSAR. Las otras variables ayudan a priorizar, pero son opcionales.

### Líneas 101 a 137: predicción individual

- Líneas 102 y 103: reúnen lo digitado en una sola ficha de sustancia.
- Líneas 104 a 107: revisan que la remoción esté entre 0 % y 100 %.
- Línea 108: carga el modelo ya entrenado.
- Líneas 109 a 113: pide predicciones por grupo y detiene el proceso si el SMILES no representa una molécula válida.
- Línea 114: busca la conclusión ecotoxicológica general.
- Líneas 116 a 125: cuando existe una conclusión confiable, presenta concentración de efecto, grupo más sensible, prioridad y los puntajes usados.
- Líneas 129 a 137: muestran la tabla por grupo y una advertencia si la molécula está fuera del dominio de aplicabilidad.

### Líneas 139 a 166: carga masiva de CSV

Esta parte permite subir un archivo con muchas sustancias.

- Líneas 141 y 152: definen las tres columnas obligatorias: nombre, identificador y SMILES.
- Líneas 144 a 159: revisan tamaño máximo, lectura correcta, columnas y máximo de 1.000 filas.
- Líneas 160 a 166: evalúan una por una las sustancias, muestran una tabla y permiten descargar el resultado.

### Líneas 168 a 174: límites científicos

Esta sección evita una interpretación incorrecta: el modelo de algas se muestra, pero no decide la prioridad general porque su validación actual no fue suficiente.

## 2. Cómo se depura la base: `build_multigroup_ecotox_dataset.py`

Este archivo recibe la descarga completa de ECOTOX y construye tablas comparables. Lea limpieza no busca aumentar artificialmente el númro de filas: elimina registros que no responden a la misma pregunta ecotoxicológica.

### Líneas 19 a 23: ubicación de archivos

Indican dónde está el archivo descargado, dónde se guardan las tablas depuradas y dónde se guardan los informes.

### Líneas 25 a 52: definición de grupos

`GROUPS` es la regla científica del proyecto. No es una clasificación automática; es una decisión metodológica explícita.

| Grupo | Filtrado aplicado | Razón |
|---|---|---|
| Crustáceos | EC50, inmovilización, 48 h, agua dulce | Compara respuestas del mismo tipo. |
| Peces | LC50, mortalidad, 96 h, agua dulce | Conserva una condición aguda comparable. |
| Algas | EC50, crecimiento o respuesta poblacional, 72 h, agua dulce | Representa productores primarios con una condición comparable. |

Los códigos como `ITX`, `MOR`, `IMBL` y `MORT` pertenecen a la codificación de ECOTOX. El programa los usa para seleccionar solo el efecto correcto.

### Líneas 55 a 61: limpieza de textos y CAS

- Líneas 55 y 56: quitan espacios innecesarios de las columnas de texto.
- Líneas 59 a 61: convierten un CAS que viene solo con números al formato usual con guiones. Si no tiene suficientes dígitos, se descarta como identificador inválido.

### Líneas 64 a 81: conversión de unidades

Todas las concentraciones deben quedar en mg/L para poder compararse. Estas líneas aceptan unidades claras como mg/L, µg/L, ng/L, ppm y ppb; una unidad ambigua se transforma en valor vacío y no pasa al entrenamiento. Esto evita comparar, por ejemplo, un valor en µg/L como si fuera mg/L.

### Líneas 84 a 113: lectura de pruebas

- Línea 85: verifica que la descarga de ECOTOX exista.
- Líneas 91 a 100: abre el archivo ZIP sin descomprimirlo manualmente y lee la tabla de pruebas por bloques de 100.000 filas. Es una forma eficiente de manejar bases grandes.
- Líneas 102 a 109: asignan cada especie a su grupo y conservan únicamente agua dulce y la duración exacta de cada grupo.
- Líneas 112 y 113: unen todas las pruebas aceptadas y guardan sus identificadores para buscar sus resultados después.

### Líneas 115 a 153: unión con resultados y control de calidad

- Líneas 115 a 128: leen resultados de toxicidad y conservan los que pertenecen a las pruebas seleccionadas.
- Líneas 130 a 137: añaden nombre químico y DTXSID usando el CAS.
- Líneas 138 a 146: comprueban endpoint, efecto y medida correcta para cada grupo.
- Líneas 147 a 152: convierten la concentración, validan el CAS, rechazan resultados censurados como `>10` o `<1` y aplican logaritmo base 10 a la concentración.
- Línea 153: guarda una descripción entendible del endpoint junto a cada registro.

### Líneas 156 a 189: consolidación y trazabilidad

Una sustancia puede tener varios ensayos válidos. Las líneas 156 a 171 calculan la mediana de los resultados en escala logarítmica por CAS y por grupo. También conservan cuántos ensayos y referencias existían. Las líneas 173 a 189 generan un informe con conteos por grupo, por lo que todo el proceso es auditable.

## 3. Qué datos se tomaron realmente para entrenar

Es importante diferenciar la **base consolidada de la tesis** de la **base final de entrenamiento**.

### Datos aprovechados del consolidado

El consolidado entregado por los grupos no contenía SMILES para las sustancias. Sí aportaba, cuando estaban disponibles, los siguientes campos:

| Campo del consolidado | Uso dentro del proyecto |
|---|---|
| `Contaminante_Normalizado` | Nombre legible y revisión de identidad. |
| `CAS` | Identificación inequívoca de la sustancia y clave para buscar estructura. |
| `Especie en el ensayo` | Reconocer el organismo evaluado. |
| `Endpoint` | Distinguir EC50, LC50, NOEC y otros resultados. |
| `Efecto_medido` | Diferenciar inmovilización, mortalidad, crecimiento, etc. |
| `Duración del ensayo` | Conservar solo exposiciones comparables. |
| `Agua dulce o marina` | Seleccionar agua dulce. |
| `EC50` y `Unidad_EC50` | Concentración experimental y su unidad. |
| `Referencia` | Trazabilidad bibliográfica. |

El consolidado fue útil para organizar el problema, reconocer los contaminantes emergentes de interés doméstico y conservar fuentes. Sin embargo, para el entrenamiento multi-organismo definitivo se usó la descarga ASCII de **ECOTOX**, porque contiene más registros estandarizados y campos técnicos verificables. Esta decisión tiene un límite importante: ECOTOX aporta la parte de **peligro ecotoxicológico**, no demuestra por sí solo el riesgo en Colombia.

Una concentración de efecto experimental puede usarse como evidencia de peligro aunque el estudio se haya realizado en otro país, siempre que se conserven comparables el organismo, endpoint, duración, medio y unidad. Lo que no puede trasladarse automáticamente a Colombia es la exposición: qué contaminantes de productos de aseo, cuidado personal o medicamentos aparecen realmente; en qué concentración; y qué tan bien los remueve una PTAR colombiana. Por eso, la prioridad ambiental completa se calcula solo al añadir ocurrencia y remoción locales documentadas. Sin ellas, la salida es una **prioridad preliminar basada en peligro**, no una clasificación de riesgo colombiano definitiva.

### Datos usados directamente en el entrenamiento multi-organismo

De ECOTOX se tomaron estas variables:

| Campo ECOTOX | Para qué se usó |
|---|---|
| `test_cas` | Reconstruir el CAS de cada sustancia. |
| `chemical_name` y `dtxsid` | Trazabilidad química. |
| `species_number` | Asignar el ensayo a crustáceos, peces o algas. |
| `media_type` | Mantener agua dulce (`FW` o `FW/`). |
| `exposure_duration_mean` y unidad | Seleccionar 48 h, 96 h o 72 h según el grupo. |
| `endpoint` | Seleccionar EC50 o LC50 según el grupo. |
| `effect` y `measurement` | Confirmar que el efecto fuera inmovilización, mortalidad o respuesta algal. |
| `conc1_mean`, `conc1_unit` y `conc1_mean_op` | Tomar la concentración, convertirla a mg/L y excluir límites como `>10` o `<1`. |
| `reference_number` | Contar referencias y conservar trazabilidad. |

Después de la depuración quedaron 1.793 combinaciones sustancia–grupo antes de validar estructuras: 853 de crustáceos, 800 de peces y 140 de algas. Tras validar la estructura química quedaron 816, 756 y 136 sustancias, respectivamente.

### Por qué el SMILES no estaba en la base y cómo se añadió

El SMILES no es un resultado ecotoxicológico; es una representación de estructura molecular. Por eso muchas bases ambientales contienen CAS y nombre, pero no SMILES.

El flujo fue:

1. Se tomó el CAS disponible en ECOTOX o en la base consolidada.
2. Se consultó PubChem usando ese CAS.
3. PubChem devolvió un SMILES y, cuando fue posible, su CID.
4. RDKit intentó leer el SMILES.
5. Solo las estructuras que RDKit pudo interpretar se marcaron como `rdkit_valid` y entraron al entrenamiento.

No se inventaron SMILES ni se asumió que sustancias con nombres parecidos eran idénticas. Las sustancias sin CAS recuperable, sin coincidencia en PubChem o con SMILES inválido se excluyeron del entrenamiento, pero pueden permanecer en los archivos de trazabilidad.

## 4. Cómo se obtienen las estructuras: `resolve_multigroup_structures.py`

El modelo no puede aprender solo de nombres; necesita estructuras químicas válidas.

- Líneas 14 a 19: definen las tablas de entrada, salida y el máximo de tres consultas simultáneas.
- Líneas 22 a 36: buscan estructuras que ya habían sido validadas en el modelo anterior para no repetir consultas.
- Líneas 39 a 46: leen las sustancias únicas y detectan cuáles todavía no tienen SMILES.
- Líneas 47 a 52: consultan PubChem por CAS. Cada resultado se valida con RDKit dentro de la función reutilizada `query_pubchem`.
- Líneas 53 a 55: juntan los SMILES obtenidos con la tabla ecotoxicológica y crean el archivo final con estructuras.
- Líneas 56 a 65: escriben un reporte con sustancias consultadas, reutilizadas y validadas.

Una estructura solo entra si RDKit puede interpretarla. Así se evita entrenar con moléculas mal escritas o indefinidas.

## 5. Cómo aprende el modelo: `train_multigroup_qsar.py`

### Líneas 22 a 26: reglas mínimas

Estas líneas fijan la entrada, el nombre del modelo guardado y un mínimo de 100 sustancias con estructura válida por grupo. No garantizan calidad por sí mismas, pero impiden entrenar con grupos demasiado pequeños.

### Líneas 29 a 72: función `train_group`

Esta función repite el mismo proceso para cada grupo.

- Línea 30: llama a `prepare_features`, que obtiene descriptores, SMILES y familias estructurales.
- Líneas 31 y 32: detienen el entrenamiento si hay menos de 100 sustancias válidas.
- Líneas 33 y 34: construyen fingerprints Morgan de 2.048 posiciones. Cada posición indica la presencia o ausencia de un patrón estructural.
- Líneas 35 a 38: unen los 10 descriptores 2D y los fingerprints para formar `X`, las variables de entrada. `y` es el resultado experimental en escala logarítmica.
- Línea 39: crea Extra Trees con 400 árboles. Cada árbol aprende reglas distintas y el resultado final promedia todos ellos.
- Línea 40: realiza validación cruzada por familias estructurales. Así, moléculas de la misma familia no quedan simultáneamente en entrenamiento y prueba.
- Línea 41: entrena finalmente con todos los datos válidos una vez terminada la evaluación.
- Línea 42: calcula el umbral de similitud química que delimita el dominio de aplicabilidad.
- Líneas 43 a 55: guardan modelo, variables usadas, rangos observados y fingerprints de referencia.
- Líneas 56 a 66: calculan métricas de desempeño: MAE, RMSE y R².
- Líneas 67 a 72: conservan las predicciones realizadas durante la validación para poder auditar los errores.

### Líneas 75 a 122: entrenamiento completo

- Líneas 79 a 85: conservan solamente filas con SMILES válido y resultado de toxicidad disponible.
- Líneas 91 a 96: separan la base por grupo y entrenan un modelo para cada uno.
- Líneas 100 a 107: construyen la “tarjeta del modelo”. La regla de la línea 103 permite integrar solo grupos con R² mayor que cero. Actualmente incluye crustáceos y peces; algas queda visible, pero complementario.
- Líneas 108 a 121: guardan métricas, predicciones de validación y el informe técnico.

## 6. Las propiedades químicas: `train_qsar.py`

Este archivo contiene funciones reutilizadas por el entrenamiento y la aplicación.

- `descriptors_from_smiles`: calcula diez propiedades desde el SMILES: peso molecular, logP, área polar, donadores y aceptores de puentes de hidrógeno, enlaces rotables, átomos pesados, anillos, anillos aromáticos y fracción de carbono sp3.
- `fingerprint_from_smiles`: convierte la estructura en un fingerprint Morgan de 2.048 bits.
- `scaffold_group`: identifica la familia estructural de la molécula.
- `cross_validate`: divide el conjunto por familias estructurales, entrena con unas familias y evalúa con otras.
- `nearest_training_similarity`: mide la similitud de Tanimoto entre una molécula nueva y las moléculas conocidas.

## 7. Cómo se predice por grupo: `multi_group_prediction.py`

### Líneas 11 a 45: `predict_by_group`

- Línea 15: recorre los tres modelos guardados.
- Línea 16: transforma el SMILES nuevo en las mismas variables usadas durante el entrenamiento.
- Líneas 18 y 19: el modelo predice primero el logaritmo de la concentración y luego se convierte a mg/L.
- Líneas 20 a 24: revisan si las propiedades de la sustancia quedan fuera de los rangos conocidos.
- Líneas 25 y 26: calculan similitud de Tanimoto y recuperan el umbral del grupo.
- Líneas 28 a 41: construyen una fila de resultado con endpoint, concentración, métricas, similitud y advertencias.

### Líneas 48 a 66: `general_ecotoxicity`

Esta es la integración conservadora:

1. Conserva grupos con validación aceptable y molécula dentro del dominio.
2. Si no queda ninguno, informa que no hay estimación general confiable.
3. Si hay varios, selecciona la menor concentración de efecto porque representa el grupo potencialmente más sensible.

No mezcla EC50 y LC50 para volverlas un dato experimental único. Las usa como salidas independientes y aplica una regla de priorización posterior.

## 8. Regla de prioridad: `prioritize_environmental_risk.py`

Este archivo define reglas simples y transparentes, no entrena otra IA.

| Variable | Puntaje alto | Sentido ambiental |
|---|---|---|
| Concentración de efecto | Menor o igual a 1 mg/L | Efecto potencial a baja concentración. |
| logKow | Mayor o igual a 4 | Mayor afinidad por fases orgánicas. |
| Ocurrencia | Mayor o igual a 1 µg/L | Mayor presencia medida en agua. |
| Remoción | Menor o igual a 50 % | Menor eliminación en tratamiento. |

Los puntajes se promedian: 1.00–1.66 es prioridad baja, 1.67–2.33 media y 2.34–3.00 alta. Si faltan ocurrencia o remoción, el resultado se denomina preliminar.

## Cómo explicarlo en una sustentación

> “Primero organizamos los ensayos por grupos acuáticos y condiciones comparables. Después validamos identidad química y unidades. Con las estructuras SMILES calculamos propiedades moleculares y entrenamos un modelo independiente para crustáceos, peces y algas. La aplicación recibe una molécula nueva, estima una concentración de efecto para cada grupo, revisa si se parece a los compuestos de entrenamiento y, con los grupos confiables, genera una prioridad ambiental que puede complementarse con ocurrencia y remoción en planta.”

## Qué no debe afirmarse

- No se deben presentar las predicciones como reemplazo de un bioensayo.
- No se debe decir que el modelo predice todas las especies acuáticas.
- No se debe afirmar una prioridad completa si faltan datos de ocurrencia y remoción.
- No se debe usar una predicción fuera del dominio de aplicabilidad como conclusión definitiva.
