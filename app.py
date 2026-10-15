"""Interfaz local para predicción QSAR y priorización ambiental."""

from __future__ import annotations

from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import streamlit as st

from src.prioritize_environmental_risk import (
    ecotox_score,
    evaluate_substance,
    logkow_score,
    occurrence_score,
    priority_label,
    removal_score,
)
from src.train_qsar import descriptors_from_smiles, nearest_training_similarity, prediction_frame_from_smiles


ROOT = Path(__file__).resolve().parent
MODEL_PATH = ROOT / "models" / "qsar_daphnia_ec50.joblib"


@st.cache_resource
def load_model(model_version: float) -> dict[str, object]:
    """Recarga el modelo si el archivo entrenado cambia."""
    return joblib.load(MODEL_PATH)


def number_or_none(value: str) -> float | None:
    value = value.strip().replace(",", ".")
    if not value:
        return None
    try:
        return float(value)
    except ValueError:
        return None


st.set_page_config(page_title="Priorización ecotoxicológica", page_icon="🧪", layout="centered")
st.title("Priorización ecotoxicológica")
st.caption("QSAR: EC50 de inmovilización en Daphnia magna, agua dulce, 48 h.")

if not MODEL_PATH.exists():
    st.error("No se encontró el modelo entrenado. Ejecuta primero src/train_qsar.py.")
    st.stop()

with st.form("substance"):
    compound_name = st.text_input("Nombre de la sustancia", placeholder="Ej. Paracetamol")
    chemical_key = st.text_input("CAS o identificador", placeholder="Ej. 103-90-2")
    smiles = st.text_area("SMILES *", placeholder="Ej. CC(=O)NC1=CC=C(O)C=C1")
    st.markdown("#### Datos ambientales opcionales")
    logkow = st.text_input("logKow medido", placeholder="Ej. 0.5")
    occurrence = st.text_input("Concentración observada (µg/L)", placeholder="Ej. 0.2")
    removal = st.text_input("Remoción en tratamiento (%)", placeholder="Ej. 85")
    submitted = st.form_submit_button("Predecir y priorizar")

if submitted:
    values = descriptors_from_smiles(smiles)
    if values is None:
        st.error("El SMILES no es válido. Verifica la estructura química.")
        st.stop()

    bundle = load_model(MODEL_PATH.stat().st_mtime)
    x, values = prediction_frame_from_smiles(smiles, bundle)
    log_ec50 = float(bundle["model"].predict(x)[0])
    ec50 = 10 ** log_ec50
    measured_logkow = number_or_none(logkow)
    observed_occurrence = number_or_none(occurrence)
    observed_removal = number_or_none(removal)
    if observed_removal is not None and not 0 <= observed_removal <= 100:
        st.error("La remoción debe estar entre 0 y 100 %.")
        st.stop()

    logkow_used = measured_logkow if measured_logkow is not None else float(values["log_p"])
    scores: dict[str, float] = {
        "Ecotoxicidad": ecotox_score(ec50),
        "logKow": logkow_score(logkow_used),
    }
    if observed_occurrence is not None and observed_occurrence >= 0:
        scores["Ocurrencia"] = occurrence_score(observed_occurrence)
    if observed_removal is not None:
        scores["Remoción"] = removal_score(observed_removal)

    priority_score = float(np.mean(list(scores.values())))
    outside = [
        feature for feature in bundle["descriptor_features"]
        if values[feature] < bundle["feature_ranges"][feature]["min"]
        or values[feature] > bundle["feature_ranges"][feature]["max"]
    ]
    similarity = nearest_training_similarity(smiles, bundle["reference_fingerprints"])
    similarity_warning = similarity < bundle["applicability_similarity_threshold"]
    completeness = len(scores) / 4
    left, middle, right = st.columns(3)
    left.metric("EC50 estimada", f"{ec50:.3g} mg/L")
    middle.metric("Prioridad", priority_label(priority_score))
    right.metric("Completitud ambiental", f"{completeness:.0%}")
    st.write(f"**Sustancia:** {compound_name or 'Sin nombre'} · **Identificador:** {chemical_key or 'No informado'}")
    st.write(f"Puntuación promedio: **{priority_score:.2f} / 3**. logKow usado: **{logkow_used:.2f}** ({'medido' if measured_logkow is not None else 'estimado por RDKit'}).")
    st.dataframe(pd.DataFrame([scores]), hide_index=True, width="stretch")
    if completeness < 1:
        st.info("Prioridad preliminar: añade concentración observada y remoción para completar la evaluación ambiental.")
    st.write(
        f"Similitud estructural máxima con el entrenamiento: **{similarity:.2f}** "
        f"(umbral mínimo: {bundle['applicability_similarity_threshold']:.2f})."
    )
    if outside or similarity_warning:
        reasons = []
        if outside:
            reasons.append("descriptores fuera de rango: " + ", ".join(outside))
        if similarity_warning:
            reasons.append("estructura poco similar a las sustancias de entrenamiento")
        st.warning("Interpreta la predicción con cautela: " + "; ".join(reasons) + ".")
    else:
        st.success("Los descriptores evaluados están dentro de los rangos observados durante el entrenamiento.")

st.divider()
st.subheader("Procesar un archivo CSV")
st.caption("Columnas obligatorias: `compound_name`, `chemical_key`, `smiles`. Opcionales: `logkow`, `occurrence_ug_l`, `removal_pct`.")
uploaded = st.file_uploader("Carga un CSV", type="csv")
if uploaded is not None:
    if uploaded.size > 10 * 1024 * 1024:
        st.error("El archivo supera el máximo de 10 MB.")
    else:
        try:
            batch = pd.read_csv(uploaded)
        except (UnicodeDecodeError, pd.errors.ParserError) as error:
            st.error(f"No se pudo leer el CSV: {error}")
        else:
            required_columns = {"compound_name", "chemical_key", "smiles"}
            missing_columns = required_columns - set(batch.columns)
            if missing_columns:
                st.error(f"Faltan columnas obligatorias: {', '.join(sorted(missing_columns))}.")
            elif len(batch) == 0:
                st.warning("El archivo no contiene filas.")
            elif len(batch) > 1_000:
                st.error("El archivo supera el máximo de 1.000 filas por ejecución.")
            elif st.button("Procesar archivo"):
                bundle = load_model(MODEL_PATH.stat().st_mtime)
                results = pd.DataFrame([evaluate_substance(bundle, row) for row in batch.to_dict("records")])
                st.success(f"Se procesaron {len(results)} sustancias.")
                st.dataframe(results, hide_index=True, width="stretch")
                st.download_button(
                    "Descargar resultados CSV",
                    data=results.to_csv(index=False).encode("utf-8-sig"),
                    file_name="priorizacion_ecotoxicologica.csv",
                    mime="text/csv",
                )

with st.expander("Alcance y limitaciones del modelo"):
    st.markdown(
        """
        - El modelo estima `log10(EC50 mg/L)` para inmovilización de *Daphnia magna* en agua dulce a 48 h.
        - Se validó por familias estructurales; su desempeño de validación cruzada fue MAE 0,963 y R² 0,423.
        - La similitud de Tanimoto y los rangos de descriptores indican si una predicción está dentro del dominio de aplicabilidad.
        - La prioridad ambiental integra la estimación QSAR con logKow y, cuando se proporcionan, ocurrencia y remoción. No sustituye ensayos experimentales ni una evaluación regulatoria.
        """
    )
