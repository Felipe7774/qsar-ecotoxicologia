"""Interfaz para priorización ecotoxicológica multi-organismo."""

from __future__ import annotations

from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import streamlit as st

from src.multi_group_prediction import general_ecotoxicity, predict_by_group
from src.prioritize_environmental_risk import ecotox_score, logkow_score, occurrence_score, priority_label, removal_score


ROOT = Path(__file__).resolve().parent
MODEL_PATH = ROOT / "models" / "qsar_multiorganismo.joblib"


@st.cache_resource
def load_model(model_version: float) -> dict[str, object]:
    return joblib.load(MODEL_PATH)


def as_number(value: object) -> float | None:
    text = str(value).strip().replace(",", ".")
    if not text or text.casefold() == "nan":
        return None
    try:
        return float(text)
    except ValueError:
        return None


def environmental_priority(concentration: float, descriptors: dict[str, float], source: dict[str, object]) -> dict[str, object]:
    measured_logkow = as_number(source.get("logkow", ""))
    occurrence = as_number(source.get("occurrence_ug_l", ""))
    removal = as_number(source.get("removal_pct", ""))
    logkow = measured_logkow if measured_logkow is not None else float(descriptors["log_p"])
    scores = {"Efecto ecotoxicológico": ecotox_score(concentration), "logKow": logkow_score(logkow)}
    if occurrence is not None and occurrence >= 0:
        scores["Ocurrencia"] = occurrence_score(occurrence)
    if removal is not None and 0 <= removal <= 100:
        scores["Remoción"] = removal_score(removal)
    score = float(np.mean(list(scores.values())))
    return {
        "scores": scores,
        "priority": priority_label(score),
        "score": score,
        "completeness": len(scores) / 4,
        "logkow": logkow,
        "logkow_source": "medido" if measured_logkow is not None else "estimado por RDKit",
    }


def risk_quotients(predictions: pd.DataFrame, occurrence_ug_l: object) -> pd.DataFrame:
    """Calcula un RQ de cribado por grupo, sin mezclar endpoints."""
    occurrence = as_number(occurrence_ug_l)
    result = predictions[[
        "grupo", "endpoint", "concentracion_efecto_mg_l", "dentro_dominio",
        "validacion_predictiva_aceptable",
    ]].copy()
    result["exposicion_local_mg_l"] = np.nan
    result["rq_cribado"] = np.nan
    result["categoria_rq"] = "Sin dato local"
    if occurrence is None or occurrence < 0:
        return result
    result["exposicion_local_mg_l"] = occurrence / 1000
    result["rq_cribado"] = result["exposicion_local_mg_l"] / result["concentracion_efecto_mg_l"]
    result["categoria_rq"] = pd.cut(
        result["rq_cribado"], bins=[-np.inf, 0.1, 1, np.inf],
        labels=["Bajo (cribado)", "Intermedio (cribado)", "Alto (cribado)"],
    ).astype(object)
    result.loc[~result["dentro_dominio"] | ~result["validacion_predictiva_aceptable"], "categoria_rq"] = "No confiable: fuera de dominio/validación"
    return result


def evaluate_row(model_card: dict[str, object], source: dict[str, object]) -> dict[str, object]:
    result = dict(source)
    try:
        predictions, descriptors = predict_by_group(model_card, str(source.get("smiles", "")).strip())
    except ValueError as error:
        result.update({"estado": str(error), "prioridad_general": "No calculable"})
        return result
    general = general_ecotoxicity(predictions)
    result.update({
        "estado": general["status"], "grupo_mas_sensible": general["group"],
        "concentracion_efecto_general_mg_l": general["concentration_mg_l"], "criterio_integracion": general["reason"],
    })
    for item in predictions.to_dict("records"):
        group = str(item["grupo"]).lower().replace("á", "a").replace("é", "e")
        result[f"concentracion_efecto_{group}_mg_l"] = item["concentracion_efecto_mg_l"]
        result[f"dominio_{group}"] = item["dentro_dominio"]
    rq = risk_quotients(predictions, source.get("occurrence_ug_l", ""))
    for item in rq.to_dict("records"):
        group = str(item["grupo"]).lower().replace("á", "a").replace("é", "e")
        result[f"rq_{group}"] = item["rq_cribado"]
        result[f"categoria_rq_{group}"] = item["categoria_rq"]
    if pd.notna(general["concentration_mg_l"]):
        priority = environmental_priority(float(general["concentration_mg_l"]), descriptors, source)
        result.update({
            "prioridad_general": priority["priority"], "puntuacion_prioridad": priority["score"],
            "completitud_ambiental_pct": priority["completeness"] * 100, "logkow_usado": priority["logkow"],
        })
    else:
        result["prioridad_general"] = "No calculable"
    return result


st.set_page_config(page_title="Priorización ecotoxicológica", page_icon="🧪", layout="centered")
st.title("Priorización ecotoxicológica multi-organismo")
st.caption("Apoyo a la priorización de contaminantes emergentes domésticos en Colombia: peligro ecotoxicológico estimado + evidencia ambiental local.")

if not MODEL_PATH.exists():
    st.error("No se encontró el modelo multi-organismo. Ejecuta src/train_multigroup_qsar.py.")
    st.stop()

with st.form("substance"):
    compound_name = st.text_input("Nombre de la sustancia", placeholder="Ej. Paracetamol")
    chemical_key = st.text_input("CAS o identificador", placeholder="Ej. 103-90-2")
    smiles = st.text_area("SMILES *", placeholder="Ej. CC(=O)NC1=CC=C(O)C=C1")
    st.markdown("#### Datos ambientales de Colombia (opcionales, pero necesarios para una prioridad local completa)")
    logkow = st.text_input("logKow medido", placeholder="Ej. 0.5")
    occurrence = st.text_input("Concentración observada local (µg/L)", placeholder="Ej. 0.2")
    matrix = st.text_input("Matriz local", placeholder="Ej. Efluente de PTAR o río")
    location = st.text_input("Ciudad/cuenca y fuente", placeholder="Ej. Cali, río Cauca, artículo 2018")
    removal = st.text_input("Remoción en tratamiento (%)", placeholder="Ej. 85")
    submitted = st.form_submit_button("Predecir y priorizar")

if submitted:
    source = {"compound_name": compound_name, "chemical_key": chemical_key, "smiles": smiles,
              "logkow": logkow, "occurrence_ug_l": occurrence, "matrix": matrix, "location": location,
              "removal_pct": removal}
    removal_value = as_number(removal)
    if removal_value is not None and not 0 <= removal_value <= 100:
        st.error("La remoción debe estar entre 0 y 100 %.")
        st.stop()
    bundle = load_model(MODEL_PATH.stat().st_mtime)
    try:
        predictions, descriptors = predict_by_group(bundle, smiles)
    except ValueError:
        st.error("El SMILES no es válido. Verifica la estructura química.")
        st.stop()
    general = general_ecotoxicity(predictions)
    st.write(f"**Sustancia:** {compound_name or 'Sin nombre'} · **Identificador:** {chemical_key or 'No informado'}")
    if pd.notna(general["concentration_mg_l"]):
        priority = environmental_priority(float(general["concentration_mg_l"]), descriptors, source)
        first, second, third = st.columns(3)
        first.metric("Concentración de efecto general", f"{general['concentration_mg_l']:.3g} mg/L")
        second.metric("Grupo más sensible", str(general["group"]))
        third.metric("Prioridad ambiental", str(priority["priority"]))
        st.write(f"Puntuación: **{priority['score']:.2f} / 3** · logKow usado: **{priority['logkow']:.2f}** ({priority['logkow_source']}).")
        st.dataframe(pd.DataFrame([priority["scores"]]), hide_index=True, width="stretch")
        if priority["completeness"] < 1:
            st.info("Prioridad preliminar: añade ocurrencia en Colombia y remoción de la PTAR o proceso de interés para completar la evaluación local.")
    else:
        st.warning(general["reason"])

    display = predictions.copy()
    display["concentracion_efecto_mg_l"] = display["concentracion_efecto_mg_l"].map(lambda value: f"{value:.3g}")
    display["similitud_tanimoto"] = display["similitud_tanimoto"].map(lambda value: f"{value:.2f}")
    display["r2_validacion"] = display["r2_validacion"].map(lambda value: f"{value:.2f}" if pd.notna(value) else "")
    st.subheader("Resultados por grupo acuático")
    st.dataframe(display[["grupo", "endpoint", "concentracion_efecto_mg_l", "r2_validacion", "similitud_tanimoto", "dentro_dominio", "validacion_predictiva_aceptable"]], hide_index=True, width="stretch")
    st.caption(general["reason"])
    if not predictions["dentro_dominio"].all():
        st.warning("Hay al menos un grupo fuera de su dominio de aplicabilidad; interpreta esa estimación con cautela.")

    rq = risk_quotients(predictions, occurrence)
    if as_number(occurrence) is not None and as_number(occurrence) >= 0:
        st.subheader("Riesgo ecotoxicológico de cribado")
        st.caption(f"RQ = concentración local / concentración de efecto predicha. Matriz: {matrix or 'no informada'} · Ubicación/fuente: {location or 'no informada'}.")
        rq_display = rq.copy()
        rq_display["exposicion_local_mg_l"] = rq_display["exposicion_local_mg_l"].map(lambda value: f"{value:.3g}")
        rq_display["concentracion_efecto_mg_l"] = rq_display["concentracion_efecto_mg_l"].map(lambda value: f"{value:.3g}")
        rq_display["rq_cribado"] = rq_display["rq_cribado"].map(lambda value: f"{value:.3g}")
        st.dataframe(rq_display, hide_index=True, width="stretch")
        st.warning("Esta categoría es de cribado. Solo es interpretable como evidencia local si la concentración corresponde a la matriz y ubicación declaradas; no es una clasificación regulatoria definitiva.")
    else:
        st.info("Para calcular riesgo ecotoxicológico de cribado, añade una concentración ambiental colombiana en µg/L y su matriz/fuente.")

st.divider()
st.subheader("Procesar un archivo CSV")
st.caption("Obligatorias: `compound_name`, `chemical_key`, `smiles`. Para RQ local añade `occurrence_ug_l`, `matrix` y `location`. Opcionales: `logkow`, `removal_pct`.")
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
            required = {"compound_name", "chemical_key", "smiles"}
            missing = required - set(batch.columns)
            if missing:
                st.error(f"Faltan columnas obligatorias: {', '.join(sorted(missing))}.")
            elif len(batch) == 0:
                st.warning("El archivo no contiene filas.")
            elif len(batch) > 1_000:
                st.error("El archivo supera el máximo de 1.000 filas por ejecución.")
            elif st.button("Procesar archivo"):
                bundle = load_model(MODEL_PATH.stat().st_mtime)
                output = pd.DataFrame([evaluate_row(bundle, row) for row in batch.to_dict("records")])
                st.success(f"Se procesaron {len(output)} sustancias.")
                st.dataframe(output, hide_index=True, width="stretch")
                st.download_button("Descargar resultados CSV", data=output.to_csv(index=False).encode("utf-8-sig"),
                                   file_name="priorizacion_ecotoxicologica_multi_organismo.csv", mime="text/csv")

with st.expander("Alcance y limitaciones"):
    st.markdown("""
    - ECOTOX aporta evidencia de peligro ecotoxicológico experimental; no representa por sí solo las concentraciones, fuentes domésticas ni procesos de tratamiento de Colombia.
    - Se entrenó un modelo independiente para cada grupo con ensayos comparables de agua dulce.
    - La integración general usa el grupo con menor concentración de efecto solo si su modelo superó la validación y la sustancia está dentro del dominio de aplicabilidad.
    - Las estimaciones de algas se muestran como evidencia complementaria: su desempeño actual no permite usarlas para decidir la prioridad general.
    - El RQ de cribado se calcula por grupo como exposición local / concentración de efecto predicha. No se mezclan EC50 y LC50 en un único endpoint.
    - La prioridad para Colombia integra efecto ecotoxicológico, logKow y, cuando se informan, ocurrencia local y remoción de la PTAR o proceso de interés. Sin datos locales de exposición, solo se muestra peligro QSAR, no riesgo local.
    """)
