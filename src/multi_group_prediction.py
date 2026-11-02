"""Predicción e integración transparente de los modelos por grupo acuático."""

from __future__ import annotations

import numpy as np
import pandas as pd

from .train_qsar import nearest_training_similarity, prediction_frame_from_smiles


def predict_by_group(model_card: dict[str, object], smiles: str) -> tuple[pd.DataFrame, dict[str, float]]:
    rows: list[dict[str, object]] = []
    descriptor_values: dict[str, float] | None = None
    metrics = model_card.get("group_metrics", {})
    for group, bundle in model_card["groups"].items():
        x, values = prediction_frame_from_smiles(smiles, bundle)
        descriptor_values = values
        log_concentration = float(bundle["model"].predict(x)[0])
        concentration = float(10 ** log_concentration)
        outside = [
            feature for feature in bundle["descriptor_features"]
            if values[feature] < bundle["feature_ranges"][feature]["min"]
            or values[feature] > bundle["feature_ranges"][feature]["max"]
        ]
        similarity = nearest_training_similarity(smiles, bundle["reference_fingerprints"])
        threshold = float(bundle["applicability_similarity_threshold"])
        group_metrics = metrics.get(group, {})
        rows.append(
            {
                "grupo": group,
                "endpoint": bundle["endpoint"],
                "log10_concentracion_efecto_mg_l": log_concentration,
                "concentracion_efecto_mg_l": concentration,
                "mae_validacion": group_metrics.get("mae_mean"),
                "r2_validacion": group_metrics.get("r2_mean"),
                "similitud_tanimoto": similarity,
                "umbral_tanimoto": threshold,
                "dentro_dominio": not outside and similarity >= threshold,
                "descriptores_fuera_rango": "; ".join(outside),
                "validacion_predictiva_aceptable": group in model_card.get("aggregation_groups", []),
            }
        )
    if descriptor_values is None:
        raise ValueError("No hay modelos cargados.")
    return pd.DataFrame(rows), descriptor_values


def general_ecotoxicity(predictions: pd.DataFrame) -> dict[str, object]:
    """Devuelve el grupo más sensible entre los modelos validados y aplicables."""
    usable = predictions.loc[
        predictions["validacion_predictiva_aceptable"] & predictions["dentro_dominio"]
    ].copy()
    if usable.empty:
        return {
            "status": "Sin estimación general confiable",
            "group": None,
            "concentration_mg_l": np.nan,
            "reason": "Ningún grupo con validación aceptable quedó dentro de su dominio de aplicabilidad.",
        }
    most_sensitive = usable.loc[usable["concentracion_efecto_mg_l"].idxmin()]
    return {
        "status": "Estimación general basada en grupos aplicables",
        "group": str(most_sensitive["grupo"]),
        "concentration_mg_l": float(most_sensitive["concentracion_efecto_mg_l"]),
        "reason": "Se seleccionó la concentración de efecto más baja entre los grupos con validación predictiva aceptable.",
    }
