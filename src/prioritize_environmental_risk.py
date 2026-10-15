"""Convierte predicción QSAR y datos ambientales en una prioridad explicable.

La prioridad es una regla de apoyo a decisión, no una nueva etiqueta con la que
se entrena el modelo QSAR. Los campos ambientales faltantes disminuyen la
completitud; nunca se sustituyen por cero.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

try:
    from .train_qsar import nearest_training_similarity, prediction_frame_from_smiles
except ImportError:  # Ejecución directa: python src/prioritize_environmental_risk.py
    from train_qsar import nearest_training_similarity, prediction_frame_from_smiles


PROJECT = Path(__file__).resolve().parents[1]
MODEL_PATH = PROJECT / "models" / "qsar_daphnia_ec50.joblib"
DEFAULT_INPUT = PROJECT / "data" / "processed" / "priority_input.csv"
DEFAULT_OUTPUT = PROJECT / "results" / "environmental_priority.csv"
REQUIRED_COLUMNS = {"compound_name", "chemical_key", "smiles"}
OPTIONAL_COLUMNS = {"logkow", "occurrence_ug_l", "removal_pct"}


def as_number(value: object) -> float:
    number = pd.to_numeric(value, errors="coerce")
    return float(number) if pd.notna(number) else np.nan


def ecotox_score(ec50_mg_l: float) -> float:
    if ec50_mg_l <= 1:
        return 3.0
    if ec50_mg_l <= 10:
        return 2.0
    return 1.0


def logkow_score(value: float) -> float:
    if value >= 4:
        return 3.0
    if value >= 2:
        return 2.0
    return 1.0


def occurrence_score(value: float) -> float:
    if value >= 1:
        return 3.0
    if value >= 0.1:
        return 2.0
    return 1.0


def removal_score(value: float) -> float:
    if value <= 50:
        return 3.0
    if value <= 80:
        return 2.0
    return 1.0


def priority_label(mean_score: float) -> str:
    if mean_score >= 2.34:
        return "Alta"
    if mean_score >= 1.67:
        return "Media"
    return "Baja"


def predict(bundle: dict[str, object], smiles: str) -> tuple[float, float, float, list[str], float]:
    row, values = prediction_frame_from_smiles(smiles, bundle)
    log_ec50 = float(bundle["model"].predict(row)[0])
    outside = [
        feature for feature in bundle["descriptor_features"]
        if values[feature] < bundle["feature_ranges"][feature]["min"]
        or values[feature] > bundle["feature_ranges"][feature]["max"]
    ]
    similarity = nearest_training_similarity(smiles, bundle["reference_fingerprints"])
    return log_ec50, 10 ** log_ec50, float(values["log_p"]), outside, similarity


def evaluate_substance(bundle: dict[str, object], source: dict[str, object]) -> dict[str, object]:
    """Predice y prioriza una sustancia conservando sus columnas de entrada."""
    result = dict(source)
    smiles = str(source.get("smiles", "")).strip()
    if not smiles:
        result.update({"priority_status": "Falta SMILES", "priority": "No calculable"})
        return result

    try:
        log_ec50, predicted_ec50, derived_logp, outside, similarity = predict(bundle, smiles)
    except ValueError as error:
        result.update({"priority_status": str(error), "priority": "No calculable"})
        return result

    measured_logkow = as_number(source.get("logkow"))
    occurrence = as_number(source.get("occurrence_ug_l"))
    removal = as_number(source.get("removal_pct"))
    logkow_used = measured_logkow if pd.notna(measured_logkow) else derived_logp
    score_items = {
        "ecotoxicity_score": ecotox_score(predicted_ec50),
        "logkow_score": logkow_score(logkow_used),
    }
    if pd.notna(occurrence) and occurrence >= 0:
        score_items["occurrence_score"] = occurrence_score(occurrence)
    if pd.notna(removal) and 0 <= removal <= 100:
        score_items["removal_score"] = removal_score(removal)

    mean_score = float(np.mean(list(score_items.values())))
    completeness = len(score_items) / 4
    result.update(
        {
            "predicted_log_ec50_mg_l": log_ec50,
            "predicted_ec50_mg_l": predicted_ec50,
            "logkow_used": logkow_used,
            "logkow_source": "medido" if pd.notna(measured_logkow) else "estimado_RDKit",
            **score_items,
            "priority_score_mean": mean_score,
            "priority": priority_label(mean_score),
            "environmental_data_completeness_pct": completeness * 100,
            "applicability_warning": bool(outside) or similarity < bundle["applicability_similarity_threshold"],
            "descriptors_outside_training_range": "; ".join(outside),
            "nearest_training_tanimoto": similarity,
            "minimum_accepted_tanimoto": bundle["applicability_similarity_threshold"],
            "similarity_warning": similarity < bundle["applicability_similarity_threshold"],
            "priority_status": "Completa" if completeness == 1 else "Preliminar: faltan datos ambientales",
        }
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description="Clasifica prioridad ambiental con la predicción QSAR.")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="CSV de entrada con compound_name, chemical_key y smiles")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help="CSV de resultados")
    args = parser.parse_args()
    if not MODEL_PATH.exists():
        raise FileNotFoundError("No existe el modelo QSAR entrenado.")
    if not args.input.exists():
        raise FileNotFoundError(
            f"No existe {args.input}. Crea un CSV con: compound_name, chemical_key, smiles, "
            "logkow (opcional), occurrence_ug_l (opcional), removal_pct (opcional)."
        )

    source = pd.read_csv(args.input)
    missing = REQUIRED_COLUMNS - set(source.columns)
    if missing:
        raise ValueError(f"Faltan columnas requeridas: {sorted(missing)}")
    for column in OPTIONAL_COLUMNS - set(source.columns):
        source[column] = np.nan

    bundle = joblib.load(MODEL_PATH)
    rows = [evaluate_substance(bundle, row.to_dict()) for _, row in source.iterrows()]

    output = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, index=False, encoding="utf-8-sig")
    print(f"Prioridades calculadas: {len(output)}")
    print(f"Salida: {args.output}")


if __name__ == "__main__":
    main()
