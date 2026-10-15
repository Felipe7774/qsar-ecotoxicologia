"""Genera una predicción QSAR para un SMILES con revisión básica de aplicabilidad."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import joblib
import pandas as pd

try:
    from .train_qsar import nearest_training_similarity, prediction_frame_from_smiles
except ImportError:  # Ejecución directa: python src/predict_qsar.py
    from train_qsar import nearest_training_similarity, prediction_frame_from_smiles


PROJECT = Path(__file__).resolve().parents[1]
MODEL_PATH = PROJECT / "models" / "qsar_daphnia_ec50.joblib"


def main() -> None:
    parser = argparse.ArgumentParser(description="Predice EC50 de Daphnia magna a 48 h a partir de un SMILES.")
    parser.add_argument("smiles", help="SMILES de la sustancia individual")
    args = parser.parse_args()
    if not MODEL_PATH.exists():
        raise FileNotFoundError("No existe el modelo entrenado. Ejecuta src/train_qsar.py primero.")

    bundle = joblib.load(MODEL_PATH)
    frame, descriptor_values = prediction_frame_from_smiles(args.smiles, bundle)
    log_prediction = float(bundle["model"].predict(frame)[0])
    outside_ranges = [
        feature for feature, value in descriptor_values.items()
        if value < bundle["feature_ranges"][feature]["min"] or value > bundle["feature_ranges"][feature]["max"]
    ]
    similarity = nearest_training_similarity(args.smiles, bundle["reference_fingerprints"])
    similarity_warning = bool(similarity < bundle["applicability_similarity_threshold"])
    result = {
        "endpoint": bundle["endpoint"],
        "predicted_log_ec50_mg_l": log_prediction,
        "predicted_ec50_mg_l": 10 ** log_prediction,
        "applicability_warning": bool(outside_ranges) or similarity_warning,
        "descriptors_outside_training_range": outside_ranges,
        "nearest_training_tanimoto": similarity,
        "minimum_accepted_tanimoto": bundle["applicability_similarity_threshold"],
        "similarity_warning": similarity_warning,
    }
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
