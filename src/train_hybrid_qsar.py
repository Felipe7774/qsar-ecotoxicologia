"""Entrena la versión QSAR seleccionada: descriptores 2D más fingerprints Morgan."""

from __future__ import annotations

import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import ExtraTreesRegressor

from train_qsar import (
    FINGERPRINT_SIZE,
    RANDOM_STATE,
    cross_validate,
    fingerprint_from_smiles,
    prepare_features,
)


PROJECT = Path(__file__).resolve().parents[1]
INPUT = PROJECT / "data" / "processed" / "qsar_training_ready.csv"
MODELS = PROJECT / "models"
RESULTS = PROJECT / "results"


def applicability_domain(fingerprints: np.ndarray) -> tuple[np.ndarray, float]:
    intersections = fingerprints @ fingerprints.T
    bit_counts = fingerprints.sum(axis=1)
    unions = bit_counts[:, None] + bit_counts[None, :] - intersections
    similarities = np.divide(intersections, unions, out=np.zeros_like(intersections, dtype=float), where=unions > 0)
    np.fill_diagonal(similarities, -np.inf)
    nearest = similarities.max(axis=1)
    return nearest, float(np.quantile(nearest, 0.05))


def main() -> None:
    if not INPUT.exists():
        raise FileNotFoundError("No existe qsar_training_ready.csv.")
    raw = pd.read_csv(INPUT)
    prepared, descriptor_features, groups = prepare_features(raw)
    fingerprints = np.stack([fingerprint_from_smiles(smiles) for smiles in prepared["smiles"]])
    fingerprint_columns = [f"morgan_{index}" for index in range(FINGERPRINT_SIZE)]
    x_fingerprints = pd.DataFrame(fingerprints, columns=fingerprint_columns)
    x = pd.concat([prepared[descriptor_features].reset_index(drop=True), x_fingerprints], axis=1)
    y = prepared["log_ec50_mg_l"].reset_index(drop=True)
    model = ExtraTreesRegressor(n_estimators=500, min_samples_leaf=2, random_state=RANDOM_STATE, n_jobs=-1)
    folds, oof = cross_validate(model, x, y, groups.reset_index(drop=True))
    model.fit(x, y)
    nearest, threshold = applicability_domain(fingerprints)

    MODELS.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    joblib.dump(
        {
            "model": model,
            "features": list(x.columns),
            "descriptor_features": descriptor_features,
            "feature_mode": "hybrid",
            "feature_ranges": {feature: {"min": float(x[feature].min()), "max": float(x[feature].max())} for feature in descriptor_features},
            "reference_fingerprints": fingerprints,
            "applicability_similarity_threshold": threshold,
            "endpoint": "Daphnia magna EC50 48 h inmovilización",
        },
        MODELS / "qsar_daphnia_ec50.joblib",
    )

    oof = oof.sort_values("row_index")
    oof["chemical_key"] = prepared.iloc[oof["row_index"]]["chemical_key"].to_numpy()
    oof["scaffold"] = prepared.iloc[oof["row_index"]]["scaffold"].to_numpy()
    oof["residual_log_ec50_mg_l"] = oof["observed_log_ec50_mg_l"] - oof["predicted_log_ec50_mg_l"]
    oof.drop(columns="row_index").to_csv(RESULTS / "out_of_fold_predictions.csv", index=False)
    folds.to_csv(RESULTS / "selected_model_cv_folds.csv", index=False)
    comparison = pd.DataFrame(
        [{"model": "extra_trees_hybrid", "mae_mean": folds["mae"].mean(), "rmse_mean": folds["rmse"].mean(), "r2_mean": folds["r2"].mean()}]
    )
    comparison.to_csv(RESULTS / "model_comparison.csv", index=False)

    importances = pd.Series(model.feature_importances_, index=x.columns)
    descriptor_importance = (
        importances.loc[descriptor_features].rename_axis("variable").reset_index(name="importance").sort_values("importance", ascending=False)
    )
    descriptor_importance.to_csv(RESULTS / "feature_importance.csv", index=False)
    pd.DataFrame(
        {
            "componente": ["descriptores_2D", "fingerprints_Morgan"],
            "importance_total": [float(importances.loc[descriptor_features].sum()), float(importances.loc[fingerprint_columns].sum())],
        }
    ).to_csv(RESULTS / "representation_importance.csv", index=False)
    report = {
        "selected_model": "extra_trees_hybrid",
        "representation": "10 descriptores 2D + fingerprint Morgan (radio 2, 2048 bits)",
        "unique_chemicals": int(len(prepared)),
        "structural_families": int(groups.nunique()),
        "metrics": comparison.to_dict(orient="records"),
        "applicability_similarity_threshold": threshold,
        "median_training_nearest_tanimoto": float(np.median(nearest)),
    }
    (RESULTS / "training_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(comparison.to_string(index=False))


if __name__ == "__main__":
    main()
