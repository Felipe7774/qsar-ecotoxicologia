"""Compara descriptores, fingerprints Morgan y una combinación de ambos."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from sklearn.ensemble import ExtraTreesRegressor

from train_qsar import (
    FINGERPRINT_SIZE,
    RANDOM_STATE,
    cross_validate,
    descriptors_from_smiles,
    fingerprint_from_smiles,
    scaffold_group,
)


PROJECT = Path(__file__).resolve().parents[1]
INPUT = PROJECT / "data" / "processed" / "qsar_training_ready.csv"
RESULTS = PROJECT / "results"


def main() -> None:
    if not INPUT.exists():
        raise FileNotFoundError("No existe el conjunto de entrenamiento final.")
    data = pd.read_csv(INPUT)
    descriptor_rows: list[dict[str, object]] = []
    fingerprints = []
    groups = []
    targets = []
    for _, row in data.iterrows():
        descriptor_values = descriptors_from_smiles(str(row["smiles"]))
        fingerprint = fingerprint_from_smiles(str(row["smiles"]))
        if descriptor_values is None or fingerprint is None:
            continue
        descriptor_rows.append(descriptor_values)
        fingerprints.append(fingerprint)
        groups.append(scaffold_group(str(row["smiles"]), str(row["chemical_key"])))
        targets.append(float(row["log_ec50_mg_l"]))

    x_descriptors = pd.DataFrame(descriptor_rows)
    x_fingerprints = pd.DataFrame(fingerprints, columns=[f"morgan_{index}" for index in range(FINGERPRINT_SIZE)])
    x_hybrid = pd.concat([x_descriptors.reset_index(drop=True), x_fingerprints.reset_index(drop=True)], axis=1)
    y = pd.Series(targets)
    group_series = pd.Series(groups)
    feature_sets = {
        "descriptores_2D": x_descriptors,
        "morgan_2048": x_fingerprints,
        "hibrido": x_hybrid,
    }

    rows = []
    for name, features in feature_sets.items():
        model = ExtraTreesRegressor(
            n_estimators=500, min_samples_leaf=2, random_state=RANDOM_STATE, n_jobs=-1
        )
        folds, _ = cross_validate(model, features, y, group_series)
        rows.append(
            {
                "representacion": name,
                "n_variables": int(features.shape[1]),
                "mae_mean": float(folds["mae"].mean()),
                "rmse_mean": float(folds["rmse"].mean()),
                "r2_mean": float(folds["r2"].mean()),
            }
        )
    benchmark = pd.DataFrame(rows).sort_values("mae_mean")
    RESULTS.mkdir(parents=True, exist_ok=True)
    benchmark.to_csv(RESULTS / "fingerprint_model_benchmark.csv", index=False)
    print(benchmark.to_string(index=False))


if __name__ == "__main__":
    main()
