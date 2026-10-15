"""Entrena y evalúa un QSAR de regresión para log10(EC50 mg/L)."""

from __future__ import annotations

import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdFingerprintGenerator, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.ensemble import ExtraTreesRegressor, RandomForestRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GroupKFold


PROJECT = Path(__file__).resolve().parents[1]
INPUT = PROJECT / "data" / "processed" / "qsar_training_ready.csv"
MODELS = PROJECT / "models"
RESULTS = PROJECT / "results"
MIN_UNIQUE_CHEMICALS = 300
RANDOM_STATE = 42
FINGERPRINT_SIZE = 2048
FINGERPRINT_GENERATOR = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=FINGERPRINT_SIZE)


def descriptors_from_smiles(smiles: str) -> dict[str, float] | None:
    """Calcula descriptores 2D interpretables; no usa la toxicidad como entrada."""
    molecule = Chem.MolFromSmiles(str(smiles))
    if molecule is None:
        return None
    return {
        "mol_weight": Descriptors.MolWt(molecule),
        "log_p": Crippen.MolLogP(molecule),
        "tpsa": rdMolDescriptors.CalcTPSA(molecule),
        "h_bond_donors": Lipinski.NumHDonors(molecule),
        "h_bond_acceptors": Lipinski.NumHAcceptors(molecule),
        "rotatable_bonds": Lipinski.NumRotatableBonds(molecule),
        "heavy_atoms": Lipinski.HeavyAtomCount(molecule),
        "rings": Lipinski.RingCount(molecule),
        "aromatic_rings": Lipinski.NumAromaticRings(molecule),
        "fraction_csp3": rdMolDescriptors.CalcFractionCSP3(molecule),
    }


def fingerprint_from_smiles(smiles: str) -> np.ndarray | None:
    """Genera un fingerprint Morgan para evaluar similitud estructural."""
    molecule = Chem.MolFromSmiles(str(smiles))
    if molecule is None:
        return None
    fingerprint = FINGERPRINT_GENERATOR.GetFingerprint(molecule)
    values = np.zeros((FINGERPRINT_SIZE,), dtype=np.int16)
    DataStructs.ConvertToNumpyArray(fingerprint, values)
    return values


def nearest_training_similarity(smiles: str, reference_fingerprints: np.ndarray) -> float:
    """Calcula la mayor similitud de Tanimoto frente al entrenamiento."""
    fingerprint = fingerprint_from_smiles(smiles)
    if fingerprint is None or len(reference_fingerprints) == 0:
        return float("nan")
    intersection = reference_fingerprints @ fingerprint
    union = reference_fingerprints.sum(axis=1) + fingerprint.sum() - intersection
    similarities = np.divide(intersection, union, out=np.zeros_like(intersection, dtype=float), where=union > 0)
    return float(similarities.max())


def prediction_frame_from_smiles(smiles: str, bundle: dict[str, object]) -> tuple[pd.DataFrame, dict[str, float]]:
    """Construye las variables correctas para un modelo de descriptores o híbrido."""
    descriptor_values = descriptors_from_smiles(smiles)
    if descriptor_values is None:
        raise ValueError("El SMILES no puede ser interpretado por RDKit.")
    descriptor_features = bundle.get("descriptor_features", bundle["features"])
    values: dict[str, float] = {feature: descriptor_values[feature] for feature in descriptor_features}
    if bundle.get("feature_mode", "descriptors") == "hybrid":
        fingerprint = fingerprint_from_smiles(smiles)
        if fingerprint is None:
            raise ValueError("No se pudo generar el fingerprint molecular.")
        for index, value in enumerate(fingerprint):
            values[f"morgan_{index}"] = int(value)
    return pd.DataFrame([{feature: values[feature] for feature in bundle["features"]}]), descriptor_values


def scaffold_group(smiles: str, fallback: str) -> str:
    molecule = Chem.MolFromSmiles(str(smiles))
    if molecule is None:
        return fallback
    scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=molecule)
    return scaffold or fallback


def prepare_features(data: pd.DataFrame) -> tuple[pd.DataFrame, list[str], pd.Series]:
    required = {"chemical_key", "smiles", "log_ec50_mg_l"}
    missing = required - set(data.columns)
    if missing:
        raise ValueError(f"Faltan columnas para entrenar: {sorted(missing)}")

    rows: list[dict[str, object]] = []
    for _, row in data.iterrows():
        descriptor_values = descriptors_from_smiles(row["smiles"])
        if descriptor_values is None or pd.isna(row["log_ec50_mg_l"]):
            continue
        descriptor_values.update(
            {
                "chemical_key": str(row["chemical_key"]),
                "smiles": str(row["smiles"]),
                "log_ec50_mg_l": float(row["log_ec50_mg_l"]),
                "scaffold": scaffold_group(str(row["smiles"]), str(row["chemical_key"])),
            }
        )
        rows.append(descriptor_values)
    prepared = pd.DataFrame(rows)
    feature_names = [column for column in prepared.columns if column not in {"chemical_key", "smiles", "log_ec50_mg_l", "scaffold"}]
    return prepared, feature_names, prepared["scaffold"]


def cross_validate(model: object, x: pd.DataFrame, y: pd.Series, groups: pd.Series) -> tuple[pd.DataFrame, pd.DataFrame]:
    folds = min(5, int(groups.nunique()))
    if folds < 3:
        raise ValueError("Se requieren al menos tres familias estructurales para validar el modelo.")
    rows: list[dict[str, float]] = []
    predictions: list[pd.DataFrame] = []
    splitter = GroupKFold(n_splits=folds)
    for fold, (train_index, test_index) in enumerate(splitter.split(x, y, groups), start=1):
        model.fit(x.iloc[train_index], y.iloc[train_index])
        predicted = model.predict(x.iloc[test_index])
        rows.append(
            {
                "fold": fold,
                "mae": mean_absolute_error(y.iloc[test_index], predicted),
                "rmse": mean_squared_error(y.iloc[test_index], predicted) ** 0.5,
                "r2": r2_score(y.iloc[test_index], predicted),
                "test_compounds": len(test_index),
            }
        )
        predictions.append(
            pd.DataFrame(
                {
                    "row_index": test_index,
                    "fold": fold,
                    "observed_log_ec50_mg_l": y.iloc[test_index].to_numpy(),
                    "predicted_log_ec50_mg_l": predicted,
                }
            )
        )
    return pd.DataFrame(rows), pd.concat(predictions, ignore_index=True)


def main() -> None:
    if not INPUT.exists():
        raise FileNotFoundError(f"No existe {INPUT}. Prepara primero los registros de ECOTOX.")
    raw = pd.read_csv(INPUT)
    unique_chemicals = int(raw.get("chemical_key", pd.Series(dtype=str)).nunique())
    if unique_chemicals < MIN_UNIQUE_CHEMICALS:
        raise ValueError(
            f"No se entrenó: hay {unique_chemicals} sustancias únicas y se requieren al menos "
            f"{MIN_UNIQUE_CHEMICALS} con endpoint y estructura verificables."
        )

    prepared, feature_names, groups = prepare_features(raw)
    if len(prepared) < MIN_UNIQUE_CHEMICALS:
        raise ValueError("Hay menos de 300 estructuras SMILES válidas después de la depuración.")
    x = prepared[feature_names]
    y = prepared["log_ec50_mg_l"]
    candidates = {
        "random_forest": RandomForestRegressor(n_estimators=500, min_samples_leaf=2, random_state=RANDOM_STATE, n_jobs=-1),
        "extra_trees": ExtraTreesRegressor(n_estimators=500, min_samples_leaf=2, random_state=RANDOM_STATE, n_jobs=-1),
    }

    scores: dict[str, pd.DataFrame] = {}
    out_of_fold: dict[str, pd.DataFrame] = {}
    for name, model in candidates.items():
        scores[name], out_of_fold[name] = cross_validate(model, x, y, groups)
    summary = pd.DataFrame(
        [
            {"model": name, "mae_mean": frame["mae"].mean(), "rmse_mean": frame["rmse"].mean(), "r2_mean": frame["r2"].mean()}
            for name, frame in scores.items()
        ]
    ).sort_values("mae_mean")

    selected_name = str(summary.iloc[0]["model"])
    selected = candidates[selected_name]
    selected.fit(x, y)
    fingerprints = np.stack([fingerprint_from_smiles(smiles) for smiles in prepared["smiles"]])
    intersections = fingerprints @ fingerprints.T
    bit_counts = fingerprints.sum(axis=1)
    unions = bit_counts[:, None] + bit_counts[None, :] - intersections
    similarities = np.divide(intersections, unions, out=np.zeros_like(intersections, dtype=float), where=unions > 0)
    np.fill_diagonal(similarities, -np.inf)
    nearest_training = similarities.max(axis=1)
    similarity_threshold = float(np.quantile(nearest_training, 0.05))
    MODELS.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    joblib.dump(
        {
            "model": selected,
            "features": feature_names,
            "descriptor_features": feature_names,
            "feature_mode": "descriptors",
            "feature_ranges": {feature: {"min": float(x[feature].min()), "max": float(x[feature].max())} for feature in feature_names},
            "reference_fingerprints": fingerprints,
            "applicability_similarity_threshold": similarity_threshold,
            "endpoint": "Daphnia magna EC50 48 h inmovilización",
        },
        MODELS / "qsar_daphnia_ec50.joblib",
    )
    summary.to_csv(RESULTS / "model_comparison.csv", index=False)
    scores[selected_name].to_csv(RESULTS / "selected_model_cv_folds.csv", index=False)
    selected_oof = out_of_fold[selected_name].copy().sort_values("row_index")
    selected_oof["chemical_key"] = prepared.iloc[selected_oof["row_index"]]["chemical_key"].to_numpy()
    selected_oof["scaffold"] = prepared.iloc[selected_oof["row_index"]]["scaffold"].to_numpy()
    selected_oof["residual_log_ec50_mg_l"] = selected_oof["observed_log_ec50_mg_l"] - selected_oof["predicted_log_ec50_mg_l"]
    selected_oof.drop(columns="row_index").to_csv(RESULTS / "out_of_fold_predictions.csv", index=False)
    importance = pd.DataFrame(
        {"variable": feature_names, "importance": selected.feature_importances_}
    ).sort_values("importance", ascending=False)
    importance.to_csv(RESULTS / "feature_importance.csv", index=False)
    (RESULTS / "training_report.json").write_text(
        json.dumps(
            {
                "selected_model": selected_name,
                "unique_chemicals": int(len(prepared)),
                "structural_families": int(groups.nunique()),
                "features": feature_names,
                "metrics": summary.to_dict(orient="records"),
                "feature_importance": importance.to_dict(orient="records"),
                "applicability_domain": {
                    "fingerprint": "Morgan radius 2, 2048 bits",
                    "minimum_nearest_neighbor_tanimoto": similarity_threshold,
                    "training_nearest_neighbor_tanimoto_median": float(np.median(nearest_training)),
                },
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    print(f"Modelo seleccionado: {selected_name}")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
