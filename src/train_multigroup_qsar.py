"""Entrena un conjunto QSAR por grupo de organismos acuáticos.

Los modelos se entrenan por separado para no tratar como equivalentes la
inmovilización de crustáceos, la mortalidad de peces y la respuesta de algas.
La capa de priorización usa después el resultado más sensible de los grupos.
"""

from __future__ import annotations

import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import ExtraTreesRegressor

from train_hybrid_qsar import applicability_domain
from train_qsar import FINGERPRINT_SIZE, RANDOM_STATE, cross_validate, fingerprint_from_smiles, prepare_features


PROJECT = Path(__file__).resolve().parents[1]
INPUT = PROJECT / "data" / "processed" / "multigrupo_training_with_structures.csv"
MODEL_PATH = PROJECT / "models" / "qsar_multiorganismo.joblib"
RESULTS = PROJECT / "results"
MIN_CHEMICALS_PER_GROUP = 100


def train_group(name: str, source: pd.DataFrame) -> tuple[dict[str, object], dict[str, object], pd.DataFrame, pd.DataFrame]:
    prepared, descriptor_features, scaffolds = prepare_features(source)
    if len(prepared) < MIN_CHEMICALS_PER_GROUP:
        raise ValueError(f"{name}: se requieren al menos {MIN_CHEMICALS_PER_GROUP} sustancias con estructura válida.")
    fingerprints = np.stack([fingerprint_from_smiles(smiles) for smiles in prepared["smiles"]])
    fingerprint_columns = [f"morgan_{index}" for index in range(FINGERPRINT_SIZE)]
    x = pd.concat(
        [prepared[descriptor_features].reset_index(drop=True), pd.DataFrame(fingerprints, columns=fingerprint_columns)], axis=1
    )
    y = prepared["log_ec50_mg_l"].reset_index(drop=True)
    model = ExtraTreesRegressor(n_estimators=400, min_samples_leaf=2, random_state=RANDOM_STATE, n_jobs=-1)
    folds, oof = cross_validate(model, x, y, scaffolds.reset_index(drop=True))
    model.fit(x, y)
    nearest, threshold = applicability_domain(fingerprints)
    bundle = {
        "model": model,
        "features": list(x.columns),
        "descriptor_features": descriptor_features,
        "feature_mode": "hybrid",
        "feature_ranges": {
            feature: {"min": float(x[feature].min()), "max": float(x[feature].max())}
            for feature in descriptor_features
        },
        "reference_fingerprints": fingerprints,
        "applicability_similarity_threshold": threshold,
        "endpoint": str(source["endpoint_descripcion"].iloc[0]),
    }
    metrics = {
        "group": name,
        "unique_chemicals": int(len(prepared)),
        "structural_families": int(scaffolds.nunique()),
        "mae_mean": float(folds["mae"].mean()),
        "rmse_mean": float(folds["rmse"].mean()),
        "r2_mean": float(folds["r2"].mean()),
        "applicability_similarity_threshold": float(threshold),
        "median_training_nearest_tanimoto": float(np.median(nearest)),
        "endpoint": bundle["endpoint"],
    }
    oof = oof.sort_values("row_index")
    oof["grupo_organismos"] = name
    oof["chemical_key"] = prepared.iloc[oof["row_index"]]["chemical_key"].to_numpy()
    oof["scaffold"] = prepared.iloc[oof["row_index"]]["scaffold"].to_numpy()
    oof["residual_log_mg_l"] = oof["observed_log_ec50_mg_l"] - oof["predicted_log_ec50_mg_l"]
    return bundle, metrics, folds.assign(grupo_organismos=name), oof.drop(columns="row_index")


def main() -> None:
    if not INPUT.exists():
        raise FileNotFoundError(f"No existe {INPUT}. Ejecuta la resolución de estructuras primero.")
    raw = pd.read_csv(INPUT, dtype={"chemical_key": str, "smiles": str})
    ready = raw.loc[
        raw["structure_status"].eq("rdkit_valid")
        & raw["smiles"].notna()
        & raw["log_ec50_mg_l"].notna()
    ].copy()
    if ready.empty:
        raise ValueError("No hay filas multi-organismo con estructura química válida.")

    bundles: dict[str, dict[str, object]] = {}
    summaries: list[dict[str, object]] = []
    all_folds: list[pd.DataFrame] = []
    all_oof: list[pd.DataFrame] = []
    for group, source in ready.groupby("grupo_organismos", sort=True):
        bundle, summary, folds, oof = train_group(str(group), source.reset_index(drop=True))
        bundles[str(group)] = bundle
        summaries.append(summary)
        all_folds.append(folds)
        all_oof.append(oof)

    MODEL_PATH.parent.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    model_card = {
        "model_type": "multi_group_qsar_extra_trees_hybrid",
        "integration_rule": "La prioridad ecotoxicológica general usa la concentración estimada más baja entre los grupos con validación predictiva aceptable y dentro de dominio de aplicabilidad.",
        "aggregation_groups": [row["group"] for row in summaries if row["r2_mean"] > 0],
        "group_metrics": {row["group"]: row for row in summaries},
        "groups": bundles,
    }
    joblib.dump(model_card, MODEL_PATH)
    summary = pd.DataFrame(summaries).sort_values("group")
    summary.to_csv(RESULTS / "multigrupo_model_metrics.csv", index=False, encoding="utf-8-sig")
    pd.concat(all_folds, ignore_index=True).to_csv(RESULTS / "multigrupo_cv_folds.csv", index=False, encoding="utf-8-sig")
    pd.concat(all_oof, ignore_index=True).to_csv(RESULTS / "multigrupo_out_of_fold_predictions.csv", index=False, encoding="utf-8-sig")
    report = {
        "model": "Extra Trees híbrido: 10 descriptores 2D + fingerprint Morgan de 2048 bits",
        "training_rows": int(len(ready)),
        "groups": summaries,
        "integration_rule": model_card["integration_rule"],
        "model_path": str(MODEL_PATH),
    }
    (RESULTS / "multigrupo_training_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
