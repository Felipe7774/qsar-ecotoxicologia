"""Publica el conjunto final de entrenamiento y su ficha de trazabilidad."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
INPUT = PROJECT / "data" / "processed" / "qsar_training_with_structures.csv"
OUTPUT = PROJECT / "data" / "processed" / "qsar_training_ready.csv"
RESULTS = PROJECT / "results"


def main() -> None:
    if not INPUT.exists():
        raise FileNotFoundError(f"No existe {INPUT}. Ejecuta resolve_training_structures.py primero.")
    data = pd.read_csv(INPUT, dtype={"chemical_key": str, "casrn": str})
    ready = data.loc[data["structure_status"].eq("rdkit_valid")].copy()
    ready = ready.loc[ready["log_ec50_mg_l"].notna() & ready["smiles"].notna()].copy()
    ready.to_csv(OUTPUT, index=False, encoding="utf-8-sig")
    report = {
        "endpoint": "Daphnia magna; agua dulce; exposición 48 h; EC50; inmovilización",
        "target": "log10(EC50_mg_L)",
        "unique_chemicals": int(ready["chemical_key"].nunique()),
        "eligible_assays_after_qc": int(ready["n_assays"].sum()),
        "records_excluded_without_rdkit_smiles": int(len(data) - len(ready)),
        "min_recommended_unique_chemicals": 300,
        "ready_to_train": bool(ready["chemical_key"].nunique() >= 300),
        "dataset": str(OUTPUT),
    }
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "dataset_card.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
