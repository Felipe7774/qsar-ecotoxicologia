"""Aplica control de calidad y consolida el objetivo QSAR desde ECOTOX."""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
INPUT = PROJECT / "data" / "processed" / "ecotox_daphnia_ec50_48h_freshwater_raw.csv"
PROCESSED = PROJECT / "data" / "processed"
RESULTS = PROJECT / "results"


def casrn(value: object) -> str:
    """Restaura guiones en un CAS numérico de ECOTOX sin validar su dígito."""
    digits = re.sub(r"\D", "", str(value))
    if len(digits) < 5:
        return ""
    return f"{digits[:-3]}-{digits[-3:-1]}-{digits[-1]}"


def to_mg_l(value: float, unit: object) -> float:
    """Convierte unidades de concentración en agua aceptadas y documentadas."""
    if pd.isna(value) or value <= 0:
        return np.nan
    text = str(unit).strip().casefold()
    factors = {
        "mg/l": 1.0,
        "ai mg/l": 1.0,
        "ug/l": 1 / 1_000,
        "ai ug/l": 1 / 1_000,
        "ae ug/l": 1 / 1_000,
        "ng/l": 1 / 1_000_000,
        "ai ng/l": 1 / 1_000_000,
        "ppm": 1.0,
        "ppb": 1 / 1_000,
        "ppt": 1 / 1_000_000,
    }
    return value * factors[text] if text in factors else np.nan


def main() -> None:
    if not INPUT.exists():
        raise FileNotFoundError(f"No existe {INPUT}. Ejecuta extract_ecotox_endpoint.py primero.")
    PROCESSED.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    raw = pd.read_csv(INPUT, dtype=str)
    raw["ec50_value"] = pd.to_numeric(raw["ec50_value"], errors="coerce")
    raw["casrn"] = raw["test_cas"].map(casrn)
    raw["ec50_mg_l"] = [to_mg_l(value, unit) for value, unit in zip(raw["ec50_value"], raw["conc1_unit"])]

    valid_measurement = raw["measurement"].fillna("").str.startswith("IMBL")
    valid_operator = raw["conc1_mean_op"].fillna("").isin(["", "="])
    valid_identity = raw["casrn"].ne("")
    valid_value = raw["ec50_mg_l"].notna() & raw["ec50_mg_l"].gt(0)
    qc = raw.loc[valid_measurement & valid_operator & valid_identity & valid_value].copy()
    qc["log_ec50_mg_l"] = np.log10(qc["ec50_mg_l"])
    qc.to_csv(PROCESSED / "ecotox_daphnia_ec50_48h_immobilization_qc.csv", index=False, encoding="utf-8-sig")

    training = (
        qc.groupby("casrn", as_index=False)
        .agg(
            chemical_key=("casrn", "first"),
            chemical_name=("chemical_name", "first"),
            dtxsid=("dtxsid", "first"),
            log_ec50_mg_l=("log_ec50_mg_l", "median"),
            ec50_mg_l=("ec50_mg_l", "median"),
            n_assays=("result_id", "count"),
            n_references=("reference_number", "nunique"),
            log_ec50_sd=("log_ec50_mg_l", "std"),
        )
    )
    training["smiles"] = ""
    training.to_csv(PROCESSED / "qsar_training_before_structures.csv", index=False, encoding="utf-8-sig")

    report = {
        "raw_ec50_records": int(len(raw)),
        "immobilization_measurement_records": int(valid_measurement.sum()),
        "excluded_non_exact_operators": int((~valid_operator).sum()),
        "accepted_qc_records": int(len(qc)),
        "unique_substances_before_structure_resolution": int(len(training)),
        "median_assays_per_substance": float(training["n_assays"].median()) if not training.empty else 0,
        "training_table": str(PROCESSED / "qsar_training_before_structures.csv"),
        "note": "ppm y ppb se trataron como concentración acuosa aproximada mg/L y ug/L; las unidades molares se excluyeron.",
    }
    (RESULTS / "qsar_training_qc_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
