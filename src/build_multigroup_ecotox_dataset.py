"""Construye conjuntos QSAR comparables para tres grupos acuáticos.

No mezcla todos los ensayos en una sola variable objetivo. Cada grupo mantiene
un endpoint agudo, una duración y un efecto explícitos; posteriormente las
predicciones se integran mediante el resultado más sensible.
"""

from __future__ import annotations

import json
import re
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
ZIP_FILE = PROJECT / "data" / "raw" / "ecotox_ascii_09_15_2026.zip"
PROCESSED = PROJECT / "data" / "processed"
RESULTS = PROJECT / "results"
PREFIX = "ecotox_ascii_09_15_2026/"

# Cada especificación constituye un endpoint independiente y documentado.
GROUPS = {
    "Crustáceos": {
        "species_ids": {"5"},  # Daphnia magna
        "endpoint": "EC50",
        "effect": "ITX",
        "measurement_prefixes": ("IMBL",),
        "hours": "48",
        "description": "EC50 de inmovilización, 48 h, agua dulce",
    },
    "Peces": {
        "species_ids": {"4", "5156"},  # trucha arcoíris y pez cebra
        "endpoint": "LC50",
        "effect": "MOR",
        "measurement_prefixes": ("MORT",),
        "hours": "96",
        "description": "LC50 de mortalidad, 96 h, agua dulce",
    },
    "Algas": {
        # Raphidocelis subcapitata y su sinónimo histórico Selenastrum capricornutum.
        "species_ids": {"58471", "58477"},
        "endpoint": "EC50",
        "effect": "POP",
        "measurement_prefixes": ("ABND", "PGRT", "CHLO"),
        "hours": "72",
        "description": "EC50 de respuesta poblacional o crecimiento, 72 h, agua dulce",
    },
}


def strip_frame(frame: pd.DataFrame) -> pd.DataFrame:
    return frame.apply(lambda column: column.str.strip() if column.dtype == "object" else column)


def casrn(value: object) -> str:
    digits = re.sub(r"\D", "", str(value))
    return f"{digits[:-3]}-{digits[-3:-1]}-{digits[-1]}" if len(digits) >= 5 else ""


def to_mg_l(value: float, unit: object) -> float:
    if pd.isna(value) or value <= 0:
        return np.nan
    factors = {
        "mg/l": 1.0,
        "mg/dm3": 1.0,
        "ai mg/l": 1.0,
        "ae mg/l": 1.0,
        "ug/l": 1 / 1_000,
        "ai ug/l": 1 / 1_000,
        "ae ug/l": 1 / 1_000,
        "ng/l": 1 / 1_000_000,
        "ai ng/l": 1 / 1_000_000,
        "ppm": 1.0,
        "ppb": 1 / 1_000,
        "ppt": 1 / 1_000_000,
    }
    return value * factors.get(str(unit).strip().casefold(), np.nan)


def main() -> None:
    if not ZIP_FILE.exists():
        raise FileNotFoundError(f"No existe la descarga ECOTOX: {ZIP_FILE}")
    PROCESSED.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)

    species_to_group = {species: group for group, config in GROUPS.items() for species in config["species_ids"]}
    with zipfile.ZipFile(ZIP_FILE) as archive:
        test_columns = [
            "test_id", "reference_number", "test_cas", "species_number", "organism_habitat",
            "exposure_duration_mean", "exposure_duration_unit", "media_type", "study_type", "test_type",
        ]
        test_parts: list[pd.DataFrame] = []
        for chunk in pd.read_csv(
            archive.open(PREFIX + "tests.txt"), sep="|", dtype=str, encoding="latin1",
            usecols=test_columns, chunksize=100_000, low_memory=False,
        ):
            chunk = strip_frame(chunk)
            chunk["grupo_organismos"] = chunk["species_number"].map(species_to_group)
            expected_hours = chunk["grupo_organismos"].map({key: value["hours"] for key, value in GROUPS.items()})
            selected = chunk.loc[
                chunk["grupo_organismos"].notna()
                & chunk["media_type"].fillna("").str.startswith("FW")
                & chunk["exposure_duration_mean"].eq(expected_hours)
                & chunk["exposure_duration_unit"].eq("h")
            ]
            if not selected.empty:
                test_parts.append(selected)
        tests = pd.concat(test_parts, ignore_index=True)
        test_ids = set(tests["test_id"])

        result_columns = [
            "result_id", "test_id", "obs_duration_mean", "obs_duration_unit", "endpoint",
            "effect", "measurement", "conc1_mean_op", "conc1_mean", "conc1_unit", "endpoint_assigned",
        ]
        result_parts: list[pd.DataFrame] = []
        for chunk in pd.read_csv(
            archive.open(PREFIX + "results.txt"), sep="|", dtype=str, encoding="latin1",
            usecols=result_columns, chunksize=100_000, low_memory=False,
        ):
            chunk = strip_frame(chunk)
            selected = chunk.loc[chunk["test_id"].isin(test_ids)]
            if not selected.empty:
                result_parts.append(selected)
        results = pd.concat(result_parts, ignore_index=True)

        chemicals = strip_frame(pd.read_csv(
            archive.open(PREFIX + "validation/chemicals.txt"), sep="|", dtype=str, encoding="latin1"
        )).rename(columns={"cas_number": "test_cas"})

    raw = results.merge(tests, on="test_id", how="inner", validate="many_to_one")
    raw = raw.merge(
        chemicals[["test_cas", "chemical_name", "dtxsid"]], on="test_cas", how="left", validate="many_to_one"
    )
    valid_by_group = pd.Series(False, index=raw.index)
    for group, config in GROUPS.items():
        measurement = raw["measurement"].fillna("")
        valid_by_group |= (
            raw["grupo_organismos"].eq(group)
            & raw["endpoint"].eq(config["endpoint"])
            & raw["effect"].eq(config["effect"])
            & measurement.str.startswith(config["measurement_prefixes"])
        )
    raw["concentration_value"] = pd.to_numeric(raw["conc1_mean"], errors="coerce")
    raw["casrn"] = raw["test_cas"].map(casrn)
    raw["ec50_mg_l"] = [to_mg_l(value, unit) for value, unit in zip(raw["concentration_value"], raw["conc1_unit"])]
    valid_operator = raw["conc1_mean_op"].fillna("").isin(["", "="])
    qc = raw.loc[valid_by_group & valid_operator & raw["casrn"].ne("") & raw["ec50_mg_l"].gt(0)].copy()
    qc["log_ec50_mg_l"] = np.log10(qc["ec50_mg_l"])
    qc["endpoint_descripcion"] = qc["grupo_organismos"].map({key: value["description"] for key, value in GROUPS.items()})
    qc.to_csv(PROCESSED / "ecotox_multigrupo_qc.csv", index=False, encoding="utf-8-sig")

    training = (
        qc.groupby(["grupo_organismos", "casrn"], as_index=False)
        .agg(
            chemical_key=("casrn", "first"),
            chemical_name=("chemical_name", "first"),
            dtxsid=("dtxsid", "first"),
            endpoint_descripcion=("endpoint_descripcion", "first"),
            log_ec50_mg_l=("log_ec50_mg_l", "median"),
            ec50_mg_l=("ec50_mg_l", "median"),
            n_assays=("result_id", "count"),
            n_references=("reference_number", "nunique"),
            log_ec50_sd=("log_ec50_mg_l", "std"),
        )
    )
    training["smiles"] = ""
    training.to_csv(PROCESSED / "multigrupo_training_before_structures.csv", index=False, encoding="utf-8-sig")

    per_group = training.groupby("grupo_organismos").agg(
        sustancias_unicas=("chemical_key", "nunique"),
        ensayos_validos=("n_assays", "sum"),
    ).to_dict(orient="index")
    report = {
        "source": str(ZIP_FILE),
        "groups": {key: value["description"] for key, value in GROUPS.items()},
        "qc_records": int(len(qc)),
        "unique_group_chemical_pairs": int(len(training)),
        "per_group": {key: {name: int(value) for name, value in row.items()} for key, row in per_group.items()},
        "raw_qc_output": str(PROCESSED / "ecotox_multigrupo_qc.csv"),
        "training_output": str(PROCESSED / "multigrupo_training_before_structures.csv"),
    }
    (RESULTS / "multigrupo_extraction_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
