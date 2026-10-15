"""Filtra una exportación ECOTOX al endpoint definido para el proyecto."""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"


def clean_name(value: object) -> str:
    text = str(value).strip().lower()
    return re.sub(r"[^a-z0-9]+", "", text)


def find_column(frame: pd.DataFrame, aliases: list[str], required: bool = True) -> str | None:
    normalized = {clean_name(col): col for col in frame.columns}
    for alias in aliases:
        if clean_name(alias) in normalized:
            return normalized[clean_name(alias)]
    if required:
        raise ValueError(f"No se encontró una columna compatible con: {aliases}. Columnas: {list(frame.columns)}")
    return None


def read_export(path: Path) -> pd.DataFrame:
    if path.suffix.lower() in {".xlsx", ".xls"}:
        return pd.read_excel(path)
    return pd.read_csv(path, sep=None, engine="python", encoding_errors="replace")


def numeric_value(value: object) -> float:
    match = re.search(r"[-+]?\d+(?:[.,]\d+)?(?:[eE][-+]?\d+)?", str(value))
    if not match:
        return np.nan
    return float(match.group().replace(",", "."))


def hours_value(value: object) -> float:
    return numeric_value(value)


def to_mg_l(value: float, unit: object) -> float:
    unit_text = str(unit).lower().replace("μ", "µ").replace(" ", "")
    if "mg/l" in unit_text or "mg/liter" in unit_text:
        return value
    if "µg/l" in unit_text or "ug/l" in unit_text or "microgram/l" in unit_text:
        return value / 1000
    if "ng/l" in unit_text or "nanogram/l" in unit_text:
        return value / 1_000_000
    if "g/l" in unit_text or "gram/l" in unit_text:
        return value * 1000
    return np.nan


def main() -> None:
    candidates = sorted(RAW.glob("ecotox_export.*"))
    if not candidates:
        raise FileNotFoundError("Guarda la exportación como data/raw/ecotox_export.csv o .xlsx antes de ejecutar el script.")
    input_path = candidates[0]
    frame = read_export(input_path)
    print(f"Archivo: {input_path.name} | filas: {len(frame)}")
    print(f"Columnas detectadas: {list(frame.columns)}")

    species_col = find_column(frame, ["Species", "Scientific Name", "Species Scientific Name", "Organism"])
    endpoint_col = find_column(frame, ["Endpoint", "Endpoint Name"])
    effect_col = find_column(frame, ["Effect", "Effect Name"])
    duration_col = find_column(frame, ["Exposure Duration", "Exposure Time", "Duration", "Test Duration"])
    value_col = find_column(frame, ["Result", "Test Result", "Concentration", "Value", "Endpoint Value"])
    unit_col = find_column(frame, ["Units", "Unit", "Result Units", "Concentration Units"])
    chemical_col = find_column(frame, ["Chemical", "Chemical Name", "Common Name", "CAS"])
    smiles_col = find_column(frame, ["SMILES", "Canonical SMILES", "DSSTox Structure SMILES"], required=False)
    cas_col = find_column(frame, ["CAS", "CASRN", "CAS Number"], required=False)

    species = frame[species_col].astype(str).str.lower()
    endpoint = frame[endpoint_col].astype(str).str.lower()
    effect = frame[effect_col].astype(str).str.lower()
    duration = frame[duration_col].map(hours_value)
    selected = (
        species.str.contains("daphnia magna", na=False)
        & endpoint.str.contains(r"\bec50\b", regex=True, na=False)
        & effect.str.contains("immobil", na=False)
        & np.isclose(duration, 48, atol=0.01, equal_nan=False)
    )
    filtered = frame.loc[selected].copy()
    filtered["ec50_value_original"] = filtered[value_col].map(numeric_value)
    filtered["ec50_mg_l"] = [to_mg_l(v, u) for v, u in zip(filtered["ec50_value_original"], filtered[unit_col])]
    filtered = filtered[filtered["ec50_mg_l"].gt(0) & filtered["ec50_mg_l"].notna()].copy()
    filtered["log_ec50_mg_l"] = np.log10(filtered["ec50_mg_l"])

    output = pd.DataFrame(
        {
            "chemical_name": filtered[chemical_col].astype(str),
            "casrn": filtered[cas_col].astype(str) if cas_col else "",
            "smiles": filtered[smiles_col].astype(str) if smiles_col else "",
            "species": filtered[species_col].astype(str),
            "exposure_hours": duration.loc[filtered.index].astype(float),
            "effect": filtered[effect_col].astype(str),
            "endpoint": filtered[endpoint_col].astype(str),
            "ec50_value": filtered["ec50_value_original"],
            "ec50_unit": filtered[unit_col].astype(str),
            "ec50_mg_l": filtered["ec50_mg_l"],
            "log_ec50_mg_l": filtered["log_ec50_mg_l"],
        }
    )
    PROCESSED.mkdir(parents=True, exist_ok=True)
    output_path = PROCESSED / "ecotox_daphnia_48h_raw.csv"
    output.to_csv(output_path, index=False)
    summary = {"input": input_path.name, "input_rows": len(frame), "endpoint_rows": len(output), "unique_chemicals": int(output["chemical_name"].nunique())}
    (PROCESSED / "ecotox_filter_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"Ensayos válidos: {len(output)}")
    print(f"Compuestos únicos: {output['chemical_name'].nunique()}")
    print(f"Salida: {output_path}")


if __name__ == "__main__":
    main()
