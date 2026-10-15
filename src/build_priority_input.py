"""Crea la entrada de priorización a partir de candidatos y el consolidado CES."""

from __future__ import annotations

import re
import unicodedata
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
THESIS = PROJECT.parent
CANDIDATES = PROJECT / "data" / "processed" / "qsar_manual_review.csv"
CONSOLIDATED = THESIS / "outputs" / "consolidado_ces_20261007" / "Consolidado_CES.xlsx"
OUTPUT = PROJECT / "data" / "processed" / "priority_input.csv"


def key(value: object) -> str:
    if pd.isna(value):
        return ""
    text = unicodedata.normalize("NFKD", str(value)).encode("ascii", "ignore").decode("ascii")
    return re.sub(r"[^a-z0-9]+", "", text.lower())


def numeric(value: object) -> float:
    if isinstance(value, (int, float, np.number)) and not pd.isna(value):
        return float(value)
    text = str(value).strip().replace(",", ".")
    return float(text) if re.fullmatch(r"[+]?[0-9]+(?:\.[0-9]+)?", text) else np.nan


def first_median(values: pd.Series) -> float:
    numeric_values = values.map(numeric).dropna()
    return float(numeric_values.median()) if not numeric_values.empty else np.nan


def main() -> None:
    if not CANDIDATES.exists() or not CONSOLIDATED.exists():
        raise FileNotFoundError("Falta el consolidado o qsar_manual_review.csv.")
    candidates = pd.read_csv(CANDIDATES, dtype=str)
    data = pd.read_excel(CONSOLIDATED, sheet_name="Consolidado", header=2)
    data["cas_key"] = data["CAS"].map(key)
    data["name_key"] = data["Contaminante_Normalizado"].map(key)
    rows: list[dict[str, object]] = []
    for _, candidate in candidates.iterrows():
        cas_key = key(candidate.get("CAS"))
        name_key = key(candidate.get("Contaminante_Normalizado"))
        matches = data.loc[(data["cas_key"].eq(cas_key) & data["cas_key"].ne("")) | data["name_key"].eq(name_key)]
        rows.append(
            {
                "compound_name": candidate.get("Contaminante_Normalizado", ""),
                "chemical_key": candidate.get("CAS") or candidate.get("ID_consolidado"),
                "smiles": candidate.get("smiles_isomeric", ""),
                "logkow": first_median(matches["logKow"]) if not matches.empty else np.nan,
                "occurrence_ug_l": first_median(matches["C_med_ugL"]) if not matches.empty else np.nan,
                "removal_pct": first_median(matches["Remocion_pct_estandarizada"]) if not matches.empty else np.nan,
                "matched_consolidated_rows": int(len(matches)),
            }
        )
    output = pd.DataFrame(rows)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(OUTPUT, index=False, encoding="utf-8-sig")
    print(f"Filas preparadas: {len(output)}")
    print(f"Con al menos un dato ambiental: {int(output[['logkow', 'occurrence_ug_l', 'removal_pct']].notna().any(axis=1).sum())}")
    print(f"Salida: {OUTPUT}")


if __name__ == "__main__":
    main()
