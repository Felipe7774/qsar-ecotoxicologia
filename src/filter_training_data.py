"""Prepara y audita el conjunto QSAR desde el consolidado de la tesis.

Este script no infiere condiciones ausentes. Conserva tres salidas:
1. todos los ensayos ecotoxicológicos encontrados;
2. candidatos que requieren verificar la duración de exposición;
3. registros que cumplen de forma estricta el endpoint QSAR definido.
"""

from __future__ import annotations

import json
import re
import unicodedata
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
THESIS = PROJECT.parent
INPUT = THESIS / "outputs" / "consolidado_ces_20261007" / "Consolidado_CES.xlsx"
PROCESSED = PROJECT / "data" / "processed"
RESULTS = PROJECT / "results"


def normalized_text(value: object) -> str:
    """Normaliza texto para comparar sin perder el valor original en la salida."""
    if pd.isna(value):
        return ""
    text = unicodedata.normalize("NFKD", str(value)).encode("ascii", "ignore").decode("ascii")
    return re.sub(r"\s+", " ", text).strip().lower()


def is_scalar_number(value: object) -> bool:
    """Acepta un número único; rechaza límites (>, <) e intervalos (a-b)."""
    if isinstance(value, (int, float, np.number)) and not pd.isna(value):
        return float(value) > 0
    text = str(value).strip().replace(",", ".")
    pattern = r"^[+]?\d+(?:\.\d+)?(?:e[+-]?\d+)?$"
    return bool(re.fullmatch(pattern, text, flags=re.IGNORECASE)) and float(text) > 0


def numeric_value(value: object) -> float:
    if not is_scalar_number(value):
        return np.nan
    return float(str(value).strip().replace(",", "."))


def concentration_to_mg_l(value: float, unit: object) -> float:
    """Convierte solamente unidades de concentración inequívocas a mg/L."""
    text = normalized_text(unit).replace(" ", "").replace("μ", "u").replace("µ", "u")
    text = text.replace("l-1", "/l")
    if text in {"mg/l", "mgperliter", "mgl"}:
        return value
    if text in {"ug/l", "mcg/l", "microg/l", "microgram/l"}:
        return value / 1_000
    if text in {"ng/l", "nanog/l", "nanogram/l"}:
        return value / 1_000_000
    if text in {"g/l", "gram/l"}:
        return value * 1_000
    return np.nan


def exact_48_hours(value: object) -> bool:
    """48 h explícitas; no convierte intervalos como 24-48 h en 48 h."""
    if isinstance(value, (int, float, np.number)) and not pd.isna(value):
        return np.isclose(float(value), 48.0)
    text = normalized_text(value)
    return bool(re.fullmatch(r"48(?:[.,]0+)?\s*(h|hr|hora|horas)?", text))


def requires_columns(frame: pd.DataFrame) -> dict[str, str]:
    aliases = {
        "record_type": "Tipo_registro",
        "chemical": "Contaminante_Normalizado",
        "casrn": "CAS",
        "species": "Especie en el ensayo",
        "endpoint": "Endpoint",
        "effect": "Efecto_medido",
        "duration": "Duración del ensayo",
        "medium": "Agua dulce o marina",
        "value": "EC50",
        "unit": "Unidad_EC50",
        "source": "Referencia",
        "row_id": "ID_consolidado",
    }
    missing = [name for name in aliases.values() if name not in frame.columns]
    if missing:
        raise ValueError(f"Faltan columnas requeridas: {missing}")
    return aliases


def export(frame: pd.DataFrame, path: Path) -> None:
    columns = [
        "ID_consolidado", "Contaminante_Normalizado", "CAS", "Especie en el ensayo",
        "Endpoint", "Efecto_medido", "Duración del ensayo", "Agua dulce o marina",
        "EC50", "Unidad_EC50", "ec50_mg_l", "log_ec50_mg_l", "Referencia",
    ]
    available = [column for column in columns if column in frame.columns]
    frame.loc[:, available].to_csv(path, index=False, encoding="utf-8-sig")


def main() -> None:
    if not INPUT.exists():
        raise FileNotFoundError(f"No se encontró el consolidado: {INPUT}")

    PROCESSED.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    data = pd.read_excel(INPUT, sheet_name="Consolidado", header=2)
    c = requires_columns(data)

    assays = data.loc[data[c["record_type"]].map(normalized_text).eq("ensayo ecotoxicologico")].copy()
    assays["ec50_mg_l"] = [
        concentration_to_mg_l(numeric_value(value), unit)
        for value, unit in zip(assays[c["value"]], assays[c["unit"]])
    ]
    assays["log_ec50_mg_l"] = np.where(
        assays["ec50_mg_l"].gt(0), np.log10(assays["ec50_mg_l"]), np.nan
    )

    species = assays[c["species"]].map(normalized_text)
    endpoint = assays[c["endpoint"]].map(normalized_text)
    effect = assays[c["effect"]].map(normalized_text)
    medium = assays[c["medium"]].map(normalized_text)
    chemical = assays[c["chemical"]].map(normalized_text)
    casrn = assays[c["casrn"]].map(normalized_text)

    exact_species = species.str.fullmatch(r"daphnia magna( \([^)]*\))?")
    exact_endpoint = endpoint.eq("ec50")
    immobilization = effect.str.contains("inmoviliz", na=False)
    freshwater = medium.str.contains("agua dulce|freshwater", na=False)
    identified = chemical.ne("") | casrn.ne("")
    numeric_concentration = assays["ec50_mg_l"].notna() & assays["ec50_mg_l"].gt(0)
    duration_48 = assays[c["duration"]].map(exact_48_hours)

    base_endpoint = assays.loc[exact_species & exact_endpoint & immobilization & freshwater & identified & numeric_concentration].copy()
    pending_duration = base_endpoint.loc[~duration_48.loc[base_endpoint.index]].copy()
    strict = base_endpoint.loc[duration_48.loc[base_endpoint.index]].copy()

    export(assays, PROCESSED / "ecotox_assays_from_consolidated.csv")
    export(pending_duration, PROCESSED / "qsar_candidates_pending_duration.csv")
    export(strict, PROCESSED / "qsar_candidates_strict.csv")

    if not strict.empty:
        strict["chemical_key"] = np.where(casrn.loc[strict.index].ne(""), casrn.loc[strict.index], chemical.loc[strict.index])
        training = (
            strict.groupby("chemical_key", as_index=False)
            .agg(
                chemical_name=(c["chemical"], "first"),
                casrn=(c["casrn"], "first"),
                log_ec50_mg_l=("log_ec50_mg_l", "median"),
                ec50_mg_l=("ec50_mg_l", "median"),
                n_assays=(c["row_id"], "count"),
            )
        )
    else:
        training = pd.DataFrame(columns=["chemical_key", "chemical_name", "casrn", "log_ec50_mg_l", "ec50_mg_l", "n_assays"])
    training.to_csv(PROCESSED / "qsar_training_ready.csv", index=False, encoding="utf-8-sig")

    report = {
        "source": str(INPUT),
        "total_rows_in_consolidated": int(len(data)),
        "ecotoxicology_assays": int(len(assays)),
        "daphnia_magna_ec50_immobilization_freshwater_with_identifier_and_numeric_value": int(len(base_endpoint)),
        "pending_manual_duration_verification": int(len(pending_duration)),
        "strict_48h_records": int(len(strict)),
        "unique_chemicals_ready_for_training": int(len(training)),
        "minimum_unique_chemicals_recommended": 300,
        "ready_to_train": bool(len(training) >= 300),
        "reason_if_not_ready": (
            "No se entrenó: faltan registros que cumplan de manera comprobable el endpoint completo "
            "y estructuras químicas (SMILES) para generar descriptores QSAR."
            if len(training) < 300 else ""
        ),
    }
    (RESULTS / "filtering_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")

    print(json.dumps(report, ensure_ascii=False, indent=2))
    print(f"\nCandidatos con duración por verificar: {PROCESSED / 'qsar_candidates_pending_duration.csv'}")
    print(f"Candidatos estrictos: {PROCESSED / 'qsar_candidates_strict.csv'}")


if __name__ == "__main__":
    main()
