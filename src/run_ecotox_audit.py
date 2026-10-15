"""Verifica la muestra de auditoría contra las tablas originales ECOTOX ASCII."""

from __future__ import annotations

import json
import re
import zipfile
from datetime import date
from pathlib import Path

import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
ZIP_FILE = PROJECT / "data" / "raw" / "ecotox_ascii_09_15_2026.zip"
SAMPLE = PROJECT / "data" / "processed" / "ecotox_audit_sample_30.csv"
OUTPUT = PROJECT / "data" / "processed" / "ecotox_audit_completed.csv"
REPORT = PROJECT / "results" / "ecotox_audit_report.json"
PREFIX = "ecotox_ascii_09_15_2026/"
VALID_UNITS = {"mg/l", "ai mg/l", "ug/l", "ai ug/l", "ae ug/l", "ng/l", "ai ng/l", "ppm", "ppb", "ppt"}


def casrn(value: object) -> str:
    digits = re.sub(r"\D", "", str(value))
    return f"{digits[:-3]}-{digits[-3:-1]}-{digits[-1]}" if len(digits) >= 5 else ""


def is_positive_number(value: object) -> bool:
    try:
        return float(value) > 0
    except (TypeError, ValueError):
        return False


def main() -> None:
    if not ZIP_FILE.exists() or not SAMPLE.exists():
        raise FileNotFoundError("Falta el archivo ASCII ECOTOX o la muestra de auditoría.")
    sample = pd.read_csv(SAMPLE, dtype={"result_id": str, "test_id": str, "casrn": str})
    result_ids = set(sample["result_id"].dropna().astype(str))
    with zipfile.ZipFile(ZIP_FILE) as archive:
        results = pd.read_csv(
            archive.open(PREFIX + "results.txt"), sep="|", dtype=str, encoding="latin1",
            usecols=["result_id", "test_id", "endpoint", "measurement", "conc1_mean_op", "conc1_mean", "conc1_unit"],
        )
        results = results.loc[results["result_id"].astype(str).isin(result_ids)].copy()
        test_ids = set(results["test_id"].dropna().astype(str))
        tests = pd.read_csv(
            archive.open(PREFIX + "tests.txt"), sep="|", dtype=str, encoding="latin1",
            usecols=["test_id", "test_cas", "species_number", "exposure_duration_mean", "exposure_duration_unit", "media_type"],
        )
        tests = tests.loc[tests["test_id"].astype(str).isin(test_ids)].copy()

    source = results.merge(tests, on="test_id", how="left", validate="one_to_one")
    source["casrn_source"] = source["test_cas"].map(casrn)
    source["valid_cas"] = source["casrn_source"].ne("")
    source["valid_species"] = source["species_number"].str.strip().eq("5")
    source["valid_duration"] = source["exposure_duration_mean"].str.strip().eq("48") & source["exposure_duration_unit"].str.strip().eq("h")
    source["valid_medium"] = source["media_type"].str.strip().eq("FW")
    source["valid_endpoint"] = source["endpoint"].str.strip().eq("EC50")
    source["valid_immobilization"] = source["measurement"].str.strip().str.startswith("IMBL", na=False)
    source["valid_operator"] = source["conc1_mean_op"].fillna("").str.strip().isin(["", "="])
    source["valid_value"] = source["conc1_mean"].map(is_positive_number)
    source["valid_unit"] = source["conc1_unit"].fillna("").str.strip().str.casefold().isin(VALID_UNITS)
    checks = [
        "valid_cas", "valid_species", "valid_duration", "valid_medium", "valid_endpoint",
        "valid_immobilization", "valid_operator", "valid_value", "valid_unit",
    ]
    source["automatic_decision"] = source[checks].all(axis=1).map({True: "Aceptar", False: "Revisar"})
    merged = sample.merge(
        source[["result_id", "automatic_decision", *checks]], on="result_id", how="left", validate="one_to_one"
    )
    merged["revisor"] = "Auditoría automática ECOTOX ASCII"
    merged["fecha_revision"] = date.today().isoformat()
    merged["casrn_verificado"] = merged["valid_cas"].map({True: "Sí", False: "No"})
    merged["especie_verificada"] = merged["valid_species"].map({True: "Sí", False: "No"})
    merged["duracion_48h_verificada"] = merged["valid_duration"].map({True: "Sí", False: "No"})
    merged["endpoint_ec50_verificado"] = merged["valid_endpoint"].map({True: "Sí", False: "No"})
    merged["inmovilizacion_verificada"] = merged["valid_immobilization"].map({True: "Sí", False: "No"})
    merged["unidad_y_valor_verificados"] = (merged["valid_unit"] & merged["valid_value"] & merged["valid_operator"]).map({True: "Sí", False: "No"})
    merged["decision"] = merged["automatic_decision"] + " (automático)"
    merged["notas"] = "Trazabilidad verificada contra la descarga ASCII ECOTOX; no sustituye la consulta manual de la publicación primaria."
    merged.to_csv(OUTPUT, index=False, encoding="utf-8-sig")

    report = {
        "sample_size": int(len(merged)),
        "accepted_automatically": int(merged["automatic_decision"].eq("Aceptar").sum()),
        "requires_review": int(merged["automatic_decision"].eq("Revisar").sum()),
        "source": str(ZIP_FILE),
        "output": str(OUTPUT),
        "scope": "Verificación automática de criterios de inclusión y trazabilidad contra tablas ECOTOX ASCII.",
    }
    REPORT.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
