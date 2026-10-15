"""Crea una muestra estratificada y reproducible para auditar datos ECOTOX."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
TRAINING = PROJECT / "data" / "processed" / "qsar_training_ready.csv"
QC_RECORDS = PROJECT / "data" / "processed" / "ecotox_daphnia_ec50_48h_immobilization_qc.csv"
OUTPUT = PROJECT / "data" / "processed" / "ecotox_audit_sample_30.csv"
SEED = 20261015
SAMPLE_SIZE = 30


def main() -> None:
    if not TRAINING.exists() or not QC_RECORDS.exists():
        raise FileNotFoundError("Faltan el conjunto de entrenamiento o los registros ECOTOX depurados.")
    training = pd.read_csv(TRAINING, dtype={"chemical_key": str, "casrn": str})
    qc = pd.read_csv(QC_RECORDS, dtype={"casrn": str, "result_id": str, "reference_number": str})
    if len(training) < SAMPLE_SIZE:
        raise ValueError(f"Se necesitan al menos {SAMPLE_SIZE} sustancias para la auditoría.")

    training = training.copy()
    training["audit_stratum"] = pd.qcut(
        training["log_ec50_mg_l"], q=3, labels=["toxicidad_alta", "toxicidad_media", "toxicidad_baja"], duplicates="drop"
    )
    rng = np.random.default_rng(SEED)
    per_stratum = SAMPLE_SIZE // training["audit_stratum"].nunique()
    sample_parts = [
        group.sample(n=min(per_stratum, len(group)), random_state=int(rng.integers(0, 2**31 - 1)))
        for _, group in training.groupby("audit_stratum", observed=True)
    ]
    sample = pd.concat(sample_parts, ignore_index=True)
    if len(sample) < SAMPLE_SIZE:
        remaining = training.drop(index=sample.index, errors="ignore")
        sample = pd.concat([sample, remaining.sample(n=SAMPLE_SIZE - len(sample), random_state=SEED)], ignore_index=True)

    evidence_columns = [
        "casrn", "result_id", "reference_number", "test_id", "chemical_name", "ec50_value",
        "conc1_unit", "measurement", "effect", "endpoint", "exposure_duration_mean",
        "exposure_duration_unit", "media_type", "species_number",
    ]
    evidence = qc.loc[:, [column for column in evidence_columns if column in qc.columns]].copy()
    evidence = evidence.sort_values("result_id").drop_duplicates("casrn")
    audit = sample.merge(evidence, on="casrn", how="left", suffixes=("_training", "_ecotox"))
    audit = audit.rename(columns={"chemical_name_training": "chemical_name", "chemical_name_ecotox": "chemical_name_ecotox"})
    audit["revisor"] = ""
    audit["fecha_revision"] = ""
    audit["casrn_verificado"] = ""
    audit["especie_verificada"] = ""
    audit["duracion_48h_verificada"] = ""
    audit["endpoint_ec50_verificado"] = ""
    audit["inmovilizacion_verificada"] = ""
    audit["unidad_y_valor_verificados"] = ""
    audit["decision"] = "Pendiente"
    audit["notas"] = ""
    ordered = [
        "chemical_key", "casrn", "chemical_name", "audit_stratum", "log_ec50_mg_l", "ec50_mg_l",
        "n_assays", "n_references", "log_ec50_sd", "result_id", "reference_number", "test_id",
        "ec50_value", "conc1_unit", "measurement", "effect", "endpoint", "exposure_duration_mean",
        "exposure_duration_unit", "media_type", "species_number", "revisor", "fecha_revision",
        "casrn_verificado", "especie_verificada", "duracion_48h_verificada", "endpoint_ec50_verificado",
        "inmovilizacion_verificada", "unidad_y_valor_verificados", "decision", "notas",
    ]
    audit = audit.loc[:, [column for column in ordered if column in audit.columns]]
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    audit.to_csv(OUTPUT, index=False, encoding="utf-8-sig")
    print(f"Muestra creada: {len(audit)} sustancias")
    print(audit["audit_stratum"].value_counts().to_string())
    print(f"Salida: {OUTPUT}")


if __name__ == "__main__":
    main()
