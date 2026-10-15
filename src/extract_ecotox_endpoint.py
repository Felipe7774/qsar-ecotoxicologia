"""Extrae del archivo ASCII de ECOTOX el endpoint QSAR definido.

Lee el ZIP por bloques, sin descomprimir su contenido completo. Conserva primero
todos los EC50 comparables de *Daphnia magna* para auditar el código de efecto
antes de aplicar restricciones adicionales.
"""

from __future__ import annotations

import json
import zipfile
from pathlib import Path

import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
ZIP_FILE = PROJECT / "data" / "raw" / "ecotox_ascii_09_15_2026.zip"
PROCESSED = PROJECT / "data" / "processed"
RESULTS = PROJECT / "results"
PREFIX = "ecotox_ascii_09_15_2026/"


def strip_frame(frame: pd.DataFrame) -> pd.DataFrame:
    return frame.apply(lambda column: column.str.strip() if column.dtype == "object" else column)


def read_lookup(archive: zipfile.ZipFile, filename: str) -> pd.DataFrame:
    return strip_frame(pd.read_csv(archive.open(PREFIX + filename), sep="|", dtype=str, encoding="latin1"))


def main() -> None:
    if not ZIP_FILE.exists():
        raise FileNotFoundError(f"No existe la descarga ECOTOX: {ZIP_FILE}")
    PROCESSED.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)

    with zipfile.ZipFile(ZIP_FILE) as archive:
        species = read_lookup(archive, "validation/species.txt")
        daphnia_ids = set(species.loc[species["latin_name"].str.casefold().eq("daphnia magna"), "species_number"])
        if not daphnia_ids:
            raise ValueError("No se encontró Daphnia magna en el catálogo ECOTOX.")

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
            selected = chunk.loc[
                chunk["species_number"].isin(daphnia_ids)
                & chunk["exposure_duration_mean"].eq("48")
                & chunk["exposure_duration_unit"].eq("h")
                & chunk["media_type"].eq("FW")
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
            selected = chunk.loc[chunk["test_id"].isin(test_ids) & chunk["endpoint"].eq("EC50")]
            if not selected.empty:
                result_parts.append(selected)
        results = pd.concat(result_parts, ignore_index=True)

        chemicals = read_lookup(archive, "validation/chemicals.txt")
        chemicals = chemicals.rename(columns={"cas_number": "test_cas", "chemical_name": "chemical_name"})
        endpoint = results.merge(tests, on="test_id", how="left", validate="many_to_one")
        endpoint = endpoint.merge(chemicals[["test_cas", "chemical_name", "dtxsid"]], on="test_cas", how="left")
        endpoint["ec50_value"] = pd.to_numeric(endpoint["conc1_mean"], errors="coerce")

    raw_output = PROCESSED / "ecotox_daphnia_ec50_48h_freshwater_raw.csv"
    endpoint.to_csv(raw_output, index=False, encoding="utf-8-sig")
    effect_counts = endpoint["effect"].fillna("(vacío)").value_counts().to_dict()
    report = {
        "zip_source": str(ZIP_FILE),
        "daphnia_magna_species_numbers": sorted(daphnia_ids),
        "tests_daphnia_magna_48h_freshwater": int(len(tests)),
        "ec50_results_before_effect_and_unit_qc": int(len(endpoint)),
        "unique_cas_before_qc": int(endpoint["test_cas"].replace("", pd.NA).nunique()),
        "effect_code_counts": {str(key): int(value) for key, value in effect_counts.items()},
        "raw_output": str(raw_output),
    }
    (RESULTS / "ecotox_endpoint_extraction_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
