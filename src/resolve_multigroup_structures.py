"""Resuelve y valida SMILES para el conjunto ecotoxicológico multi-organismo."""

from __future__ import annotations

import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

from resolve_training_structures import query_pubchem


PROJECT = Path(__file__).resolve().parents[1]
INPUT = PROJECT / "data" / "processed" / "multigrupo_training_before_structures.csv"
OUTPUT = PROJECT / "data" / "processed" / "multigrupo_training_with_structures.csv"
PREVIOUS = PROJECT / "data" / "processed" / "qsar_training_ready.csv"
REPORT = PROJECT / "results" / "multigrupo_structure_resolution_report.json"
MAX_WORKERS = 3


def cached_structures() -> dict[str, dict[str, str]]:
    if not PREVIOUS.exists():
        return {}
    previous = pd.read_csv(PREVIOUS, dtype={"chemical_key": str, "smiles": str})
    valid = previous.loc[previous.get("structure_status", pd.Series("", index=previous.index)).eq("rdkit_valid")]
    return {
        str(row.chemical_key): {
            "chemical_key": str(row.chemical_key),
            "pubchem_cid": str(getattr(row, "pubchem_cid", "")),
            "smiles": str(row.smiles),
            "structure_status": "rdkit_valid",
        }
        for row in valid.itertuples(index=False)
        if pd.notna(row.smiles) and str(row.smiles).strip()
    }


def main() -> None:
    if not INPUT.exists():
        raise FileNotFoundError(f"No existe {INPUT}. Ejecuta build_multigroup_ecotox_dataset.py primero.")
    source = pd.read_csv(INPUT, dtype={"chemical_key": str, "casrn": str})
    unique = source.drop_duplicates("chemical_key").copy()
    cache = cached_structures()
    missing = [key for key in unique["chemical_key"].dropna().astype(str) if key not in cache]
    resolved_rows = list(cache.values())
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(query_pubchem, key): key for key in missing}
        for position, future in enumerate(as_completed(futures), start=1):
            resolved_rows.append(future.result())
            if position % 50 == 0 or position == len(missing):
                print(f"Estructuras nuevas consultadas: {position}/{len(missing)}", flush=True)
    resolved = pd.DataFrame(resolved_rows).drop_duplicates("chemical_key", keep="last")
    final = source.drop(columns=["smiles"], errors="ignore").merge(resolved, on="chemical_key", how="left", validate="many_to_one")
    final.to_csv(OUTPUT, index=False, encoding="utf-8-sig")
    report = {
        "input_group_chemical_pairs": int(len(source)),
        "distinct_chemicals": int(unique["chemical_key"].nunique()),
        "reused_rdkit_valid_structures": int(len(cache)),
        "new_pubchem_queries": int(len(missing)),
        "rdkit_valid_group_chemical_pairs": int(final["structure_status"].eq("rdkit_valid").sum()),
        "output": str(OUTPUT),
    }
    REPORT.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
