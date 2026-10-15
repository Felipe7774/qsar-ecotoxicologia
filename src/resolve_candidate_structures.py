"""Resuelve estructuras químicas de candidatos mediante PubChem PUG-REST.

Las coincidencias quedan marcadas como ``pendiente_revision``: una búsqueda por
nombre o CAS ayuda a obtener el SMILES, pero la inclusión QSAR requiere que la
persona responsable confirme la identidad contra la fuente ecotoxicológica.
"""

from __future__ import annotations

import csv
import json
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen

import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
INPUT = PROJECT / "data" / "processed" / "qsar_candidates_pending_duration.csv"
OUTPUT = PROJECT / "data" / "processed" / "qsar_candidates_with_structures.csv"
REVIEW = PROJECT / "data" / "processed" / "qsar_manual_review.csv"
REPORT = PROJECT / "results" / "structure_resolution_report.json"
PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{query}/property/CanonicalSMILES,IsomericSMILES,MolecularFormula,IUPACName/JSON"


def first_text(*values: object) -> str:
    for value in values:
        if pd.notna(value) and str(value).strip():
            return str(value).strip()
    return ""


def get_pubchem_properties(query: str) -> tuple[dict[str, object], str]:
    """Devuelve una propiedad PubChem o un estado legible para auditoría."""
    url = PUBCHEM.format(query=quote(query, safe=""))
    request = Request(url, headers={"User-Agent": "QSAR-ecotoxicologia-tesis/1.0"})
    try:
        with urlopen(request, timeout=30) as response:
            payload = json.loads(response.read().decode("utf-8"))
    except HTTPError as error:
        return {}, f"http_{error.code}"
    except URLError as error:
        return {}, f"network_error: {error.reason}"
    except TimeoutError:
        return {}, "timeout"
    except json.JSONDecodeError:
        return {}, "invalid_json"

    properties = payload.get("PropertyTable", {}).get("Properties", [])
    if not properties:
        return {}, "no_match"
    return properties[0], "resolved"


def main() -> None:
    if not INPUT.exists():
        raise FileNotFoundError(
            f"No existe {INPUT}. Ejecuta primero src/filter_training_data.py."
        )

    candidates = pd.read_csv(INPUT)
    rows: list[dict[str, object]] = []
    for _, candidate in candidates.iterrows():
        casrn = first_text(candidate.get("CAS"))
        chemical_name = first_text(candidate.get("Contaminante_Normalizado"))
        queries = [item for item in [casrn, chemical_name] if item]
        properties: dict[str, object] = {}
        status = "no_identifier"
        used_query = ""
        for query in queries:
            properties, status = get_pubchem_properties(query)
            used_query = query
            if status == "resolved":
                break
            time.sleep(0.25)

        row = candidate.to_dict()
        row.update(
            {
                "pubchem_query": used_query,
                "pubchem_status": status,
                "pubchem_cid": properties.get("CID", ""),
                "smiles_connectivity": first_text(
                    properties.get("ConnectivitySMILES"), properties.get("CanonicalSMILES")
                ),
                "smiles_isomeric": first_text(properties.get("SMILES"), properties.get("IsomericSMILES")),
                "molecular_formula": first_text(properties.get("MolecularFormula")),
                "iupac_name": first_text(properties.get("IUPACName")),
                "identity_status": "pendiente_revision" if status == "resolved" else "sin_resolver",
            }
        )
        rows.append(row)
        time.sleep(0.25)

    resolved = pd.DataFrame(rows)
    resolved.to_csv(OUTPUT, index=False, encoding="utf-8-sig", quoting=csv.QUOTE_MINIMAL)

    review_columns = [
        "ID_consolidado", "Contaminante_Normalizado", "CAS", "EC50", "Unidad_EC50",
        "Referencia", "pubchem_query", "pubchem_status", "pubchem_cid",
        "smiles_isomeric", "molecular_formula", "iupac_name", "identity_status",
    ]
    review = resolved.loc[:, [column for column in review_columns if column in resolved.columns]].copy()
    review["duracion_verificada_h"] = ""
    review["fuente_duracion"] = ""
    review["smiles_confirmado"] = ""
    review["incluir_en_qsar"] = ""
    review["notas_revision"] = ""
    review.to_csv(REVIEW, index=False, encoding="utf-8-sig", quoting=csv.QUOTE_MINIMAL)

    report = {
        "input_candidates": int(len(candidates)),
        "pubchem_resolved": int(resolved["pubchem_status"].eq("resolved").sum()),
        "pending_identity_review": int(resolved["identity_status"].eq("pendiente_revision").sum()),
        "unresolved": int(resolved["identity_status"].eq("sin_resolver").sum()),
        "review_file": str(REVIEW),
        "note": "Las estructuras resueltas automáticamente requieren confirmación antes de entrenar.",
    }
    REPORT.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
