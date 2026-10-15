"""Obtiene SMILES desde PubChem para las sustancias ECOTOX y valida su lectura RDKit."""

from __future__ import annotations

import json
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from threading import Lock
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen

import pandas as pd
from rdkit import Chem


PROJECT = Path(__file__).resolve().parents[1]
INPUT = PROJECT / "data" / "processed" / "qsar_training_before_structures.csv"
OUTPUT = PROJECT / "data" / "processed" / "qsar_training_with_structures.csv"
REPORT = PROJECT / "results" / "structure_resolution_training_report.json"
PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{query}/property/IsomericSMILES,CanonicalSMILES/JSON"
MAX_WORKERS = 3
REQUEST_INTERVAL_SECONDS = 0.25
_rate_lock = Lock()
_last_request = 0.0


def wait_for_request_slot() -> None:
    """Limita el conjunto de solicitudes a cuatro por segundo para PubChem."""
    global _last_request
    with _rate_lock:
        delay = REQUEST_INTERVAL_SECONDS - (time.monotonic() - _last_request)
        if delay > 0:
            time.sleep(delay)
        _last_request = time.monotonic()


def query_pubchem(casrn: str) -> dict[str, str]:
    url = PUBCHEM.format(query=quote(casrn, safe=""))
    request = Request(url, headers={"User-Agent": "QSAR-ecotoxicologia-tesis/1.0"})
    try:
        wait_for_request_slot()
        with urlopen(request, timeout=30) as response:
            payload = json.loads(response.read().decode("utf-8"))
    except HTTPError as error:
        return {"chemical_key": casrn, "structure_status": f"http_{error.code}", "smiles": ""}
    except URLError as error:
        return {"chemical_key": casrn, "structure_status": f"network_error: {error.reason}", "smiles": ""}
    except TimeoutError:
        return {"chemical_key": casrn, "structure_status": "timeout", "smiles": ""}

    properties = payload.get("PropertyTable", {}).get("Properties", [])
    if not properties:
        return {"chemical_key": casrn, "structure_status": "no_match", "smiles": ""}
    record = properties[0]
    smiles = str(record.get("SMILES") or record.get("IsomericSMILES") or record.get("ConnectivitySMILES") or record.get("CanonicalSMILES") or "")
    molecule = Chem.MolFromSmiles(smiles) if smiles else None
    return {
        "chemical_key": casrn,
        "pubchem_cid": str(record.get("CID", "")),
        "smiles": smiles,
        "structure_status": "rdkit_valid" if molecule is not None else "rdkit_invalid",
    }


def main() -> None:
    if not INPUT.exists():
        raise FileNotFoundError(f"No existe {INPUT}. Ejecuta build_qsar_training_table.py primero.")
    source = pd.read_csv(INPUT, dtype={"chemical_key": str, "casrn": str})
    resolved_rows: list[dict[str, str]] = []
    keys = source["chemical_key"].dropna().astype(str).tolist()
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(query_pubchem, key): key for key in keys}
        for position, future in enumerate(as_completed(futures), start=1):
            resolved_rows.append(future.result())
            if position % 50 == 0 or position == len(keys):
                print(f"Estructuras consultadas: {position}/{len(keys)}", flush=True)
            time.sleep(0.05)

    resolved = pd.DataFrame(resolved_rows)
    final = source.drop(columns=["smiles"], errors="ignore").merge(resolved, on="chemical_key", how="left", validate="one_to_one")
    final.to_csv(OUTPUT, index=False, encoding="utf-8-sig")
    report = {
        "input_substances": int(len(source)),
        "pubchem_rdkit_valid": int(final["structure_status"].eq("rdkit_valid").sum()),
        "not_usable_for_rdkit": int((~final["structure_status"].eq("rdkit_valid")).sum()),
        "output": str(OUTPUT),
        "note": "La identidad se obtuvo por CASRN; los casos sin estructura o con SMILES no interpretable se excluyen del entrenamiento.",
    }
    REPORT.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
