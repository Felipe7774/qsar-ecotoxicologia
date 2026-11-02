"""Verifica que la versión pública contenga sus artefactos esenciales."""

from __future__ import annotations

import json
from pathlib import Path

import joblib


PROJECT = Path(__file__).resolve().parents[1]
MODEL_PATH = PROJECT / "models" / "qsar_multiorganismo.joblib"
AUDIT_PATH = PROJECT / "results" / "ecotox_audit_report.json"
OUTPUT_PATH = PROJECT / "results" / "public_release_check.json"


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Falta el artefacto requerido: {path.relative_to(PROJECT)}")


def main() -> None:
    required_files = [
        PROJECT / "app.py",
        PROJECT / "requirements.txt",
        PROJECT / "Dockerfile",
        MODEL_PATH,
        PROJECT / "docs" / "informe_tecnico_qsar.md",
        PROJECT / "docs" / "guia_usuario_local.md",
        PROJECT / "docs" / "protocolo_auditoria_ecotox.md",
        AUDIT_PATH,
    ]
    for path in required_files:
        require(path)

    bundle = joblib.load(MODEL_PATH)
    groups = bundle.get("groups", {})
    if not {"Crustáceos", "Peces", "Algas"}.issubset(groups):
        raise ValueError("Faltan grupos ecotoxicológicos en el modelo público.")
    for name, group_bundle in groups.items():
        if group_bundle.get("feature_mode") != "hybrid":
            raise ValueError(f"{name}: se esperaba un modelo híbrido.")
        if len(group_bundle.get("descriptor_features", [])) != 10:
            raise ValueError(f"{name}: faltan descriptores estructurales.")
        if len(group_bundle.get("reference_fingerprints", [])) == 0:
            raise ValueError(f"{name}: faltan fingerprints de referencia.")

    audit = json.loads(AUDIT_PATH.read_text(encoding="utf-8"))
    if audit.get("sample_size") != 30 or audit.get("accepted_automatically") != 30:
        raise ValueError("La auditoría ECOTOX no cumple el resultado documentado de 30/30.")

    report = {
        "status": "ready_for_public_hosting",
        "model_type": bundle["model_type"],
        "aggregation_groups": bundle["aggregation_groups"],
        "groups": {name: {"descriptor_count": len(group["descriptor_features"]), "applicability_similarity_threshold": group["applicability_similarity_threshold"]} for name, group in groups.items()},
        "ecotox_audit": audit,
    }
    OUTPUT_PATH.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"Comprobación superada: {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
