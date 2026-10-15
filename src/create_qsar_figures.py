"""Crea figuras reproducibles de validación para el manuscrito de tesis."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
RESULTS = PROJECT / "results"


def main() -> None:
    oof_path = RESULTS / "out_of_fold_predictions.csv"
    importance_path = RESULTS / "feature_importance.csv"
    if not oof_path.exists() or not importance_path.exists():
        raise FileNotFoundError("Ejecuta src/train_qsar.py antes de generar las figuras.")
    oof = pd.read_csv(oof_path)
    importance = pd.read_csv(importance_path).sort_values("importance")

    plt.style.use("seaborn-v0_8-whitegrid")
    fig, axis = plt.subplots(figsize=(6.6, 5.6), dpi=180)
    axis.scatter(
        oof["observed_log_ec50_mg_l"], oof["predicted_log_ec50_mg_l"],
        alpha=0.72, color="#1769aa", edgecolors="white", linewidths=0.35,
    )
    lower = min(oof["observed_log_ec50_mg_l"].min(), oof["predicted_log_ec50_mg_l"].min())
    upper = max(oof["observed_log_ec50_mg_l"].max(), oof["predicted_log_ec50_mg_l"].max())
    axis.plot([lower, upper], [lower, upper], "--", color="#4a4a4a", linewidth=1, label="Predicción perfecta")
    axis.set_xlabel("EC50 observada, log10(mg/L)")
    axis.set_ylabel("EC50 predicha fuera de muestra, log10(mg/L)")
    axis.legend(frameon=True)
    fig.tight_layout()
    fig.savefig(RESULTS / "observed_vs_predicted.png", bbox_inches="tight")
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(7.2, 5.6), dpi=180)
    axis.barh(importance["variable"], importance["importance"], color="#e57a44")
    axis.set_xlabel("Importancia en Extra Trees")
    axis.set_ylabel("Descriptor molecular")
    fig.tight_layout()
    fig.savefig(RESULTS / "feature_importance.png", bbox_inches="tight")
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(6.6, 4.6), dpi=180)
    axis.scatter(
        oof["predicted_log_ec50_mg_l"], oof["residual_log_ec50_mg_l"],
        alpha=0.72, color="#198754", edgecolors="white", linewidths=0.35,
    )
    axis.axhline(0, color="#4a4a4a", linestyle="--", linewidth=1)
    axis.set_xlabel("EC50 predicha fuera de muestra, log10(mg/L)")
    axis.set_ylabel("Residuo: observada − predicha")
    fig.tight_layout()
    fig.savefig(RESULTS / "residuals_vs_predicted.png", bbox_inches="tight")
    plt.close(fig)
    print(f"Figuras creadas en: {RESULTS}")


if __name__ == "__main__":
    main()
