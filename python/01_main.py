"""Fit the national primary models and direct road-user comparisons."""

from __future__ import annotations

import pandas as pd
from functions.allocation import contrast, estimate
from functions.data import ensure_results, national_counts, write_csv
from functions.specifications import national_specs


def run() -> dict[str, object]:
    counts = national_counts()
    paths = ensure_results()

    fitted: dict[str, dict[str, object]] = {}
    rows: list[dict[str, object]] = []

    for spec in national_specs():
        result = estimate(
            counts,
            band=spec.band,
            reference=spec.reference,
            endpoint=spec.endpoint,
            geo=spec.geography,
            width=spec.width_minutes,
            years=spec.years,
            seasons=spec.seasons,
            corrected=spec.corrected,
        )
        fitted[spec.name] = result
        rows.append(
            {
                "specification": spec.name,
                "band": spec.band,
                "reference": spec.reference,
                "endpoint": spec.endpoint,
                "geography": spec.geography,
                "width_minutes": spec.width_minutes,
                "corrected": spec.corrected,
                **{key: result[key] for key in ("ratio", "low", "high", "se", "df")},
            }
        )

    comparisons = [
        ("primary_evening", "primary_morning", "evening_vs_morning"),
        ("primary_evening", "other_fatal_crashes_evening", "pedestrian_vs_other_evening"),
        ("primary_morning", "other_fatal_crashes_morning", "pedestrian_vs_other_morning"),
        ("early_evening", "late_evening", "early_vs_late_evening"),
        ("east_evening", "west_evening", "east_vs_west_evening"),
        ("early_morning", "late_morning", "early_vs_late_morning"),
        ("east_morning", "west_morning", "east_vs_west_morning"),
    ]
    direct_rows = [contrast(fitted[left], fitted[right], name) for left, right, name in comparisons]

    estimates = pd.DataFrame(rows)
    direct = pd.DataFrame(direct_rows)
    write_csv(estimates, paths["estimates"] / "national_models.csv")
    write_csv(direct, paths["estimates"] / "direct_contrasts.csv")

    primary = estimates[estimates.specification.isin(["primary_evening", "primary_morning"])]
    print("\nNational primary estimates")
    print(primary[["specification", "ratio", "low", "high"]].to_string(index=False))
    print("\nDirect road-user comparisons")
    print(direct[["contrast", "ratio", "low", "high"]].to_string(index=False))

    return {"models": fitted, "estimates": estimates, "direct": direct}


if __name__ == "__main__":
    run()
