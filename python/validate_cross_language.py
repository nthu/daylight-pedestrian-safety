"""Compare the independently recomputed R and Python outputs.

Run both complete workflows before this script.  The comparison uses explicit
keys and estimands; files may contain additional language-specific diagnostics.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

REPOSITORY = Path(__file__).resolve().parents[1]
R_RESULTS = REPOSITORY / "results" / "R" / "estimates"
PYTHON_RESULTS = REPOSITORY / "results" / "python" / "estimates"
REPORT = REPOSITORY / "results" / "cross_language_validation.json"


@dataclass(frozen=True)
class Comparison:
    file: str
    keys: tuple[str, ...]
    fields: tuple[str, ...]
    tolerance: float = 3e-6


COMPARISONS = (
    Comparison("national_models.csv", ("specification",), ("ratio", "low", "high", "se", "df")),
    Comparison("direct_contrasts.csv", ("contrast",), ("ratio", "low", "high", "se", "df")),
    Comparison("reference_sensitivity.csv", ("specification",), ("ratio", "low", "high", "se", "df")),
    Comparison("calendar_sensitivity.csv", ("specification", "band"), ("ratio", "low", "high", "se", "df")),
    Comparison("dynamic_estimates.csv", ("week", "band", "season"), ("ratio", "low", "high", "se", "df")),
    Comparison("conditional_blocks.csv", ("band", "block_size_years"), ("ratio", "low", "high", "se", "df")),
    Comparison("pre2022.csv", ("band", "endpoint"), ("ratio", "low", "high", "se", "df")),
    Comparison("annual_slope_conditional.csv", ("band", "endpoint"), ("ratio", "low", "high", "se", "df")),
    Comparison("common_time.csv", ("band", "endpoint"), ("ratio", "low", "high", "se", "df")),
    Comparison("whole_day.csv", ("endpoint", "season", "method"), ("ratio", "low", "high")),
    Comparison("clock_shift.csv", ("clock_shift_minutes", "band", "endpoint"), ("ratio", "low", "high", "se", "df")),
    Comparison(
        "five_center_curvature.csv",
        ("model", "endpoint", "band"),
        (
            "ratio",
            "low",
            "high",
            "se",
            "df",
            "quadratic_wald",
            "quadratic_p_value",
            "autumn_quadratic",
            "spring_quadratic",
        ),
    ),
    Comparison(
        "five_center_deletions.csv",
        ("endpoint", "band", "omitted_policy_year"),
        ("log_ratio", "ratio"),
    ),
    Comparison(
        "wild_cluster_calibration.csv",
        ("comparison",),
        (
            "ratio",
            "cr1_low",
            "cr1_high",
            "cr1_p_value",
            "wild_p_value",
            "wild_critical_95",
            "wild_low",
            "wild_high",
            "observed_t",
            "clusters",
            "unique_sign_patterns",
            "equivalent_full_sign_patterns",
        ),
    ),
    Comparison(
        "season_heterogeneity.csv",
        ("band", "result_type", "label"),
        ("ratio", "low", "high", "se", "p_value", "df", "shared_years"),
    ),
    Comparison(
        "leave_one_year.csv",
        ("band", "omitted_policy_year"),
        ("ratio", "low", "high", "full_ratio", "full_low", "full_high"),
    ),
    Comparison(
        "policy_era_interaction.csv",
        ("band", "result_type", "label"),
        ("ratio", "low", "high", "se", "p_value", "df", "years_in_era"),
    ),
    Comparison(
        "period_interactions.csv",
        ("model", "band", "result_type", "label"),
        ("ratio", "low", "high", "se", "p_value", "df", "years_in_era"),
    ),
    Comparison(
        "alaska_inclusion.csv",
        ("geography", "endpoint", "band"),
        ("ratio", "low", "high", "se", "df", "informative_strata", "raw_strata"),
    ),
    Comparison(
        "arizona_negative_control.csv",
        ("endpoint", "band"),
        (
            "focal_events",
            "reference_events",
            "candidate_tables",
            "informative_tables",
            "all_zero_tables",
            "three_center_adjusted_ratio",
            "three_center_adjusted_low",
            "three_center_adjusted_high",
            "reduced_actual_date_only_ratio",
            "reduced_actual_date_only_low",
            "reduced_actual_date_only_high",
        ),
        tolerance=1e-5,
    ),
    Comparison(
        "sham_decomposition.csv",
        ("band", "endpoint", "sham"),
        ("gain_fraction", "loss_fraction", "log_contrast", "contrast", "contribution"),
    ),
    Comparison("triplet_support.csv", ("band", "endpoint", "year", "season"), ("support_size", "informative")),
    Comparison(
        "standardization_support.csv",
        ("band", "endpoint", "year", "season", "sham", "period"),
        ("total", "fraction", "gain_side", "weight"),
    ),
    Comparison("state_effects.csv", ("endpoint", "band"), ("ratio", "low", "high", "se", "df")),
    Comparison("state_direct.csv", ("band",), ("ratio", "low", "high", "se", "df")),
    Comparison(
        "state_deletions.csv",
        ("endpoint", "band", "omitted_policy_year"),
        ("log_ratio",),
    ),
    Comparison(
        "state_support.csv",
        ("endpoint", "band", "state"),
        ("appearances", "tables", "informative_tables", "all_zero_tables"),
        tolerance=0,
    ),
    Comparison("mississippi_effects.csv", ("endpoint", "band"), ("ratio", "low", "high", "se", "df", "holm_p")),
    Comparison("mississippi_direct.csv", ("band", "contrast"), ("ratio", "low", "high", "se", "df")),
    Comparison("mississippi_variance.csv", ("band", "method"), ("ratio", "low", "high", "se")),
    Comparison(
        "traffic_support.csv",
        ("metric", "band", "minimum_sites"),
        (
            "jurisdictions",
            "sites",
            "ratio_equal_jurisdiction",
            "low",
            "high",
            "ratio_median_jurisdiction",
            "ratio_site_weighted",
        ),
    ),
    Comparison("legacy_weekly_estimates.csv", ("band", "season", "week"), ("ratio", "low", "high", "se", "df")),
    Comparison("legacy_weekly_pretests.csv", ("band", "season"), ("F", "p_value", "df_num", "df_denom")),
)


def _load(folder: Path, filename: str) -> pd.DataFrame:
    path = folder / filename
    if not path.is_file():
        raise FileNotFoundError(f"Run the complete workflow first; missing {path}")
    return pd.read_csv(path)


def _normalise_boolean(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.astype(int)
    if series.dtype == object:
        mapped = series.astype(str).str.lower().map({"true": 1, "false": 0})
        if mapped.notna().all():
            return mapped.astype(int)
    return series


def compare(specification: Comparison) -> dict[str, object]:
    r_frame = _load(R_RESULTS, specification.file)
    python_frame = _load(PYTHON_RESULTS, specification.file)
    required = set(specification.keys + specification.fields)
    for label, frame in (("R", r_frame), ("Python", python_frame)):
        missing = required - set(frame.columns)
        if missing:
            raise AssertionError(f"{specification.file}: {label} missing columns {sorted(missing)}")
        if frame.duplicated(list(specification.keys)).any():
            raise AssertionError(f"{specification.file}: duplicate {label} keys")

    r_selected = r_frame[list(specification.keys + specification.fields)].copy()
    python_selected = python_frame[list(specification.keys + specification.fields)].copy()
    merged = r_selected.merge(
        python_selected,
        on=list(specification.keys),
        how="outer",
        suffixes=("_R", "_python"),
        indicator=True,
        validate="one_to_one",
    )
    if not merged["_merge"].eq("both").all():
        unmatched = merged.loc[~merged["_merge"].eq("both"), list(specification.keys) + ["_merge"]]
        raise AssertionError(f"{specification.file}: unmatched rows\n{unmatched}")

    maximum_difference = 0.0
    field_differences: dict[str, float] = {}
    for field in specification.fields:
        left = _normalise_boolean(merged[f"{field}_R"])
        right = _normalise_boolean(merged[f"{field}_python"])
        left_numeric = pd.to_numeric(left, errors="raise").to_numpy(dtype=float)
        right_numeric = pd.to_numeric(right, errors="raise").to_numpy(dtype=float)
        if not np.all(np.isfinite(left_numeric) == np.isfinite(right_numeric)):
            raise AssertionError(f"{specification.file}: finite-value mismatch in {field}")
        difference = np.abs(left_numeric - right_numeric)
        finite = difference[np.isfinite(difference)]
        field_maximum = float(finite.max()) if finite.size else 0.0
        field_differences[field] = field_maximum
        maximum_difference = max(maximum_difference, field_maximum)
        if specification.tolerance == 0:
            passed = np.array_equal(left_numeric, right_numeric, equal_nan=True)
        else:
            passed = np.allclose(
                left_numeric,
                right_numeric,
                rtol=specification.tolerance,
                atol=1e-9,
                equal_nan=True,
            )
        if not passed:
            raise AssertionError(f"{specification.file}: {field} differs; max absolute difference {field_maximum:.6g}")

    return {
        "file": specification.file,
        "rows": int(len(merged)),
        "keys": list(specification.keys),
        "fields": list(specification.fields),
        "maximum_absolute_difference": maximum_difference,
        "field_differences": field_differences,
        "passed": True,
    }


def main() -> None:
    checks = [compare(specification) for specification in COMPARISONS]
    figure_counts = {}
    for language in ("R", "python"):
        folder = REPOSITORY / "results" / language / "figures"
        pdfs = sorted(folder.glob("Figure*.pdf"))
        pngs = sorted(folder.glob("Figure*.png"))
        csvs = sorted(folder.glob("Figure*.csv"))
        if not (len(pdfs) == len(pngs) == len(csvs) == 10):
            raise AssertionError(f"{language}: expected ten PDF, PNG and CSV figure outputs")
        figure_counts[language] = {"pdf": len(pdfs), "png": len(pngs), "csv": len(csvs)}

    report = {
        "passed": True,
        "independent_implementations": ["R", "Python"],
        "comparisons": checks,
        "figure_outputs": figure_counts,
        "files_compared": len(checks),
        "rows_compared": sum(item["rows"] for item in checks),
    }
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(report, indent=2) + "\n")
    print(
        f"Cross-language validation passed: {report['files_compared']} files, "
        f"{report['rows_compared']} keyed rows, ten figures per implementation."
    )


if __name__ == "__main__":
    main()
