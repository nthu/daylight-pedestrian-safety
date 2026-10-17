"""Create manuscript tables and all ten figures from computed results."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats
from functions.data import ensure_results, write_csv
from functions.legacy import estimate_legacy_weeks
from functions.plots import run as render_figures


def _main_tables(paths: dict[str, object]) -> None:
    estimates = pd.read_csv(paths["estimates"] / "national_models.csv").set_index("specification")
    direct = pd.read_csv(paths["estimates"] / "direct_contrasts.csv").set_index("contrast")

    rows = []
    for endpoint, prefix in [
        ("Fatal pedestrian crashes", "primary"),
        ("Other fatal crashes", "other_fatal_crashes"),
        ("All fatal crashes", "all_fatal_crashes"),
    ]:
        row = {"endpoint_or_contrast": endpoint}
        for band in ("evening", "morning"):
            fitted = estimates.loc[f"{prefix}_{band}"]
            row[f"{band}_ratio"] = fitted.ratio
            row[f"{band}_low"] = fitted.low
            row[f"{band}_high"] = fitted.high
        rows.append(row)

    direct_row = {"endpoint_or_contrast": "Pedestrian / other: ratio of RCRs"}
    for band in ("evening", "morning"):
        fitted = direct.loc[f"pedestrian_vs_other_{band}"]
        direct_row[f"{band}_ratio"] = fitted.ratio
        direct_row[f"{band}_low"] = fitted.low
        direct_row[f"{band}_high"] = fitted.high
    rows.append(direct_row)
    write_csv(pd.DataFrame(rows), paths["tables"] / "table_2_primary_contrasts.csv")

    leave_year = pd.read_csv(paths["estimates"] / "leave_one_year.csv")
    conditional = pd.read_csv(paths["estimates"] / "conditional_blocks.csv")
    wild = pd.read_csv(paths["estimates"] / "wild_cluster_calibration.csv").set_index("comparison")
    five = pd.read_csv(paths["estimates"] / "five_center_curvature.csv")
    calendar = pd.read_csv(paths["estimates"] / "calendar_sensitivity.csv")

    def jackknife_row(band: str) -> tuple[float, float, float]:
        values = np.log(leave_year.loc[leave_year.band.eq(band), "ratio"].to_numpy())
        blocks = len(values)
        standard_error = np.sqrt((blocks - 1) / blocks * np.sum((values - values.mean()) ** 2))
        beta = float(np.log(estimates.loc[f"primary_{band}", "ratio"]))
        critical = stats.t.ppf(0.975, blocks - 1)
        return (
            float(np.exp(beta)),
            float(np.exp(beta - critical * standard_error)),
            float(np.exp(beta + critical * standard_error)),
        )

    diagnostic_rows = []
    for label, source, scope in [
        ("Main three-center allocation, CR1", "main", "2001–2024"),
        ("Exact policy-year multiplier calibration", "wild", "2001–2024"),
        ("Main allocation, policy-year jackknife", "jackknife", "2001–2024"),
        ("Fixed-margin conditional, policy-year jackknife", "conditional", "2001–2024"),
        ("Fixed comparison-center solar threshold", "fixed", "2001–2024"),
        ("Nonoverlapping ±42-day calendar", "nonoverlap", "2002–2023"),
        ("Five-center quadratic background", "five", "2002–2023"),
    ]:
        row = {"analysis": label, "scope": scope}
        for band in ("evening", "morning"):
            if source == "main":
                fitted = estimates.loc[f"primary_{band}"]
                values = (fitted.ratio, fitted.low, fitted.high)
            elif source == "wild":
                fitted = wild.loc[f"pedestrian_{band}"]
                values = (fitted.ratio, fitted.wild_low, fitted.wild_high)
            elif source == "jackknife":
                values = jackknife_row(band)
            elif source == "conditional":
                fitted = conditional.query("band == @band and block_size_years == 1").iloc[0]
                values = (fitted.ratio, fitted.low, fitted.high)
            elif source == "fixed":
                fitted = estimates.loc[f"fixed_clock_{band}"]
                values = (fitted.ratio, fitted.low, fitted.high)
            elif source == "nonoverlap":
                fitted = calendar.query("specification == 'nonoverlap_14d_gap42' and band == @band").iloc[0]
                values = (fitted.ratio, fitted.low, fitted.high)
            else:
                fitted = five.query(
                    "model == 'five_center_quadratic' and endpoint == 'fatal_pedestrian_crashes' and band == @band"
                ).iloc[0]
                values = (fitted.ratio, fitted.low, fitted.high)
            row[f"{band}_ratio"], row[f"{band}_low"], row[f"{band}_high"] = map(float, values)
        diagnostic_rows.append(row)
    write_csv(pd.DataFrame(diagnostic_rows), paths["tables"] / "table_3_key_diagnostics.csv")


def run() -> dict[str, object]:
    paths = ensure_results()

    legacy, pretests = estimate_legacy_weeks()
    write_csv(legacy, paths["estimates"] / "legacy_weekly_estimates.csv")
    write_csv(pretests, paths["estimates"] / "legacy_weekly_pretests.csv")

    _main_tables(paths)
    figure_result = render_figures(paths["figures"])
    print(f"\nGenerated {figure_result['figures']} figures in {paths['figures']}")
    return {"legacy_rows": len(legacy), "figure_result": figure_result}


if __name__ == "__main__":
    run()
