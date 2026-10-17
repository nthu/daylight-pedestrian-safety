"""Run the fixed-state and complementary Mississippi registry analyses."""

from __future__ import annotations

import numpy as np
import pandas as pd
from functions.allocation import contrast, estimate
from functions.data import ensure_results, mississippi_counts, state_counts, write_csv
from functions.inference import jackknife_interval
from functions.registry import information_support, mantel_haenszel_log_contrast
from functions.specifications import BANDS, OTHER_FATAL, PEDESTRIAN
from statsmodels.stats.multitest import multipletests


def _state_analysis() -> dict[str, pd.DataFrame]:
    counts = state_counts()
    years = sorted(counts.year.unique())
    rows: list[dict[str, object]] = []
    deletions: list[dict[str, object]] = []
    support_frames: list[pd.DataFrame] = []
    fitted: dict[tuple[str, str], tuple[float, np.ndarray]] = {}

    endpoints = (PEDESTRIAN, OTHER_FATAL, "all_fatal_crashes")
    for endpoint in endpoints:
        for band in BANDS:
            log_ratio, table = mantel_haenszel_log_contrast(
                counts,
                unit="state",
                band=band,
                endpoint=endpoint,
                years=years,
            )
            deleted = np.asarray(
                [
                    mantel_haenszel_log_contrast(
                        counts,
                        unit="state",
                        band=band,
                        endpoint=endpoint,
                        years=[year for year in years if year != omitted],
                    )[0]
                    for omitted in years
                ]
            )
            interval = jackknife_interval(log_ratio, deleted)
            rows.append({"endpoint": endpoint, "band": band, **interval, "states": 47})
            fitted[endpoint, band] = log_ratio, deleted
            deletions.extend(
                {
                    "endpoint": endpoint,
                    "band": band,
                    "omitted_policy_year": int(year),
                    "log_ratio": float(value),
                }
                for year, value in zip(years, deleted, strict=False)
            )
            support = information_support(table, band)
            support["endpoint"] = endpoint
            support["band"] = band
            support_frames.append(support)

    direct_rows = []
    for band in BANDS:
        pedestrian, pedestrian_deleted = fitted[PEDESTRIAN, band]
        other, other_deleted = fitted[OTHER_FATAL, band]
        direct_rows.append(
            {
                "band": band,
                **jackknife_interval(
                    pedestrian - other,
                    pedestrian_deleted - other_deleted,
                ),
            }
        )

    support = pd.concat(support_frames, ignore_index=True)
    summary = (
        support.groupby(["endpoint", "band", "state"], observed=True)
        .agg(
            appearances=("appearances", "sum"),
            tables=("informative", "size"),
            informative_tables=("informative", "sum"),
            all_zero_tables=("all_zero", "sum"),
        )
        .reset_index()
    )
    return {
        "state_effects": pd.DataFrame(rows),
        "state_direct": pd.DataFrame(direct_rows),
        "state_deletions": pd.DataFrame(deletions),
        "state_support": summary,
    }


def _mississippi_analysis() -> dict[str, pd.DataFrame]:
    counts = mississippi_counts()
    years = sorted(counts.year.unique())
    endpoints = (
        "all_crashes",
        "injury_or_fatal_crashes",
        "nonfatal_injury_crashes",
        "property_damage_only_crashes",
    )

    rows: list[dict[str, object]] = []
    for endpoint in endpoints:
        for band in BANDS:
            log_ratio, _ = mantel_haenszel_log_contrast(
                counts,
                unit="county_key",
                band=band,
                endpoint=endpoint,
                years=years,
            )
            deleted = [
                mantel_haenszel_log_contrast(
                    counts,
                    unit="county_key",
                    band=band,
                    endpoint=endpoint,
                    years=[year for year in years if year != omitted],
                )[0]
                for omitted in years
            ]
            rows.append({"endpoint": endpoint, "band": band, **jackknife_interval(log_ratio, deleted)})

    adjusted = multipletests([row["p_value"] for row in rows], method="holm")[1]
    for row, p_value in zip(rows, adjusted, strict=False):
        row["holm_p"] = float(p_value)

    national = (
        counts.groupby(
            ["year", "season", "sham", "period", "band", "width_minutes", "endpoint"],
            observed=True,
        )
        .events.sum()
        .reset_index()
    )
    national["geo"] = "MS_fixed_county"

    direct_rows = []
    for band in BANDS:
        injury = estimate(
            national,
            band=band,
            endpoint="nonfatal_injury_crashes",
            geo="MS_fixed_county",
        )
        property_only = estimate(
            national,
            band=band,
            endpoint="property_damage_only_crashes",
            geo="MS_fixed_county",
        )
        direct_rows.append(
            {
                "band": band,
                **contrast(injury, property_only, "nonfatal_injury_vs_property_damage_only"),
            }
        )

    variance_rows = []
    for band in BANDS:
        full = estimate(national, band=band, endpoint="all_crashes", geo="MS_fixed_county")
        deleted = np.asarray(
            [
                estimate(
                    national,
                    band=band,
                    endpoint="all_crashes",
                    geo="MS_fixed_county",
                    years=[year for year in years if year != omitted],
                )["log_ratio"]
                for omitted in years
            ]
        )
        jackknife = jackknife_interval(full["log_ratio"], deleted)
        variance_rows.extend(
            [
                {
                    "band": band,
                    "method": "allocation_CR1",
                    "ratio": full["ratio"],
                    "se": full["se"],
                    "low": full["low"],
                    "high": full["high"],
                },
                {
                    "band": band,
                    "method": "allocation_year_jackknife",
                    **jackknife,
                },
            ]
        )

    return {
        "mississippi_effects": pd.DataFrame(rows),
        "mississippi_direct": pd.DataFrame(direct_rows),
        "mississippi_variance": pd.DataFrame(variance_rows),
    }


def run() -> dict[str, pd.DataFrame]:
    paths = ensure_results()
    outputs = {**_state_analysis(), **_mississippi_analysis()}
    for name, frame in outputs.items():
        write_csv(frame, paths["estimates"] / f"{name}.csv")

    print("\nFixed-state pedestrian estimates")
    print(
        outputs["state_effects"]
        .query("endpoint == @PEDESTRIAN")[["band", "ratio", "low", "high"]]
        .to_string(index=False)
    )
    print("\nMississippi all-crash estimates")
    print(
        outputs["mississippi_effects"]
        .query("endpoint == 'all_crashes'")[["band", "ratio", "low", "high", "holm_p"]]
        .to_string(index=False)
    )
    return outputs


if __name__ == "__main__":
    run()
