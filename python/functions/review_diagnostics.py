"""Targeted diagnostics requested during methodological review."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

from .allocation import contrast, estimate
from .data import national_counts, sensitivity_table
from .specifications import BANDS, OTHER_FATAL, PEDESTRIAN


def seasonal_heterogeneity(counts: pd.DataFrame | None = None) -> pd.DataFrame:
    """Return gain-oriented seasonal estimates and formal spring/autumn contrasts."""

    counts = national_counts() if counts is None else counts
    rows: list[dict[str, object]] = []
    for band in BANDS:
        fitted: dict[str, dict[str, object]] = {}
        for season in ("spring", "autumn"):
            result = estimate(counts, band=band, seasons=[season])
            fitted[season] = result
            rows.append(
                {
                    "band": band,
                    "result_type": "season_estimate",
                    "label": season,
                    "ratio": result["ratio"],
                    "low": result["low"],
                    "high": result["high"],
                    "se": result["se"],
                    "p_value": result["p_value"],
                    "df": result["df"],
                    "shared_years": len(result["years"]),
                }
            )
        comparison = contrast(fitted["spring"], fitted["autumn"], f"spring_vs_autumn_{band}")
        rows.append(
            {
                "band": band,
                "result_type": "season_contrast",
                "label": "spring_vs_autumn",
                **{key: comparison[key] for key in ("ratio", "low", "high", "se", "p_value", "df", "shared_years")},
            }
        )
    return pd.DataFrame(rows)


def leave_one_year_influence(counts: pd.DataFrame | None = None) -> pd.DataFrame:
    """Refit the main pedestrian model after omitting each policy year."""

    counts = national_counts() if counts is None else counts
    rows: list[dict[str, object]] = []
    years = list(range(2001, 2025))
    for band in BANDS:
        full = estimate(counts, band=band)
        for omitted in years:
            fitted = estimate(counts, band=band, years=[year for year in years if year != omitted])
            rows.append(
                {
                    "band": band,
                    "omitted_policy_year": omitted,
                    "ratio": fitted["ratio"],
                    "low": fitted["low"],
                    "high": fitted["high"],
                    "full_ratio": full["ratio"],
                    "full_low": full["low"],
                    "full_high": full["high"],
                }
            )
    return pd.DataFrame(rows)


def _era_model(
    counts: pd.DataFrame,
    *,
    band: str,
    groups: tuple[tuple[str, Callable[[int], bool]], ...],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    subset = counts[
        counts.geo.eq("DST_main")
        & counts.width_minutes.eq(40)
        & counts.endpoint.eq(PEDESTRIAN)
        & counts.band.isin([band, "night"])
    ].copy()
    table = subset.pivot(index=["year", "season", "sham"], columns=["period", "band"], values="events").sort_index()
    a = table[("post", band)]
    b = table[("post", "night")]
    c = table[("pre", band)]
    d = table[("pre", "night")]
    total_2x2 = a + b + c + d
    focal_total = a + c
    post_total = a + b
    informative = np.maximum(0, post_total - (total_2x2 - focal_total)) < np.minimum(post_total, focal_total)
    table = table[informative]
    meta = table.index.to_frame(index=False)
    strata = len(meta)

    columns: list[np.ndarray] = []
    names: list[str] = []
    for year, season in sorted(set(meta[["year", "season"]].itertuples(index=False, name=None))):
        columns.append((meta.year.eq(year) & meta.season.eq(season)).to_numpy(float))
        names.append(f"post_{year}_{season}")
    for season in sorted(meta.season.unique()):
        columns.append(meta.season.eq(season).to_numpy(float) * meta.sham.to_numpy(float) / 35.0)
        names.append(f"post_{season}_centre_slope")
    direction = np.where(meta.season.eq("spring"), 1.0, -1.0) * (1.0 if band == "evening" else -1.0)
    for label, selector in groups:
        active = np.asarray([float(selector(int(year))) for year in meta.year], dtype=float)
        columns.append(meta.sham.eq(0).to_numpy(float) * direction * active)
        names.append(f"effect_{label}")
    extras = np.column_stack(columns)

    design: list[np.ndarray] = []
    success: list[int] = []
    failure: list[int] = []
    policy_year: list[int] = []
    eye = np.eye(strata)
    for index, row in enumerate(meta.itertuples(index=False)):
        for post, period in enumerate(("pre", "post")):
            design.append(np.r_[eye[index], post * extras[index]])
            success.append(int(table.iloc[index][(period, band)]))
            failure.append(int(table.iloc[index][(period, "night")]))
            policy_year.append(int(row.year))

    X = np.asarray(design)
    successes = np.asarray(success)
    failures = np.asarray(failure)
    year = np.asarray(policy_year)
    fitted = sm.GLM(np.column_stack([successes, failures]), X, family=sm.families.Binomial()).fit(
        maxiter=500,
        tol=1e-11,
    )
    if not fitted.converged or not np.isfinite(fitted.params).all():
        raise RuntimeError("Era-interacted grouped-binomial model failed")

    totals = successes + failures
    probability = fitted.fittedvalues
    scores = X * (successes - totals * probability)[:, None]
    information = X.T @ (X * (totals * probability * (1.0 - probability))[:, None])
    eigenvalues = np.linalg.eigvalsh((information + information.T) / 2)
    if eigenvalues.min() <= eigenvalues.max() * 1e-12:
        raise ValueError("Era-interacted information matrix is not full rank")
    bread = np.linalg.inv(information)
    clusters = np.unique(year)
    cluster_scores = np.vstack([scores[year == cluster].sum(axis=0) for cluster in clusters])
    G = len(clusters)
    N = len(successes)
    K = X.shape[1]
    correction = G / (G - 1) * (N - 1) / (N - K)
    covariance = bread @ cluster_scores.T @ cluster_scores @ bread * correction
    full_names = [f"stratum_{index}" for index in range(strata)] + names
    df = G - 1
    critical = stats.t.ppf(0.975, df)

    effects: list[dict[str, object]] = []
    for label, selector in groups:
        index = full_names.index(f"effect_{label}")
        beta = float(fitted.params[index])
        se = float(np.sqrt(covariance[index, index]))
        effects.append(
            {
                "band": band,
                "result_type": "era_coefficient",
                "label": label,
                "log_ratio": beta,
                "ratio": float(np.exp(beta)),
                "se": se,
                "low": float(np.exp(beta - critical * se)),
                "high": float(np.exp(beta + critical * se)),
                "p_value": float(2 * stats.t.sf(abs(beta / se), df)),
                "df": df,
                "years_in_era": sum(selector(int(value)) for value in clusters),
                "from": "",
                "to": "",
            }
        )

    comparisons: list[dict[str, object]] = []
    for first in range(len(groups)):
        for second in range(first + 1, len(groups)):
            label_a = groups[first][0]
            label_b = groups[second][0]
            index_a = full_names.index(f"effect_{label_a}")
            index_b = full_names.index(f"effect_{label_b}")
            beta = float(fitted.params[index_b] - fitted.params[index_a])
            se = float(
                np.sqrt(covariance[index_b, index_b] + covariance[index_a, index_a] - 2 * covariance[index_b, index_a])
            )
            comparisons.append(
                {
                    "band": band,
                    "result_type": "interaction_difference",
                    "label": f"{label_b}_vs_{label_a}",
                    "log_ratio": beta,
                    "ratio": float(np.exp(beta)),
                    "se": se,
                    "low": float(np.exp(beta - critical * se)),
                    "high": float(np.exp(beta + critical * se)),
                    "p_value": float(2 * stats.t.sf(abs(beta / se), df)),
                    "df": df,
                    "years_in_era": np.nan,
                    "from": label_a,
                    "to": label_b,
                }
            )
    return pd.DataFrame(effects), pd.DataFrame(comparisons)


def policy_era_interaction(counts: pd.DataFrame | None = None) -> pd.DataFrame:
    """Fit the 2001-2006, 2007-2019 and 2020-2024 effects jointly."""

    counts = national_counts() if counts is None else counts
    groups: tuple[tuple[str, Callable[[int], bool]], ...] = (
        ("policy_2001_2006", lambda year: year <= 2006),
        ("policy_2007_2019", lambda year: 2007 <= year <= 2019),
        ("pandemic_recovery_2020_2024", lambda year: year >= 2020),
    )
    rows: list[pd.DataFrame] = []
    for band in BANDS:
        effects, comparisons = _era_model(counts, band=band, groups=groups)
        rows.extend([effects, comparisons])
    return pd.concat(rows, ignore_index=True)


def period_interactions(counts: pd.DataFrame | None = None) -> pd.DataFrame:
    """Fit the pandemic/recovery and statutory-calendar period interactions."""

    counts = national_counts() if counts is None else counts
    specifications: tuple[tuple[str, tuple[tuple[str, Callable[[int], bool]], ...]], ...] = (
        (
            "three_eras",
            (
                ("pre_2001_2019", lambda year: year <= 2019),
                ("covid_2020_2021", lambda year: 2020 <= year <= 2021),
                ("recovery_2022_2024", lambda year: year >= 2022),
            ),
        ),
        (
            "policy_eras",
            (
                ("policy_2001_2006", lambda year: year <= 2006),
                ("policy_2007_2019", lambda year: 2007 <= year <= 2019),
                ("pandemic_recovery_2020_2024", lambda year: year >= 2020),
            ),
        ),
    )
    rows: list[pd.DataFrame] = []
    for model, groups in specifications:
        for band in BANDS:
            effects, comparisons = _era_model(counts, band=band, groups=groups)
            for frame in (effects, comparisons):
                exported = frame.copy()
                exported.insert(0, "model", model)
                rows.append(exported)
    return pd.concat(rows, ignore_index=True)


def arizona_negative_control_support(counts: pd.DataFrame | None = None) -> pd.DataFrame:
    """Report support and model status for the no-switch Arizona diagnostic."""

    counts = national_counts() if counts is None else counts
    rows: list[dict[str, object]] = []
    for endpoint in (PEDESTRIAN, OTHER_FATAL, "all_fatal_crashes"):
        for band in BANDS:
            subset = counts[
                counts.geo.eq("AZ_Phoenix")
                & counts.width_minutes.eq(40)
                & counts.endpoint.eq(endpoint)
                & counts.band.isin([band, "night"])
            ]
            table = subset.pivot(
                index=["year", "season", "sham"],
                columns=["period", "band"],
                values="events",
            ).sort_index()
            a = table[("post", band)]
            b = table[("post", "night")]
            c = table[("pre", band)]
            d = table[("pre", "night")]
            total = a + b + c + d
            focal = a + c
            post = a + b
            informative = np.maximum(0, post - (total - focal)) < np.minimum(post, focal)
            row: dict[str, object] = {
                "endpoint": endpoint,
                "band": band,
                "focal_events": int(subset[subset.band.eq(band)].events.sum()),
                "reference_events": int(subset[subset.band.eq("night")].events.sum()),
                "candidate_tables": len(table),
                "informative_tables": int(informative.sum()),
                "all_zero_tables": int(total.eq(0).sum()),
            }
            for corrected, label in ((True, "three_center_adjusted"), (False, "reduced_actual_date_only")):
                try:
                    fitted = estimate(
                        counts,
                        band=band,
                        endpoint=endpoint,
                        geo="AZ_Phoenix",
                        corrected=corrected,
                    )
                    row[f"{label}_status"] = "fit"
                    row[f"{label}_ratio"] = fitted["ratio"]
                    row[f"{label}_low"] = fitted["low"]
                    row[f"{label}_high"] = fitted["high"]
                except Exception as error:  # noqa: BLE001 - failures are the diagnostic
                    row[f"{label}_status"] = f"unsupported: {type(error).__name__}: {error}"
                    row[f"{label}_ratio"] = np.nan
                    row[f"{label}_low"] = np.nan
                    row[f"{label}_high"] = np.nan
            rows.append(row)
    return pd.DataFrame(rows)


def alaska_inclusion_sensitivity(counts: pd.DataFrame | None = None) -> pd.DataFrame:
    """Compare the prespecified contiguous-US target with an Alaska-inclusive fit."""

    counts = national_counts() if counts is None else counts
    alaska = sensitivity_table("alaska_counts.csv.gz")
    keys = [
        "year",
        "season",
        "sham",
        "period",
        "width_minutes",
        "band",
        "endpoint",
        "window_days",
        "slot_minutes",
    ]
    main = counts[counts.geo.eq("DST_main")][keys + ["events"]]
    combined = (
        pd.concat([main, alaska[keys + ["events"]]], ignore_index=True)
        .groupby(keys, as_index=False, sort=True)
        .events.sum()
    )
    combined.insert(0, "geo", "DST_main_plus_Alaska")
    expanded = pd.concat([counts, combined[counts.columns]], ignore_index=True)

    rows: list[dict[str, object]] = []
    for geography in ("DST_main", "DST_main_plus_Alaska"):
        for band in BANDS:
            fitted: dict[str, dict[str, object]] = {}
            for endpoint in (PEDESTRIAN, OTHER_FATAL, "all_fatal_crashes"):
                result = estimate(expanded, band=band, endpoint=endpoint, geo=geography)
                fitted[endpoint] = result
                rows.append(
                    {
                        "geography": geography,
                        "endpoint": endpoint,
                        "band": band,
                        "ratio": result["ratio"],
                        "low": result["low"],
                        "high": result["high"],
                        "se": result["se"],
                        "df": result["df"],
                        "informative_strata": result["informative_strata"],
                        "raw_strata": result["raw_strata"],
                    }
                )
            direct = contrast(
                fitted[PEDESTRIAN],
                fitted[OTHER_FATAL],
                f"{geography}_{band}",
            )
            rows.append(
                {
                    "geography": geography,
                    "endpoint": "direct",
                    "band": band,
                    "ratio": direct["ratio"],
                    "low": direct["low"],
                    "high": direct["high"],
                    "se": direct["se"],
                    "df": direct["df"],
                    "informative_strata": float("nan"),
                    "raw_strata": float("nan"),
                }
            )
    return pd.DataFrame(rows)


def review_diagnostics(counts: pd.DataFrame | None = None) -> dict[str, pd.DataFrame]:
    counts = national_counts() if counts is None else counts
    return {
        "season_heterogeneity": seasonal_heterogeneity(counts),
        "leave_one_year": leave_one_year_influence(counts),
        "policy_era_interaction": policy_era_interaction(counts),
        "period_interactions": period_interactions(counts),
        "arizona_negative_control": arizona_negative_control_support(counts),
        "alaska_inclusion": alaska_inclusion_sensitivity(counts),
    }
