"""Five-center curvature and policy-year multiplier-bootstrap diagnostics."""

from __future__ import annotations

from collections.abc import Iterable, Mapping

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

from .data import national_counts, sensitivity_table

PEDESTRIAN = "fatal_pedestrian_crashes"
OTHER_FATAL = "other_fatal_crashes"
BANDS = ("evening", "morning")
ENDPOINTS = (PEDESTRIAN, OTHER_FATAL)


def five_center_counts() -> pd.DataFrame:
    """Assemble k=-70,-35,0,+35,+70 counts on the common 2002-2023 support.

    The two outer centers are the actual-center rows of the archived -70 and
    +70 placebo calendars.  Their overlapping +/-35 rows are required to match
    the corresponding rows in the main calendar exactly before fitting.
    """

    base = national_counts()
    early = sensitivity_table("counts_placebo_-70.csv")
    late = sensitivity_table("counts_placebo_+70.csv")
    keys = ["year", "season", "period", "width_minutes", "band", "endpoint"]

    main = base[base.geo.eq("DST_main") & base.year.between(2002, 2023) & base.sham.isin([-35, 0, 35])][
        ["geo", *keys, "sham", "events"]
    ].copy()

    for label, frame, placebo_center, main_center in (
        ("early", early, 35, -35),
        ("late", late, -35, 35),
    ):
        placebo_overlap = frame.loc[frame.sham.eq(placebo_center), [*keys, "events"]]
        main_overlap = main.loc[main.sham.eq(main_center), [*keys, "events"]]
        merged = placebo_overlap.merge(
            main_overlap,
            on=keys,
            suffixes=("_placebo", "_main"),
            validate="one_to_one",
        )
        if len(merged) != len(placebo_overlap) or not np.array_equal(
            merged.events_placebo.to_numpy(), merged.events_main.to_numpy()
        ):
            raise AssertionError(f"{label} placebo calendar does not recover the shared main-calendar center")

    outer_early = early.loc[early.sham.eq(0), ["geo", *keys, "events"]].copy()
    outer_early["sham"] = -70
    outer_late = late.loc[late.sham.eq(0), ["geo", *keys, "events"]].copy()
    outer_late["sham"] = 70
    columns = ["geo", "year", "season", "sham", "period", "width_minutes", "band", "endpoint", "events"]
    combined = pd.concat([outer_early[columns], main[columns], outer_late[columns]], ignore_index=True)
    identity = ["geo", "year", "season", "sham", "period", "width_minutes", "band", "endpoint"]
    if combined.duplicated(identity).any():
        raise AssertionError("Duplicate five-center count identity")
    if sorted(combined.sham.unique()) != [-70, -35, 0, 35, 70]:
        raise AssertionError("Five-center support is incomplete")
    return combined.sort_values(identity).reset_index(drop=True)


def _fit_grouped_binomial(
    counts: pd.DataFrame,
    *,
    band: str,
    endpoint: str,
    quadratic: bool,
    years: Iterable[int] | None = None,
) -> dict[str, object]:
    selected = counts[
        counts.geo.eq("DST_main")
        & counts.width_minutes.eq(40)
        & counts.endpoint.eq(endpoint)
        & counts.band.isin([band, "night"])
    ].copy()
    if years is not None:
        selected = selected[selected.year.isin(list(years))]

    table = selected.pivot(
        index=["year", "season", "sham"],
        columns=["period", "band"],
        values="events",
    ).sort_index()
    a = table[("post", band)]
    b = table[("post", "night")]
    c = table[("pre", band)]
    d = table[("pre", "night")]
    total_2x2 = a + b + c + d
    focal_total = a + c
    post_total = a + b
    informative = np.maximum(0, post_total - (total_2x2 - focal_total)) < np.minimum(post_total, focal_total)
    zero_information = table.index[~informative].tolist()
    table = table[informative]
    meta = table.index.to_frame(index=False)
    strata = len(meta)
    if strata < 10:
        raise ValueError("Insufficient informative five-center strata")

    nuisance: list[np.ndarray] = []
    names: list[str] = []
    for year in sorted(meta.year.unique()):
        for season in ("autumn", "spring"):
            nuisance.append((meta.year.eq(year) & meta.season.eq(season)).to_numpy(float))
            names.append(f"post_{year}_{season}")

    center = meta.sham.to_numpy(float) / 35.0
    for season in ("autumn", "spring"):
        mask = meta.season.eq(season).to_numpy(float)
        nuisance.append(mask * center)
        names.append(f"{season}_linear")
        if quadratic:
            nuisance.append(mask * center**2)
            names.append(f"{season}_quadratic")

    direction = np.where(meta.season.eq("spring"), 1.0, -1.0) * (1.0 if band == "evening" else -1.0)
    nuisance.append(meta.sham.eq(0).to_numpy(float) * direction)
    names.append("light_gain")
    extras = np.column_stack(nuisance)

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
    model = sm.GLM(np.column_stack([successes, failures]), X, family=sm.families.Binomial())
    fitted = model.fit(maxiter=500, tol=1e-11)
    if not fitted.converged or not np.isfinite(fitted.params).all():
        raise RuntimeError("Five-center grouped-binomial model failed to converge")

    totals = successes + failures
    probabilities = fitted.fittedvalues
    score_rows = X * (successes - totals * probabilities)[:, None]
    information = X.T @ (X * (totals * probabilities * (1.0 - probabilities))[:, None])
    eigenvalues = np.linalg.eigvalsh((information + information.T) / 2)
    if eigenvalues.min() <= eigenvalues.max() * 1e-12:
        raise ValueError("Five-center information matrix is not full rank")
    bread = np.linalg.inv(information)
    clusters = np.unique(year)
    cluster_scores = np.vstack([score_rows[year == cluster].sum(axis=0) for cluster in clusters])
    G = len(clusters)
    N = len(successes)
    K = X.shape[1]
    correction = G / (G - 1) * (N - 1) / (N - K)
    covariance = bread @ cluster_scores.T @ cluster_scores @ bread * correction

    full_names = [f"stratum_{index}" for index in range(strata)] + names
    target = full_names.index("light_gain")
    beta = float(fitted.params[target])
    standard_error = float(np.sqrt(covariance[target, target]))
    degrees_freedom = G - 1
    critical = stats.t.ppf(0.975, degrees_freedom)
    influence = cluster_scores @ bread[:, target] * np.sqrt(correction)
    if not np.isclose(float(influence @ influence), standard_error**2, rtol=1e-8, atol=1e-12):
        raise AssertionError("Five-center policy-year influence decomposition mismatch")

    quadratic_indices = (
        [full_names.index(f"{season}_quadratic") for season in ("autumn", "spring")] if quadratic else []
    )
    if quadratic_indices:
        coefficients = fitted.params[quadratic_indices]
        variance = covariance[np.ix_(quadratic_indices, quadratic_indices)]
        wald = float(coefficients @ np.linalg.solve(variance, coefficients))
        quadratic_p = float(stats.f.sf(wald / len(quadratic_indices), len(quadratic_indices), degrees_freedom))
    else:
        wald = np.nan
        quadratic_p = np.nan

    return {
        "band": band,
        "endpoint": endpoint,
        "model": "five_center_quadratic" if quadratic else "five_center_linear",
        "log_ratio": beta,
        "ratio": float(np.exp(beta)),
        "se": standard_error,
        "low": float(np.exp(beta - critical * standard_error)),
        "high": float(np.exp(beta + critical * standard_error)),
        "p_value": float(2 * stats.t.sf(abs(beta / standard_error), degrees_freedom)),
        "df": degrees_freedom,
        "clusters": G,
        "allocation_observations": N,
        "parameters": K,
        "informative_strata": strata,
        "zero_information_strata": len(zero_information),
        "cr1_factor": float(correction),
        "information_condition": float(np.linalg.cond(information)),
        "quadratic_wald": wald,
        "quadratic_df_num": 2 if quadratic else np.nan,
        "quadratic_df_denom": degrees_freedom if quadratic else np.nan,
        "quadratic_p_value": quadratic_p,
        "autumn_quadratic": float(fitted.params[full_names.index("autumn_quadratic")]) if quadratic else np.nan,
        "spring_quadratic": float(fitted.params[full_names.index("spring_quadratic")]) if quadratic else np.nan,
        "year_influence": {int(cluster): float(value) for cluster, value in zip(clusters, influence, strict=True)},
    }


def five_center_analysis(full: bool = True) -> dict[str, pd.DataFrame]:
    """Fit linear and quadratic five-center models plus paired year deletions."""

    counts = five_center_counts()
    estimates: list[dict[str, object]] = []
    deletions: list[dict[str, object]] = []
    fitted: dict[tuple[str, str, str], dict[str, object]] = {}

    for quadratic in (False, True):
        model = "five_center_quadratic" if quadratic else "five_center_linear"
        for endpoint in ENDPOINTS:
            for band in BANDS:
                result = _fit_grouped_binomial(counts, band=band, endpoint=endpoint, quadratic=quadratic)
                fitted[model, endpoint, band] = result
                estimates.append({key: value for key, value in result.items() if key != "year_influence"})

        for band in BANDS:
            pedestrian = fitted[model, PEDESTRIAN, band]
            other = fitted[model, OTHER_FATAL, band]
            years = sorted(set(pedestrian["year_influence"]) & set(other["year_influence"]))
            contrast_influence = np.asarray(
                [pedestrian["year_influence"][year] - other["year_influence"][year] for year in years]
            )
            beta = pedestrian["log_ratio"] - other["log_ratio"]
            standard_error = float(np.sqrt(contrast_influence @ contrast_influence))
            df = len(years) - 1
            critical = stats.t.ppf(0.975, df)
            estimates.append(
                {
                    "band": band,
                    "endpoint": "direct",
                    "model": model,
                    "log_ratio": beta,
                    "ratio": float(np.exp(beta)),
                    "se": standard_error,
                    "low": float(np.exp(beta - critical * standard_error)),
                    "high": float(np.exp(beta + critical * standard_error)),
                    "p_value": float(2 * stats.t.sf(abs(beta / standard_error), df)),
                    "df": df,
                    "clusters": len(years),
                    "allocation_observations": np.nan,
                    "parameters": np.nan,
                    "informative_strata": np.nan,
                    "zero_information_strata": np.nan,
                    "cr1_factor": np.nan,
                    "information_condition": np.nan,
                    "quadratic_wald": np.nan,
                    "quadratic_df_num": np.nan,
                    "quadratic_df_denom": np.nan,
                    "quadratic_p_value": np.nan,
                    "autumn_quadratic": np.nan,
                    "spring_quadratic": np.nan,
                }
            )

    if full:
        quadratic = True
        for endpoint in ENDPOINTS:
            for band in BANDS:
                full_result = fitted["five_center_quadratic", endpoint, band]
                deleted = []
                for omitted in range(2002, 2024):
                    result = _fit_grouped_binomial(
                        counts,
                        band=band,
                        endpoint=endpoint,
                        quadratic=quadratic,
                        years=[year for year in range(2002, 2024) if year != omitted],
                    )
                    deleted.append(result["log_ratio"])
                    deletions.append(
                        {
                            "band": band,
                            "endpoint": endpoint,
                            "omitted_policy_year": omitted,
                            "log_ratio": result["log_ratio"],
                            "ratio": result["ratio"],
                        }
                    )
                values = np.asarray(deleted)
                blocks = len(values)
                se = float(np.sqrt((blocks - 1) / blocks * np.sum((values - values.mean()) ** 2)))
                critical = stats.t.ppf(0.975, blocks - 1)
                estimates.append(
                    {
                        "band": band,
                        "endpoint": endpoint,
                        "model": "five_center_quadratic_year_jackknife",
                        "log_ratio": full_result["log_ratio"],
                        "ratio": full_result["ratio"],
                        "se": se,
                        "low": float(np.exp(full_result["log_ratio"] - critical * se)),
                        "high": float(np.exp(full_result["log_ratio"] + critical * se)),
                        "p_value": float(2 * stats.t.sf(abs(full_result["log_ratio"] / se), blocks - 1)),
                        "df": blocks - 1,
                        "clusters": blocks,
                    }
                )

        for band in BANDS:
            pedestrian_full = fitted["five_center_quadratic", PEDESTRIAN, band]
            other_full = fitted["five_center_quadratic", OTHER_FATAL, band]
            ped = (
                pd.DataFrame(deletions)
                .query("band == @band and endpoint == @PEDESTRIAN")
                .sort_values("omitted_policy_year")
            )
            oth = (
                pd.DataFrame(deletions)
                .query("band == @band and endpoint == @OTHER_FATAL")
                .sort_values("omitted_policy_year")
            )
            values = ped.log_ratio.to_numpy() - oth.log_ratio.to_numpy()
            blocks = len(values)
            beta = pedestrian_full["log_ratio"] - other_full["log_ratio"]
            se = float(np.sqrt((blocks - 1) / blocks * np.sum((values - values.mean()) ** 2)))
            critical = stats.t.ppf(0.975, blocks - 1)
            estimates.append(
                {
                    "band": band,
                    "endpoint": "direct",
                    "model": "five_center_quadratic_year_jackknife",
                    "log_ratio": beta,
                    "ratio": float(np.exp(beta)),
                    "se": se,
                    "low": float(np.exp(beta - critical * se)),
                    "high": float(np.exp(beta + critical * se)),
                    "p_value": float(2 * stats.t.sf(abs(beta / se), blocks - 1)),
                    "df": blocks - 1,
                    "clusters": blocks,
                }
            )

    return {
        "five_center_curvature": pd.DataFrame(estimates),
        "five_center_deletions": pd.DataFrame(deletions),
    }


def exact_rademacher_multiplier(
    log_ratio: float,
    influences: Mapping[int | str, float] | Iterable[float],
    *,
    chunk_size: int = 1 << 18,
) -> dict[str, float | int]:
    """Exact two-sided Rademacher wild-cluster multiplier bootstrap-t.

    The input values are the policy-year influence contributions whose squared
    sum equals the fitted CR1 variance.  A shared sign vector can therefore be
    used for endpoint contrasts without discarding their covariance.  Global
    sign reversal leaves the absolute statistic unchanged, so only 2^(G-1)
    unique sign patterns are enumerated.
    """

    values = np.asarray(list(influences.values()) if isinstance(influences, Mapping) else list(influences), dtype=float)
    clusters = len(values)
    if clusters < 3 or not np.isfinite(values).all():
        raise ValueError("At least three finite policy-year influences are required")
    if abs(values.sum()) > 1e-8:
        raise AssertionError("Influence contributions must sum to zero at the fitted model")
    variance = float(values @ values)
    standard_error = float(np.sqrt(variance))
    observed_t = abs(log_ratio / standard_error)

    enumerated = 1 << (clusters - 1)
    bits = np.arange(1, clusters, dtype=np.uint64)
    statistics = np.empty(enumerated, dtype=np.float32)
    position = 0
    for start in range(0, enumerated, chunk_size):
        identifiers = np.arange(start, min(start + chunk_size, enumerated), dtype=np.uint64)
        signs = np.ones((len(identifiers), clusters), dtype=float)
        signs[:, 1:] = 1.0 - 2.0 * ((identifiers[:, None] >> (bits - 1)) & 1).astype(float)
        sums = signs @ values
        denominator = np.sqrt(np.maximum(variance - sums**2 / clusters, 1e-14))
        current = np.abs(sums / denominator)
        statistics[position : position + len(current)] = current
        position += len(current)

    critical = float(np.quantile(statistics, 0.95, method="higher"))
    p_value = float(np.mean(statistics >= observed_t - 1e-12))
    return {
        "wild_p_value": p_value,
        "wild_critical_95": critical,
        "wild_low": float(np.exp(log_ratio - critical * standard_error)),
        "wild_high": float(np.exp(log_ratio + critical * standard_error)),
        "observed_t": float(observed_t),
        "clusters": clusters,
        "unique_sign_patterns": enumerated,
        "equivalent_full_sign_patterns": 1 << clusters,
    }


def wild_cluster_calibration(counts: pd.DataFrame | None = None) -> pd.DataFrame:
    """Calibrate the four core comparisons with shared policy-year signs."""

    from .allocation import estimate

    counts = national_counts() if counts is None else counts
    fitted: dict[tuple[str, str], dict[str, object]] = {}
    rows: list[dict[str, object]] = []
    for endpoint in ENDPOINTS:
        for band in BANDS:
            fitted[endpoint, band] = estimate(counts, band=band, endpoint=endpoint)

    for band in BANDS:
        pedestrian = fitted[PEDESTRIAN, band]
        calibrated = exact_rademacher_multiplier(pedestrian["log_ratio"], pedestrian["year_influence"])
        rows.append(
            {
                "comparison": f"pedestrian_{band}",
                "band": band,
                "endpoint": PEDESTRIAN,
                "ratio": pedestrian["ratio"],
                "cr1_low": pedestrian["low"],
                "cr1_high": pedestrian["high"],
                "cr1_p_value": pedestrian["p_value"],
                **calibrated,
            }
        )

        other = fitted[OTHER_FATAL, band]
        years = sorted(set(pedestrian["year_influence"]) & set(other["year_influence"]))
        influence = {year: pedestrian["year_influence"][year] - other["year_influence"][year] for year in years}
        beta = pedestrian["log_ratio"] - other["log_ratio"]
        standard_error = float(np.sqrt(sum(value * value for value in influence.values())))
        df = len(years) - 1
        cr1_p = float(2 * stats.t.sf(abs(beta / standard_error), df))
        cr1_q = stats.t.ppf(0.975, df)
        calibrated = exact_rademacher_multiplier(beta, influence)
        rows.append(
            {
                "comparison": f"pedestrian_vs_other_{band}",
                "band": band,
                "endpoint": "direct",
                "ratio": float(np.exp(beta)),
                "cr1_low": float(np.exp(beta - cr1_q * standard_error)),
                "cr1_high": float(np.exp(beta + cr1_q * standard_error)),
                "cr1_p_value": cr1_p,
                **calibrated,
            }
        )
    return pd.DataFrame(rows)
