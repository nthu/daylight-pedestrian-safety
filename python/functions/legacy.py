"""Real-transition-only event-week diagnostic used in Supplementary Fig. S5."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from .conditional_poisson import fit
from .data import sensitivity_table


def estimate_legacy_weeks() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Fit all 28 Monday-Saturday event-week coefficients.

    Week -1 is the reference and is therefore not estimated.  This diagnostic
    is distinct from the sham-adjusted dynamic analysis.
    """

    counts = sensitivity_table("legacy_weeks.csv")
    if len(counts) != 1152 or not counts.events.ge(0).all():
        raise ValueError("Unexpected legacy-week input")

    estimates: list[dict[str, object]] = []
    pretests: list[dict[str, object]] = []
    for band in ("evening", "morning"):
        frame = (
            counts[counts.band.isin([band, "night"])]
            .sort_values(["year", "season", "week", "band"])
            .reset_index(drop=True)
        )
        units = frame.year.astype(str) + "|" + frame.season
        unit_codes = pd.Categorical(units)
        target = frame.band.eq(band).to_numpy(float)

        design = np.eye(len(unit_codes.categories))[unit_codes.codes] * target[:, None]
        names = [f"baseline_{category}" for category in unit_codes.categories]
        extra = []
        for season in ("spring", "autumn"):
            for week in (-4, -3, -2, 0, 1, 2, 3):
                extra.append((frame.season.eq(season) & frame.week.eq(week)).to_numpy(float) * target)
                names.append(f"{season}_week_{week}")
        design = np.column_stack([design, *extra])
        groups = units + "|" + frame.week.astype(str)

        fitted = fit(
            frame.events,
            design,
            groups,
            np.repeat("national", len(frame)),
            frame.year.astype(str),
        )
        for row in fitted["rows"]:
            term = names[row["term_index"]]
            if row["method"] != "year" or not term.startswith(("spring_", "autumn_")):
                continue
            estimate = dict(row)
            estimate.update(
                {
                    "band": band,
                    "term": term,
                    "season": term.split("_week_")[0],
                    "week": int(term.split("_week_")[1]),
                }
            )
            estimates.append(estimate)

        for season in ("spring", "autumn"):
            indices = [names.index(f"{season}_week_{week}") for week in (-4, -3, -2)]
            beta = fitted["beta"][indices]
            covariance = fitted["cov"]["year"][np.ix_(indices, indices)]
            if np.linalg.matrix_rank(covariance) != 3:
                raise ValueError("Pre-period covariance is not full rank")
            statistic = float(beta @ np.linalg.solve(covariance, beta) / 3)
            pretests.append(
                {
                    "band": band,
                    "season": season,
                    "F": statistic,
                    "df_num": 3,
                    "df_denom": 23,
                    "p_value": float(stats.f.sf(statistic, 3, 23)),
                }
            )

    return pd.DataFrame(estimates), pd.DataFrame(pretests)
