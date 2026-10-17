"""Stratified national and Mississippi registry estimators."""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np
import pandas as pd


def mantel_haenszel_log_contrast(
    frame: pd.DataFrame,
    *,
    unit: str,
    band: str,
    endpoint: str,
    years: Iterable[int],
) -> tuple[float, pd.DataFrame]:
    """Return the gain-oriented symmetric log contrast within fixed units.

    The input must contain zero-inclusive counts for a focal band and the
    midnight reference.  The fixed unit is either a state or county.
    """

    selected = frame[frame.band.isin([band, "night"]) & frame.endpoint.eq(endpoint) & frame.year.isin(list(years))]
    index = [unit, "year", "season", "sham"]
    table = selected.pivot(index=index, columns=["period", "band"], values="events")

    post_focal = table["post", band].astype(float)
    post_reference = table["post", "night"].astype(float)
    pre_focal = table["pre", band].astype(float)
    pre_reference = table["pre", "night"].astype(float)
    total = post_focal + post_reference + pre_focal + pre_reference

    numerator = (post_focal * pre_reference / total).fillna(0)
    denominator = (post_reference * pre_focal / total).fillna(0)
    components = pd.DataFrame({"numerator": numerator, "denominator": denominator})
    sums = components.groupby(level=["season", "sham"]).sum()
    if (sums <= 0).any().any():
        raise ValueError(f"Zero aggregate support for {endpoint}, {band}")

    log_odds = np.log(sums.numerator / sums.denominator)
    seasonal = []
    for season in ("spring", "autumn"):
        direction = (1 if season == "spring" else -1) * (1 if band == "evening" else -1)
        seasonal.append(direction * (log_odds[season, 0] - 0.5 * (log_odds[season, -35] + log_odds[season, 35])))
    return float(np.mean(seasonal)), table


def information_support(table: pd.DataFrame, band: str) -> pd.DataFrame:
    """Describe appearances and fixed-margin information for each table."""

    a = table["post", band]
    b = table["post", "night"]
    c = table["pre", band]
    d = table["pre", "night"]
    total = a + b + c + d
    post_total = a + b
    focal_total = a + c
    informative = np.maximum(0, post_total - (total - focal_total)) < np.minimum(post_total, focal_total)

    support = table.index.to_frame(index=False)
    support["appearances"] = total.to_numpy()
    support["informative"] = informative.to_numpy()
    support["all_zero"] = total.eq(0).to_numpy()
    return support
