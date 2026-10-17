"""Small-sample inference helpers used across the analysis scripts."""

from __future__ import annotations

from collections.abc import Callable, Iterable, Sequence

import numpy as np
from scipy import stats


def jackknife_interval(log_ratio: float, deleted: Sequence[float]) -> dict[str, float | int]:
    """Return a paired deletion-jackknife interval on the ratio scale."""

    values = np.asarray(deleted, dtype=float)
    blocks = len(values)
    if blocks < 3 or not np.isfinite(values).all():
        raise ValueError("At least three finite deletion estimates are required")

    variance = (blocks - 1) / blocks * np.sum((values - values.mean()) ** 2)
    standard_error = float(np.sqrt(variance))
    if standard_error <= 0:
        raise ValueError("Deletion variance is zero")

    degrees_freedom = blocks - 1
    critical = stats.t.ppf(0.975, degrees_freedom)
    return {
        "log_ratio": float(log_ratio),
        "ratio": float(np.exp(log_ratio)),
        "se": standard_error,
        "low": float(np.exp(log_ratio - critical * standard_error)),
        "high": float(np.exp(log_ratio + critical * standard_error)),
        "p_value": float(2 * stats.t.sf(abs(log_ratio / standard_error), degrees_freedom)),
        "df": degrees_freedom,
        "blocks": blocks,
    }


def delete_one_year(
    years: Iterable[int],
    fit: Callable[[list[int]], float],
) -> tuple[list[int], np.ndarray]:
    """Fit one model per omitted policy year and return log-ratio estimates."""

    ordered = sorted(int(year) for year in years)
    estimates = np.asarray([fit([year for year in ordered if year != omitted]) for omitted in ordered], dtype=float)
    if not np.isfinite(estimates).all():
        raise ValueError("A deletion fit returned a non-finite estimate")
    return ordered, estimates


def t_interval(log_ratio: float, standard_error: float, degrees_freedom: int) -> dict[str, float | int]:
    """Return a Student-t interval for a log ratio."""

    critical = stats.t.ppf(0.975, degrees_freedom)
    return {
        "log_ratio": float(log_ratio),
        "ratio": float(np.exp(log_ratio)),
        "se": float(standard_error),
        "low": float(np.exp(log_ratio - critical * standard_error)),
        "high": float(np.exp(log_ratio + critical * standard_error)),
        "p_value": float(2 * stats.t.sf(abs(log_ratio / standard_error), degrees_freedom)),
        "df": int(degrees_freedom),
    }
