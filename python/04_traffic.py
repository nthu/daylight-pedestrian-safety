"""Summarize matched FHWA traffic changes and site-support sensitivities."""

from __future__ import annotations

import numpy as np
import pandas as pd
from functions.data import ensure_results, traffic_counts, write_csv
from functions.specifications import BANDS
from scipy import stats


def run() -> pd.DataFrame:
    traffic = traffic_counts()
    rows: list[dict[str, object]] = []

    supported_metrics = ("log_day_level", "log_focal_daily_share", "log_focal_vs_midnight")
    for metric in supported_metrics:
        for band in BANDS:
            for minimum_sites in (1, 2, 5, 10, 20, 30):
                selected = traffic[
                    traffic.metric.eq(metric) & traffic.band.eq(band) & traffic.n_sites.ge(minimum_sites)
                ]
                values = selected.shock_log_ratio_2020_vs_2019.to_numpy()
                jurisdictions = len(selected)
                if jurisdictions < 2:
                    continue

                mean = float(values.mean())
                standard_error = float(values.std(ddof=1) / np.sqrt(jurisdictions))
                critical = stats.t.ppf(0.975, jurisdictions - 1)
                rows.append(
                    {
                        "metric": metric,
                        "band": band,
                        "minimum_sites": minimum_sites,
                        "jurisdictions": jurisdictions,
                        "sites": int(selected.n_sites.sum()),
                        "ratio_equal_jurisdiction": float(np.exp(mean)),
                        "low": float(np.exp(mean - critical * standard_error)),
                        "high": float(np.exp(mean + critical * standard_error)),
                        "ratio_median_jurisdiction": float(np.exp(np.median(values))),
                        "ratio_site_weighted": float(np.exp(np.average(values, weights=selected.n_sites))),
                    }
                )

    results = pd.DataFrame(rows)
    paths = ensure_results()
    write_csv(results, paths["estimates"] / "traffic_support.csv")

    print("\nMatched-traffic sensitivity: daily volume, evening")
    print(
        results.query("metric == 'log_day_level' and band == 'evening'")[
            ["minimum_sites", "jurisdictions", "sites", "ratio_equal_jurisdiction", "low", "high"]
        ].to_string(index=False)
    )
    return results


if __name__ == "__main__":
    run()
