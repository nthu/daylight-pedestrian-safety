"""Rebuild all ten submission figures from recomputed results and declared source metadata.
Run after `python reproduce.py --full`. No expected estimate is used for drawing.
Each figure has its own axes, vector PDF and 300-dpi PNG; numerical plot inputs
are exported next to the figures. Source-derived and fitted figures are labelled
in FIGURE_DATA_MANIFEST.json. No individual event data are required.
"""

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter

from .data import RESULTS, national_counts, support_table

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["font.size"] = 10


def run(outdir=None):
    root = Path(__file__).resolve().parents[2]
    o = RESULTS / "estimates"
    out = Path(outdir) if outdir else RESULTS / "figures"
    out.mkdir(parents=True, exist_ok=True)
    meta = {}

    def sha(path):
        return hashlib.sha256(Path(path).read_bytes()).hexdigest()

    def resolve_source(label):
        path = Path(label)
        if path.parts[:2] == ("results", "python"):
            return RESULTS.joinpath(*path.parts[2:])
        return root / path

    def read(name):
        p = o / name
        if not p.exists():
            raise FileNotFoundError(f"Run reproduce.py --full first: {p}")
        return pd.read_csv(p)

    def save(fig, name, data, sources, kind):
        if {"low", "high"}.issubset(data.columns):
            axis = fig.axes[0]
            lo, hi = (
                axis.get_ylim()
                if name in ["Figure_2_aligned_dynamics", "Figure_S5_original_event_weeks"]
                else axis.get_xlim()
            )
            assert data.low.min() >= lo and data.high.max() <= hi, (name, "interval clipped", lo, hi)
        fig.savefig(
            out / (name + ".pdf"),
            bbox_inches="tight",
            metadata={"Title": name, "Author": "", "Creator": "USFARS reproducible figure code"},
        )
        fig.savefig(out / (name + ".png"), dpi=300, bbox_inches="tight")
        plt.close(fig)
        data.to_csv(out / (name + ".csv"), index=False)
        meta[name] = {
            "rows": len(data),
            "data_sha256": sha(out / (name + ".csv")),
            "source_kind": kind,
            "sources": {p: sha(resolve_source(p)) for p in sources},
            "axes_scale": "explicit in caption/axis",
            "figure_pdf": name + ".pdf",
            "figure_png": name + ".png",
        }

    def clean(ax):
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="x", alpha=0.16)

    def ratio_axis(ax, lo, hi, ticks):
        ax.set_xscale("log")
        ax.set_xlim(lo, hi)
        ax.set_xticks(ticks)
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, p: f"{v:g}"))
        ax.minorticks_off()
        ax.axvline(1, ls=":", lw=1.2)
        clean(ax)

    def forest(
        name,
        rows,
        sources,
        groups=None,
        xlim=(0.04, 2.0),
        figheight=4.6,
        xlabel="Ratio of count ratios (log scale)",
        kind="recomputed model estimates",
    ):
        fig, ax = plt.subplots(figsize=(10.2, figheight))
        fig.subplots_adjust(left=0.40, right=0.72, bottom=0.16, top=0.93)
        ys = np.arange(len(rows), 0, -1, dtype=float)
        if groups:
            ys = np.array(groups, dtype=float)
        r = pd.DataFrame(rows)
        r["display_y"] = ys
        ax.errorbar(r.ratio, ys, xerr=np.array([r.ratio - r.low, r.high - r.ratio]), fmt="o", capsize=3, ms=5, lw=1.35)
        ratio_axis(ax, *xlim, [t for t in [0.05, 0.1, 0.2, 0.5, 1, 2] if xlim[0] <= t <= xlim[1]])
        ax.set_ylim(min(ys) - 0.65, max(ys) + 0.95)
        ax.set_yticks(ys, r.label, fontsize=10)
        ax.tick_params(axis="y", length=0, pad=10)
        ax.set_xlabel(xlabel)
        ax.text(1.06, 1.015, "Estimate (95% CI)", transform=ax.transAxes, fontsize=10, weight="bold")
        for y, t in zip(ys, r.itertuples(), strict=False):
            ax.text(
                1.06,
                y,
                f"{t.ratio:.3f} ({t.low:.3f}–{t.high:.3f})",
                transform=ax.get_yaxis_transform(),
                va="center",
                fontsize=10,
            )
        save(fig, name, r, sources, kind)

    # Figure 1: a declared numerical example and exact inclusive date windows.
    fig, ax = plt.subplots(figsize=(9.3, 7.3))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")
    ax.text(2, 97, "A  Solar-time assignment: one sunset example", fontsize=12, weight="bold")
    ax.text(
        2,
        91,
        "Calculated sunset: 18:00 standard time. Recorded civil-clock target: 18:10–18:50.",
        fontsize=10,
    )
    for y, label, suntxt in [
        (79, "Standard time", "Sunset at 18:00 civil time"),
        (64, "Daylight-saving time", "Sunset at 19:00 civil time"),
    ]:
        ax.text(2, y + 4, label, fontsize=10, weight="bold")
        ax.hlines(y, 29, 95, lw=1)
        for val in [18, 18.5, 19]:
            x = 29 + (val - 17.5) / 2 * 66
            ax.vlines(x, y - 1, y + 1, lw=1)
            ax.text(x, y - 4, f"{int(val):02d}:{30 if val % 1 else 0:02d}", ha="center", fontsize=9)
        x1 = 29 + (18 + 10 / 60 - 17.5) / 2 * 66
        x2 = 29 + (18 + 50 / 60 - 17.5) / 2 * 66
        ax.add_patch(Rectangle((x1, y - 1.5), x2 - x1, 3, fill=False, hatch="////", lw=1.2))
        sun = 18 if label == "Standard time" else 19
        sx = 29 + (sun - 17.5) / 2 * 66
        ax.annotate(
            suntxt,
            xy=(sx, y),
            xytext=(sx, y + 7),
            ha="center",
            fontsize=9,
            arrowprops={"arrowstyle": "->", "lw": 0.8},
        )
        ax.text(97, y, "Below horizon" if sun == 18 else "Above horizon", ha="left", va="center", fontsize=9)
    ax.text(2, 52, "B  Matched calendar: days relative to the statutory change", fontsize=12, weight="bold")

    def xmap(day):
        return 24 + (day + 63) / 126 * 69

    calrows = []
    for y, label, center, pre, post in [
        (40, "Early sham", -35, (-56, -43), (-28, -15)),
        (29, "Actual switch", 0, (-21, -8), (7, 20)),
        (18, "Late sham", 35, (14, 27), (42, 55)),
    ]:
        ax.text(2, y, label, va="center", fontsize=10)
        ax.hlines(y, xmap(-63), xmap(63), lw=0.6, ls=":")
        for per, (a, b), hatch in [("Before", pre, "///"), ("After", post, "...")]:
            ax.add_patch(Rectangle((xmap(a), y - 2.3), xmap(b + 1) - xmap(a), 4.6, fill=False, hatch=hatch, lw=1))
            ax.text((xmap(a) + xmap(b + 1)) / 2, y + 4, f"{a:+d} to {b:+d}", ha="center", fontsize=9)
            calrows.append({"center_days": center, "period": per, "start_day": a, "end_day": b, "days": 14})
        ax.scatter([xmap(center)], [y], marker="D", s=16)
    ax.plot([xmap(0), xmap(0)], [11, 46], ls="--", lw=1)
    ax.text(xmap(0), 8, "Day 0", ha="center", fontsize=9)
    ax.text(
        2,
        2,
        "Each band: 40 minutes. Midnight reference: 00:10–00:50. Each before/after window: 14 days.",
        fontsize=9,
    )
    save(
        fig,
        "Figure_1_design",
        pd.DataFrame(calrows),
        ["data/calendar.csv"],
        "exact calendar and illustrative clock example",
    )
    # Figure 2: no line is interpolated through the omitted reference weeks.
    dyn = read("dynamic_estimates.csv")
    assert len(dyn) == 24
    fig, ax = plt.subplots(figsize=(9.2, 5.4))
    fig.subplots_adjust(left=0.10, right=0.98, top=0.85, bottom=0.15)
    ax.set_yscale("log")
    ax.set_ylim(0.008, 16)
    ax.set_xlim(-4.6, 3.6)
    ax.axhline(1, ls=":", lw=1.1)
    ax.axvline(-0.5, ls="--", lw=1)
    ax.axvspan(-3.45, -1.55, fill=False, hatch="////", alpha=0.22)
    ax.text(-2.5, 7.5, "Shared 14-day baseline\nweeks −3 and −2", ha="center", fontsize=9)
    series = [
        ("evening", "spring", "Spring evening", "o"),
        ("evening", "autumn", "Autumn evening", "s"),
        ("morning", "spring", "Spring morning", "^"),
        ("morning", "autumn", "Autumn morning", "D"),
    ]
    for j, (band, season, label, mark) in enumerate(series):
        z = dyn[dyn.band.eq(band) & dyn.season.eq(season)].sort_values("week")
        xx = z.week + (j - 1.5) * 0.095
        ax.errorbar(
            xx,
            z.ratio,
            yerr=np.array([z.ratio - z.low, z.high - z.ratio]),
            fmt=mark,
            capsize=2.4,
            ms=5,
            lw=1.1,
            label=label,
        )
    ax.set_xticks(range(-4, 4), ["−4", "−3", "−2", "−1", "0", "+1", "+2", "+3"])
    ax.set_yticks([0.01, 0.03, 0.1, 0.3, 1, 3, 10])
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, p: f"{v:g}"))
    ax.minorticks_off()
    ax.set_xlabel("Relative week at each real or sham center")
    ax.set_ylabel("Gain-oriented ratio of count ratios (log scale)")
    ax.legend(ncol=2, loc="lower center", bbox_to_anchor=(0.5, 1.01), frameon=False, fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", alpha=0.15)
    save(
        fig,
        "Figure_2_aligned_dynamics",
        dyn,
        ["results/python/estimates/dynamic_estimates.csv"],
        "24 re-estimated dynamic contrasts",
    )
    primary = read("national_models.csv").set_index("specification")
    direct = read("direct_contrasts.csv").set_index("contrast")
    rows = []
    for band in ["evening", "morning"]:
        for endpoint, spec in [
            ("Pedestrian", "primary_" + band),
            ("Other fatal", "other_fatal_crashes_" + band),
        ]:
            r = primary.loc[spec]
            rows.append({"label": f"{band.title()} · {endpoint} (RCR)", **{k: r[k] for k in ["ratio", "low", "high"]}})
    for band in ["evening", "morning"]:
        r = direct.loc["pedestrian_vs_other_" + band]
        rows.append({"label": f"{band.title()} · Pedestrian / other", **{k: r[k] for k in ["ratio", "low", "high"]}})
    forest(
        "Figure_3_primary_contrasts",
        rows,
        ["results/python/estimates/national_models.csv", "results/python/estimates/direct_contrasts.csv"],
        groups=[7, 6, 5, 4, 2, 1],
        figheight=5.1,
        xlabel="RCR, or direct ratio of RCRs (log scale)",
    )
    cond = read("conditional_blocks.csv")
    rows = []
    for band in ["evening", "morning"]:
        r = primary.loc["primary_" + band]
        rows.append(
            {
                "label": band.title() + " · Main allocation / CR1",
                **{k: r[k] for k in ["ratio", "low", "high"]},
            }
        )
        for block, text in [(1, "Fixed margins / year jackknife"), (4, "Fixed margins / 4-year blocks")]:
            r = cond[cond.band.eq(band) & cond.block_size_years.eq(block)].iloc[0]
            rows.append({"label": band.title() + " · " + text, **{k: r[k] for k in ["ratio", "low", "high"]}})
    forest(
        "Figure_S1_estimator_sensitivity",
        rows,
        ["results/python/estimates/national_models.csv", "results/python/estimates/conditional_blocks.csv"],
        groups=[7, 6, 5, 3, 2, 1],
        xlim=(0.05, 1.2),
        figheight=4.9,
    )
    traffic = read("traffic_support.csv")
    rows = []
    for th, label in [
        (1, "All matched: 44 jurisdictions / 3,892 sites"),
        (2, "≥2 or ≥10 sites: 42 jurisdictions / 3,890 sites"),
        (20, "≥20 sites: 40 jurisdictions / 3,865 sites"),
        (30, "≥30 sites: 38 jurisdictions / 3,821 sites"),
    ]:
        z = traffic[
            traffic.metric.eq("log_day_level") & traffic.band.eq("evening") & traffic.minimum_sites.eq(th)
        ].iloc[0]
        rows.append(
            {
                "label": label,
                "threshold": th,
                "ratio": z.ratio_equal_jurisdiction,
                "low": z.low,
                "high": z.high,
            }
        )
    forest(
        "Figure_S2_traffic_support",
        rows,
        ["results/python/estimates/traffic_support.csv"],
        xlim=(0.38, 1.06),
        figheight=3.9,
        xlabel="2020/2019 ratio of post/pre vehicle-volume changes (log scale)",
        kind="recomputed jurisdiction summaries",
    )
    # Flow: counts are stage-specific; the two fitted samples overlap.
    scope = support_table("model_support.csv")
    ped = scope[scope.endpoint.eq("fatal_pedestrian_crashes")].set_index("band")
    flow = support_table("source_flow.csv")
    assert int(flow.iloc[0]["n"]) == 843530 and int(flow.iloc[-1]["n"]) == 791667
    assert 791667 - 18708 == 772959
    fig, ax = plt.subplots(figsize=(9.0, 7.0))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")

    def box(x, y, w, h, text):
        ax.add_patch(Rectangle((x, y), w, h, fill=False, lw=1))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=10)

    def arrow(x, y, xx, yy):
        ax.annotate("", xy=(xx, yy), xytext=(x, y), arrowprops={"arrowstyle": "->", "lw": 1})

    box(23, 83, 54, 12, "National crash master\n843,530 crashes (2001–2024)")
    arrow(50, 83, 50, 75)
    box(23, 62, 54, 13, "Source-eligible records, including Arizona\n791,667 crashes")
    arrow(50, 62, 50, 54)
    box(23, 41, 54, 13, "Main study geography: 47 jurisdictions\n772,959 crashes")
    ax.text(80, 68, "51,863 records\nnot retained", fontsize=9, va="center")
    ax.text(80, 47, "18,708 Arizona\ncomparison records\nnot in the main cohort", fontsize=9, va="center")
    arrow(40, 41, 24, 31)
    arrow(60, 41, 76, 31)
    for x, band in [(1, "evening"), (52, "morning")]:
        r = ped.loc[band]
        txt = "\n".join(
            [
                f"{band.title()} pedestrian model",
                f"{int(r.unique_crashes_all_centers):,} unique crashes",
                f"{int(r.event_appearances):,} appearances",
                f"{int(r.informative_tables)} / {int(r.total_tables)} informative tables",
            ]
        )
        box(x, 11, 47, 20, txt)
    ax.text(
        50,
        4,
        "The two fitted samples share reference events; they are not additive.",
        ha="center",
        fontsize=10,
    )
    save(
        fig,
        "Figure_S3_sample_flow",
        scope,
        ["data/support/model_support.csv", "data/support/source_flow.csv"],
        "source-stage and model-support metadata; not a second source extraction",
    )
    c = national_counts()
    c = c[c.geo.eq("DST_main") & c.width_minutes.eq(40) & c.endpoint.eq("fatal_pedestrian_crashes") & c.sham.eq(0)]
    order = [("spring", "evening"), ("autumn", "evening"), ("spring", "morning"), ("autumn", "morning")]
    raw = []
    for season, band in order:
        d = c[c.season.eq(season) & c.band.eq(band)].groupby("period").events.sum()
        raw.append({"season": season, "band": band, "before": int(d["pre"]), "after": int(d["post"])})
    raw = pd.DataFrame(raw)
    fig, ax = plt.subplots(figsize=(9.2, 4.7))
    x = np.arange(4)
    w = 0.33
    for label, key, offset, h in [("Before", "before", -w / 2, None), ("After", "after", w / 2, "//")]:
        bars = ax.bar(x + offset, raw[key], width=w, label=label, hatch=h)
        ax.bar_label(bars, padding=4, fontsize=11)
    ax.set_ylim(0, 735)
    ax.set_xticks(
        x,
        [
            "Spring evening\nDaylight gain",
            "Autumn evening\nDaylight loss",
            "Spring morning\nDaylight loss",
            "Autumn morning\nDaylight gain",
        ],
    )
    ax.set_ylabel("Pedestrian-fatal-crash count")
    ax.legend(frameon=False, ncol=2)
    ax.spines[["top", "right"]].set_visible(False)
    fig.subplots_adjust(bottom=0.19)
    save(
        fig,
        "Figure_S4_seasonal_counts",
        raw,
        ["data/national.csv.gz"],
        "raw aggregated actual-center focal counts",
    )
    legacy = read("legacy_weekly_estimates.csv")
    assert len(legacy) == 28
    fig, ax = plt.subplots(figsize=(9.2, 5.6))
    fig.subplots_adjust(left=0.10, right=0.98, top=0.85, bottom=0.16)
    ax.set_yscale("log")
    ax.set_ylim(0.015, 32)
    ax.set_xlim(-4.6, 3.6)
    ax.axhline(1, ls=":", lw=1.1)
    ax.axvline(-0.5, ls="--", lw=1)
    ax.axvspan(-1.45, -0.55, fill=False, hatch="////", alpha=0.23)
    ax.text(-1, 8, "Reference\nweek −1", ha="center", fontsize=9)
    for j, (band, season, label, mark) in enumerate(series):
        d = legacy[legacy.band.eq(band) & legacy.season.eq(season)].sort_values("week")
        xx = d.week + (j - 1.5) * 0.095
        ax.errorbar(
            xx,
            d.ratio,
            yerr=np.array([d.ratio - d.low, d.high - d.ratio]),
            fmt=mark,
            ms=5,
            capsize=2.4,
            lw=1.1,
            label=label,
        )
    ax.set_xticks(range(-4, 4), ["−4", "−3", "−2", "−1", "0", "+1", "+2", "+3"])
    ax.set_yticks([0.03, 0.1, 0.3, 1, 3, 10, 30])
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, p: f"{v:g}"))
    ax.minorticks_off()
    ax.set_ylabel("Focal/reference ratio relative to week −1 (log scale)")
    ax.set_xlabel("Relative week at the actual clock change (Monday–Saturday)")
    ax.legend(frameon=False, ncol=2, loc="lower center", bbox_to_anchor=(0.5, 1.01))
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", alpha=0.15)
    save(
        fig,
        "Figure_S5_original_event_weeks",
        legacy,
        ["data/sensitivity/legacy_weeks.csv", "results/python/estimates/legacy_weekly_estimates.csv"],
        "28 re-estimated legacy diagnostic coefficients; reference not estimated",
    )
    cal = read("calendar_sensitivity.csv")
    names = [
        ("matched_14d_gap35", "14 days; centers ±35 days"),
        ("nonoverlap_14d_gap42", "14 days; centers ±42 days (no overlap)"),
        ("short_7d_gap35", "7 days; centers ±35 days"),
        ("long_21d_gap35", "21 days; centers ±35 days"),
    ]
    fig, ax = plt.subplots(figsize=(9.0, 4.0))
    fig.subplots_adjust(left=0.40, right=0.97, bottom=0.19, top=0.83)
    for j, (band, marker) in enumerate([("evening", "o"), ("morning", "s")]):
        z = cal[cal.band.eq(band)].set_index("specification").reindex([name for name, _ in names])
        ys = np.arange(4)[::-1] + (0.13 if j == 0 else -0.13)
        ax.errorbar(
            z.ratio,
            ys,
            xerr=np.array([z.ratio - z.low, z.high - z.ratio]),
            fmt=marker,
            ms=5,
            capsize=3,
            lw=1.2,
            label=band.title(),
        )
    ratio_axis(ax, 0.03, 1.5, [0.03, 0.1, 0.3, 1])
    ax.set_ylim(-0.55, 3.55)
    ax.set_yticks(np.arange(4)[::-1], [label for _, label in names])
    ax.tick_params(axis="y", length=0, pad=10)
    ax.set_xlabel("Gain-oriented ratio of count ratios (log scale)")
    ax.legend(frameon=False, loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=2)
    save(
        fig,
        "Figure_S6_calendar",
        cal,
        ["results/python/estimates/calendar_sensitivity.csv"],
        "8 re-estimated calendar sensitivities, same 2002–2023 cohort",
    )
    influence = read("leave_one_year.csv")
    fig, axes = plt.subplots(2, 1, figsize=(9.2, 6.6), sharex=True)
    fig.subplots_adjust(left=0.10, right=0.98, top=0.94, bottom=0.11, hspace=0.20)
    for ax, band, marker in zip(axes, ["evening", "morning"], ["o", "s"], strict=True):
        z = influence[influence.band.eq(band)].sort_values("omitted_policy_year")
        ax.axhline(z.full_ratio.iloc[0], ls="--", lw=1.1)
        ax.scatter(z.omitted_policy_year, z.ratio, marker=marker, s=25)
        ax.set_ylabel(f"{band.title()} RCR")
        ax.set_ylim(0.12 if band == "evening" else 0.16, 0.25 if band == "evening" else 0.26)
        ax.grid(axis="y", alpha=0.18)
        ax.spines[["top", "right"]].set_visible(False)
    axes[-1].set_xlabel("Omitted policy-transition year")
    axes[-1].set_xticks(range(2001, 2025, 2))
    save(
        fig,
        "Figure_S7_policy_year_influence",
        influence.drop(columns=["low", "high"]),
        ["results/python/estimates/leave_one_year.csv"],
        "leave-one-policy-year main-model estimates; dashed lines are full-model estimates",
    )
    (out / "FIGURE_DATA_MANIFEST.json").write_text(json.dumps(meta, indent=2) + "\n")
    assert len(meta) == 10
    return {
        "pass": True,
        "figures": 10,
        "legacy_fitted_rows": 28,
        "aligned_fitted_rows": 24,
        "use_of_expected_estimates_for_drawing": False,
        "source_metadata_figures": ["Figure_S3_sample_flow"],
        "plot_input_rows": {n: r["rows"] for n, r in meta.items()},
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", type=Path)
    a = ap.parse_args()
    print(json.dumps(run(a.outdir), indent=2))
