"""Bounded post-result tests of classification periods, nuisance constraints,
common temporal weights, and calendar decomposition.

Uses supplied aggregate counts only. Primary model and data are unchanged.
The three-center conditional model conditions on each year-season's linear
center-score sufficient statistic; overlapping tables imply composite, not
independent joint, likelihood. All interval calculations use policy years.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import optimize, stats
from scipy.special import gammaln, logit, logsumexp

ENDPOINTS = ("fatal_pedestrian_crashes", "other_fatal_crashes")


def choose(n, k):
    return gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)


@dataclass
class Triplet:
    year: int
    season: str
    q: float
    observed: int
    support: np.ndarray
    log_weight: np.ndarray
    T: int
    U: int

    def moments(self, beta: float):
        lw = self.log_weight + self.q * beta * self.support
        norm = logsumexp(lw)
        prob = np.exp(lw - norm)
        mean = float(prob @ self.support)
        var = float(prob @ ((self.support - mean) ** 2))
        ix = np.flatnonzero(self.support == self.observed)
        if len(ix) != 1:
            raise ValueError("Observed count absent from conditional support")
        ll = float(lw[ix[0]] - norm)
        return ll, self.q * (self.observed - mean), self.q * self.q * var


def make_triplet(cells, year=2001, season="spring", q=1.0):
    """cells: three centers (-35,0,35) by a,b,c,d = post-focal,post-ref,pre-focal,pre-ref."""
    a = np.asarray(cells)
    if a.shape != (3, 4) or not np.isfinite(a).all() or (a < 0).any() or not np.equal(a, np.floor(a)).all():
        raise ValueError("Three nonnegative integer 2x2 tables required")
    a = a.astype(np.int64)
    P = a[:, 0] + a[:, 1]
    Q = a[:, 2] + a[:, 3]
    F = a[:, 0] + a[:, 2]
    lo = np.maximum(0, F - Q)
    hi = np.minimum(P, F)
    T = int(a[:, 0].sum())
    U = int(a[2, 0] - a[0, 0])
    values = []
    lw = []
    for z in range(int(lo[1]), int(hi[1]) + 1):
        if (T - z - U) % 2:
            continue
        minus = (T - z - U) // 2
        plus = (T - z + U) // 2
        xx = np.array([minus, z, plus], dtype=int)
        if np.any(xx < lo) or np.any(xx > hi):
            continue
        values.append(z)
        lw.append(float(np.sum(choose(P, xx) + choose(Q, F - xx))))
    if not values:
        raise ValueError("Empty conditioned triplet")
    return Triplet(
        int(year), str(season), float(q), int(a[1, 0]), np.asarray(values, dtype=float), np.asarray(lw), T, U
    )


def build_triplets(counts, band, endpoint, years=None):
    d = counts[
        counts.geo.eq("DST_main")
        & counts.width_minutes.eq(40)
        & counts.endpoint.eq(endpoint)
        & counts.band.isin([band, "night"])
    ]
    if years is not None:
        d = d[d.year.isin(years)]
    p = d.pivot(index=["year", "season", "sham"], columns=["period", "band"], values="events").sort_index()
    out = []
    for (y, s), g in p.groupby(level=["year", "season"]):
        g = g.droplevel(["year", "season"]).reindex([-35, 0, 35])
        if g.isna().any().any():
            raise ValueError("Incomplete triplet")
        cells = g[[("post", band), ("post", "night"), ("pre", band), ("pre", "night")]].to_numpy()
        q = (1 if s == "spring" else -1) * (1 if band == "evening" else -1)
        out.append(make_triplet(cells, y, s, q))
    return out


def fit_triplets(triplets):
    info_triplets = [t for t in triplets if len(t.support) > 1]
    if len(info_triplets) < 3:
        raise ValueError("Too few informative triplets")

    def vals(beta):
        return np.asarray([t.moments(float(beta)) for t in info_triplets]).sum(axis=0)

    left, right = vals(-40)[1], vals(40)[1]
    if not (left > 0 and right < 0):
        raise ValueError("Finite effect not supported (boundary likelihood)")
    beta = float(optimize.brentq(lambda b: vals(b)[1], -40, 40, xtol=1e-12))
    ll, score, info = vals(beta)
    if abs(score) > 1e-7 or info <= 0:
        raise ValueError("Score/information check failed")
    return {
        "log_ratio": beta,
        "ratio": float(np.exp(beta)),
        "score": float(score),
        "information": float(info),
        "informative_triplets": len(info_triplets),
        "total_triplets": len(triplets),
        "degenerate_triplets": len(triplets) - len(info_triplets),
        "log_likelihood": float(ll),
    }


def interval(beta, deletions):
    v = np.asarray(deletions, float)
    G = len(v)
    if G < 3 or not np.isfinite(v).all():
        raise ValueError("Incomplete deletion estimates")
    se = float(np.sqrt((G - 1) / G * np.sum((v - v.mean()) ** 2)))
    if se <= 0:
        raise ValueError("Degenerate deletion variance")
    q = stats.t.ppf(0.975, G - 1)
    return {
        "ratio": float(np.exp(beta)),
        "log_ratio": float(beta),
        "se": se,
        "low": float(np.exp(beta - q * se)),
        "high": float(np.exp(beta + q * se)),
        "p_value": float(2 * stats.t.sf(abs(beta / se), G - 1)),
        "df": G - 1,
        "policy_years": G,
    }


def standardized_inputs(counts, band):
    d = counts[
        counts.geo.eq("DST_main")
        & counts.width_minutes.eq(40)
        & counts.endpoint.isin(ENDPOINTS)
        & counts.band.isin([band, "night"])
    ]
    p = d.pivot(index=["year", "season", "endpoint", "sham", "period"], columns="band", values="events").reset_index()
    p["total"] = p[band] + p["night"]
    p["fraction"] = p[band] / p.total
    # Define exactly the same complete paired-season years for both endpoints, centers and sides.
    good = p.groupby("year").apply(lambda g: len(g) == 24 and g.total.gt(0).all(), include_groups=False)
    ys = good[good].index.tolist()
    excluded = good[~good].index.tolist()
    p = p[p.year.isin(ys)].copy()
    p["q"] = np.where(p.season.eq("spring"), 1.0, -1.0) * (1 if band == "evening" else -1)
    p["gain_side"] = (p.period.eq("post") & p.q.eq(1)) | (p.period.eq("pre") & p.q.eq(-1))
    return p, ys, excluded


def standardized_estimate(p, years):
    z = p[p.year.isin(years)]
    n = len(years) * 2
    avg = z.groupby(["endpoint", "sham", "gain_side"]).fraction.agg(["mean", "size"])
    if not avg["size"].eq(n).all():
        raise ValueError("Unequal common temporal support")
    D = {}
    rows = []
    for ep in ENDPOINTS:
        for k in [-35, 0, 35]:
            pg = float(avg.loc[(ep, k, True), "mean"])
            pl = float(avg.loc[(ep, k, False), "mean"])
            if not (0 < pg < 1 and 0 < pl < 1):
                raise ValueError("Boundary marginal probability")
            D[ep, k] = float(logit(pg) - logit(pl))
            rows.append(
                {
                    "endpoint": ep,
                    "sham": k,
                    "gain_fraction": pg,
                    "loss_fraction": pl,
                    "log_contrast": D[ep, k],
                    "contrast": float(np.exp(D[ep, k])),
                    "coefficient_in_symmetric_sum": 1.0 if k == 0 else -0.5,
                    "contribution": D[ep, k] * (1.0 if k == 0 else -0.5),
                    "weight_per_year_season": 1 / n,
                }
            )
    b = {e: D[e, 0] - 0.5 * (D[e, -35] + D[e, 35]) for e in ENDPOINTS}
    b["direct"] = b[ENDPOINTS[0]] - b[ENDPOINTS[1]]
    return b, rows
