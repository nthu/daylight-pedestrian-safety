"""Vectorized exact conditional Poisson, full-rank single-tolerance information path.

Grouping intercepts are conditioned out. No pseudoinverse or eigenvalue clipping.
Covariance sums are sensitivity estimators, not arbitrary-dependence guarantees.
"""

import numpy as np
from scipy.stats import norm, t


def _exp(value: float) -> float:
    """Exponentiate a log-scale value without hiding an infinite result."""

    with np.errstate(over="ignore"):
        return float(np.exp(value))


def evaluate(beta, X, y, g, offset):
    eta = X @ beta + offset
    unique, starts, n = np.unique(g, return_index=True, return_counts=True)
    tot = np.add.reduceat(y, starts)
    mx = np.maximum.reduceat(eta, starts)
    ee = np.exp(eta - np.repeat(mx, n))
    den = np.add.reduceat(ee, starts)
    p = ee / np.repeat(den, n)
    mu = p * np.repeat(tot, n)
    means = np.add.reduceat(X * p[:, None], starts, axis=0)
    ll = float(y @ eta - tot @ (mx + np.log(den)))
    res = y - mu
    score = X.T @ res
    H = X.T @ (X * mu[:, None]) - (means * tot[:, None]).T @ means
    gs = np.add.reduceat(X * res[:, None], starts, axis=0)
    return ll, score, (H + H.T) / 2, gs, starts


def fit(y, X, groups, county, year, offset=None, maxiter=100):
    y = np.asarray(y, float)
    X = np.asarray(X, float)
    groups = np.asarray(groups, str)
    county = np.asarray(county, str)
    year = np.asarray(year, str)
    o = np.zeros(len(y)) if offset is None else np.asarray(offset, float)
    if (
        X.ndim != 2
        or X.shape[0] != len(y)
        or o.shape != y.shape
        or not np.isfinite(X).all()
        or not np.isfinite(y).all()
        or not np.isfinite(o).all()
        or (y < 0).any()
        or not np.allclose(y, np.round(y))
    ):
        raise ValueError("Invalid inputs")
    ug, idx = np.unique(groups, return_inverse=True)
    tot = np.bincount(idx, weights=y)
    sizes = np.bincount(idx)
    keep = (tot[idx] > 0) & (sizes[idx] > 1)
    order = np.flatnonzero(keep)[np.argsort(groups[keep], kind="stable")]
    X = X[order]
    y = y[order]
    o = o[order]
    g = groups[order]
    county = county[order]
    year = year[order]
    if len(y) <= X.shape[1]:
        raise ValueError("No information")
    beta = np.zeros(X.shape[1])
    history = []
    converged = False
    for iteration in range(maxiter):
        ll, sc, H, gs, starts = evaluate(beta, X, y, g, o)
        eig = np.linalg.eigvalsh(H)
        if eig.min() <= eig.max() * 1e-12:
            raise ValueError("Information not full rank under one fixed tolerance")
        err = float(abs(sc).max())
        history.append({"iteration": iteration, "loglik_without_constants": ll, "score_max": err})
        if err < 1e-7:
            converged = True
            break
        step = np.linalg.solve(H, sc)
        scale = 1.0
        for _ in range(35):
            cand = beta + scale * step
            nl, ns, _, _, _ = evaluate(cand, X, y, g, o)
            if nl >= ll - 1e-9 and abs(ns).max() < err:
                beta = cand
                break
            scale *= 0.5
        else:
            raise ArithmeticError("No improving Newton step")
    if not converged:
        raise ArithmeticError("Iteration cap reached without convergence")
    ll, sc, H, gs, starts = evaluate(beta, X, y, g, o)
    bread = np.linalg.inv(H)
    inf = gs @ bread
    ng = len(starts)
    p = X.shape[1]
    cov = {"model": bread}
    nc = {}
    for label, z in [
        ("county", county[starts]),
        ("year", year[starts]),
        ("intersection", np.char.add(np.char.add(county[starts], "|"), year[starts])),
    ]:
        codes, ix = np.unique(z, return_inverse=True)
        n = len(codes)
        if n < 2:
            continue
        s = np.zeros((n, p))
        np.add.at(s, ix, inf)
        factor = n / (n - 1) * (ng - 1) / (ng - p)
        cov[label] = factor * (s.T @ s)
        nc[label] = n
    if "county" in cov and "year" in cov:
        cov["county_plus_year"] = cov["county"] + cov["year"]
        nc["county_plus_year"] = min(nc["county"], nc["year"])
        if "intersection" in cov:
            cov["two_way_subtraction"] = cov["county_plus_year"] - cov["intersection"]
            nc["two_way_subtraction"] = nc["county_plus_year"]
    summaries = []
    for method, V in cov.items():
        if method == "intersection":
            continue
        eig = np.linalg.eigvalsh((V + V.T) / 2)
        df = None if method == "model" else nc[method] - 1
        crit = norm.ppf(0.975) if df is None else t.ppf(0.975, df)
        for j, b in enumerate(beta):
            vv = float(V[j, j])
            se = np.sqrt(vv) if vv >= 0 else None
            rad = float(crit * se) if se is not None else None
            summaries.append(
                {
                    "term_index": j,
                    "method": method,
                    "log_ratio": float(b),
                    "ratio": _exp(b),
                    "se": se,
                    "low": _exp(b - rad) if rad is not None else None,
                    "high": _exp(b + rad) if rad is not None else None,
                    "p_value": float(2 * (norm.sf(abs(b / se)) if df is None else t.sf(abs(b / se), df)))
                    if se and se > 0
                    else None,
                    "df": df,
                    "full_covariance_psd": bool(eig.min() >= -1e-10 * max(float(abs(eig).max()), 1e-12)),
                    "minimum_eigenvalue": float(eig.min()),
                }
            )
    return {
        "beta": beta,
        "cov": cov,
        "rows": summaries,
        "used_rows": order,
        "information": H,
        "audit": {
            "iterations": iteration + 1,
            "converged": converged,
            "score_max": float(abs(sc).max()),
            "n_rows": len(y),
            "groups": ng,
            "events": int(y.sum()),
            "information_condition": float(np.linalg.cond(H)),
            "same_tolerance_rank_and_inverse": True,
            "history": history,
        },
    }
