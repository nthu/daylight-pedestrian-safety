"""Independent finite-support conditional composite likelihood for 2x2 counts.

No Q6/Q12 likelihood code is imported. Exact probability evaluation does not
imply exact joint sampling, confidence coverage, or randomized assignment.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import gammaln, logsumexp


def logchoose(n, k):
    n, k = np.asarray(n), np.asarray(k)
    return gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)


@dataclass
class ExactTables:
    cells: np.ndarray  # columns a=post focal,b=post ref,c=pre focal,d=pre ref
    X: np.ndarray
    years: np.ndarray
    names: list[str]

    def __post_init__(self):
        c = np.asarray(self.cells)
        if c.ndim != 2 or c.shape[1] != 4 or not np.isfinite(c).all():
            raise ValueError("Expected finite n-by-4 counts")
        if (c < 0).any() or not np.equal(c, np.floor(c)).all():
            raise ValueError("Counts must be nonnegative integers")
        self.cells = c.astype(np.int64)
        self.X = np.asarray(self.X, dtype=float)
        self.years = np.asarray(self.years)
        if self.X.shape[0] != len(c) or self.X.shape[1] != len(self.names):
            raise ValueError("Design dimensions do not match counts or names")
        if len(self.years) != len(c) or not np.isfinite(self.X).all():
            raise ValueError("Invalid design or cluster labels")
        a, b, cc, d = self.cells.T
        self.a = a
        self.P, self.Q, self.F = a + b, cc + d, a + cc
        self.lo = np.maximum(0, self.F - self.Q)
        self.hi = np.minimum(self.P, self.F)
        self.informative = self.hi > self.lo
        self.support = self.lo[:, None] + np.arange(int((self.hi - self.lo).max()) + 1)[None, :]
        valid = self.support <= self.hi[:, None]
        safe = np.minimum(self.support, self.hi[:, None])
        lw = logchoose(self.P[:, None], safe) + logchoose(self.Q[:, None], self.F[:, None] - safe)
        self.logw = np.where(valid, lw, -np.inf)
        self.logw_obs = logchoose(self.P, a) + logchoose(self.Q, self.F - a)

    def evaluate(self, beta):
        beta = np.asarray(beta, dtype=float)
        eta = self.X @ beta
        logz = self.logw + self.support * eta[:, None]
        norm = logsumexp(logz, axis=1)
        prob = np.exp(logz - norm[:, None])
        mean = np.sum(prob * self.support, axis=1)
        var = np.sum(prob * (self.support - mean[:, None]) ** 2, axis=1)
        residual = self.a - mean
        scores = self.X * residual[:, None]
        ll = float(np.sum(self.logw_obs + self.a * eta - norm))
        return ll, scores.sum(axis=0), self.X.T @ (self.X * var[:, None]), scores, mean, var

    def fit(self, start=None, tol=2e-8, maxiter=100):
        beta = np.zeros(self.X.shape[1]) if start is None else np.asarray(start, float).copy()
        history = []
        for it in range(maxiter):
            ll, score, info, *_ = self.evaluate(beta)
            history.append({"iteration": it, "loglik": ll, "score_max": float(abs(score).max())})
            if abs(score).max() < tol:
                break
            eig = np.linalg.eigvalsh(info)
            if eig[0] <= eig[-1] * 1e-12:
                raise RuntimeError("Information rank failure; no eigenvalue clipping")
            step = np.linalg.solve(info, score)
            t = 1.0
            for _ in range(35):
                trial = beta + t * step
                ll2 = self.evaluate(trial)[0]
                if ll2 >= ll - 1e-10:
                    beta = trial
                    break
                t *= 0.5
            else:
                raise RuntimeError("Newton line search did not improve likelihood")
        else:
            raise RuntimeError("Newton iteration limit")
        ll, score, info, scores, mean, var = self.evaluate(beta)
        eig = np.linalg.eigvalsh(info)
        if abs(score).max() >= tol or eig[0] <= eig[-1] * 1e-12:
            raise RuntimeError("Final score or information validation failed")
        inv = np.linalg.inv(info)
        unique = sorted(set(self.years.tolist()))
        G = len(unique)
        if G < 2:
            raise ValueError("At least two clusters required for reported uncertainty")
        U = np.stack([scores[self.years == yy].sum(axis=0) for yy in unique])
        cov = inv @ (U.T @ U) @ inv * G / (G - 1)
        j = self.names.index("light_gain")

        def interval(se, df):
            q = stats.t.ppf(0.975, df)
            return [float(np.exp(beta[j] - q * se)), float(np.exp(beta[j] + q * se))]

        se = float(np.sqrt(cov[j, j]))
        result = {
            "log_ratio": float(beta[j]),
            "ratio": float(np.exp(beta[j])),
            "year_sandwich_se_G_only": se,
            "year_sandwich_CI_G_only": interval(se, G - 1),
            "working_information_se": float(np.sqrt(inv[j, j])),
            "working_information_CI_t": interval(float(np.sqrt(inv[j, j])), G - 1),
            "years": list(map(int, unique)),
            "year_clusters": G,
            "tables_total": len(self.cells),
            "tables_informative": int(self.informative.sum()),
            "tables_degenerate": int((~self.informative).sum()),
            "n_parameters": len(beta),
            "score_max": float(abs(score).max()),
            "information_min_eigenvalue": float(eig[0]),
            "information_condition": float(eig[-1] / eig[0]),
            "loglik_composite": ll,
            "iterations": len(history),
            "coefficients": dict(zip(self.names, map(float, beta), strict=False)),
            "covariance_jackknife_not_yet_computed": True,
        }
        return result, beta


def make_model(
    frame: pd.DataFrame, band: str, years: Sequence[int] | None = None, seasons: Sequence[str] | None = None
):
    f = frame[frame.band.eq(band)].copy()
    if years is not None:
        f = f[f.year.isin(years)]
    if seasons is not None:
        f = f[f.season.isin(seasons)]
    if f.duplicated(["year", "season", "sham"]).any():
        raise ValueError("Duplicated table identity")
    f = f.sort_values(["year", "season", "sham"]).reset_index(drop=True)
    if set(f.sham.unique()) != {-35, 0, 35}:
        raise ValueError("Real and both sham centres are required")
    groups = f.groupby(["year", "season"]).sham.apply(lambda v: sorted(v.tolist()))
    if not all(v == [-35, 0, 35] for v in groups):
        raise ValueError("Incomplete real/sham triplet")
    ys = sorted(set(zip(f.year, f.season, strict=False)))
    names = [f"yearseason_{y}_{s}" for y, s in ys]
    columns = [(f.year.eq(y) & f.season.eq(s)).to_numpy(float) for y, s in ys]
    for ss in sorted(f.season.unique()):
        names.append("slope_" + ss)
        columns.append(f.season.eq(ss).to_numpy(float) * f.sham.to_numpy(float) / 35)
    direction = np.where(f.season.eq("spring"), 1.0, -1.0) * (1.0 if band.startswith("evening") else -1.0)
    names.append("light_gain")
    columns.append(f.sham.eq(0).to_numpy(float) * direction)
    c = f[["focal_post", "dark_post", "focal_pre", "dark_pre"]].to_numpy()
    return ExactTables(c, np.column_stack(columns), f.year.to_numpy(), names), f


def jackknife_interval(full_beta: float, deleted_betas: Sequence[float]):
    b = np.asarray(deleted_betas, dtype=float)
    K = len(b)
    if K < 3 or not np.isfinite(b).all():
        raise ValueError("Need at least three complete deletion fits")
    se = float(np.sqrt((K - 1) / K * np.sum((b - b.mean()) ** 2)))
    q = stats.t.ppf(0.975, K - 1)
    return {
        "blocks": K,
        "df": K - 1,
        "log_ratio": full_beta,
        "ratio": float(np.exp(full_beta)),
        "se": se,
        "low": float(np.exp(full_beta - q * se)),
        "high": float(np.exp(full_beta + q * se)),
        "p_value": float(2 * stats.t.sf(abs(full_beta / se), K - 1)),
        "deleted_ratio_min": float(np.exp(b.min())),
        "deleted_ratio_max": float(np.exp(b.max())),
    }


def scipy_probability_check(model: ExactTables, beta):
    ll, score, info, scores, means, variances = model.evaluate(beta)
    errors = []
    for i in range(len(model.cells)):
        eta = float(model.X[i] @ beta)
        x = np.arange(model.lo[i], model.hi[i] + 1)
        if len(x) == 1:
            errors.append(0.0)
            continue
        p = stats.nchypergeom_fisher.pmf(x, int(model.P[i] + model.Q[i]), int(model.P[i]), int(model.F[i]), np.exp(eta))
        lcustom = model.logw[i, : len(x)] + eta * x
        q = np.exp(lcustom - logsumexp(lcustom))
        errors.append(float(np.max(abs(p - q))))
    return float(max(errors))
