from __future__ import annotations

import numpy as np
from scipy import stats

from .conditional_poisson import evaluate, fit


def estimate(
    c,
    band="evening",
    reference="night",
    endpoint="fatal_pedestrian_crashes",
    geo="DST_main",
    width=40,
    years=None,
    seasons=None,
    corrected=True,
    independent=False,
):
    s = c[c.geo.eq(geo) & c.width_minutes.eq(width) & c.endpoint.eq(endpoint) & c.band.isin([band, reference])].copy()
    if years is not None:
        s = s[s.year.isin(years)]
    if seasons is not None:
        s = s[s.season.isin(seasons)]
    if not corrected:
        s = s[s.sham.eq(0)]
    tab = s.pivot(index=["year", "season", "sham"], columns=["period", "band"], values="events").sort_index()
    raw_n = len(tab)
    a = tab[("post", band)]
    b = tab[("post", reference)]
    cc = tab[("pre", band)]
    dd = tab[("pre", reference)]
    M = a + b + cc + dd
    F = a + cc
    P = a + b
    informative = np.maximum(0, P - (M - F)) < np.minimum(P, F)
    table_zero_information = tab.index[~informative].tolist()
    tab = tab[informative]
    meta = tab.index.to_frame(index=False)
    k = len(meta)
    if k < 6:
        raise ValueError("Insufficient informative strata")
    seasons_u = sorted(meta.season.unique())
    nuisance = []
    extra_names = []
    if corrected:
        for yy in sorted(meta.year.unique()):
            for ss in seasons_u:
                nuisance.append((meta.year.eq(yy) & meta.season.eq(ss)).to_numpy(float))
                extra_names.append(f"post_{yy}_{ss}")
        for ss in seasons_u:
            nuisance.append(meta.season.eq(ss).to_numpy(float) * meta.sham.to_numpy(float) / 35)
            extra_names.append("post_" + ss + "_centre_slope")
    direction = np.where(meta.season.eq("spring"), 1.0, -1.0) * (1 if band.startswith("evening") else -1)
    target = meta.sham.eq(0).to_numpy(float) * direction
    nuisance.append(target)
    extra_names.append("light_gain")
    E = np.column_stack(nuisance)
    # 2 periods x 2 categories per stratum; reference category's regressors are zero.
    X = []
    y = []
    groups = []
    yr = []
    gx = []
    succ = []
    fail = []
    for i, m in enumerate(meta.itertuples(index=False)):
        for post, period in enumerate(["pre", "post"]):
            vec = np.r_[np.eye(k)[i], post * E[i]]
            gx.append(vec)
            succ.append(int(tab.iloc[i][(period, band)]))
            fail.append(int(tab.iloc[i][(period, reference)]))
            for focal, name in [(0, reference), (1, band)]:
                X.append(vec * focal)
                y.append(int(tab.iloc[i][(period, name)]))
                groups.append(f"{i:04d}|{post}")
                yr.append(str(m.year))
    X = np.array(X)
    y = np.array(y)
    groups = np.array(groups)
    yr = np.array(yr)
    names = [f"stratum_{i}" for i in range(k)] + extra_names
    j = names.index("light_gain")
    fitres = fit(y, X, groups, np.repeat("national", len(y)), yr, maxiter=100)
    result = next(dict(v) for v in fitres["rows"] if v["method"] == "year" and v["term_index"] == j)
    use = fitres["used_rows"]
    Xu = X[use]
    yu = y[use]
    gu = groups[use]
    yu_year = yr[use]
    _, _, information, group_scores, starts = evaluate(fitres["beta"], Xu, yu, gu, np.zeros(len(use)))
    inverse_information = np.linalg.inv(information)
    inf = group_scores @ inverse_information[:, j]
    yg = yu_year[starts]
    unique = np.unique(yg)
    G = len(unique)
    ng = len(starts)
    p = X.shape[1]
    factor = G / (G - 1) * (ng - 1) / (ng - p)
    byyear = {str(v): float(inf[yg == v].sum() * np.sqrt(factor)) for v in unique}
    if not np.isclose(sum(v * v for v in byyear.values()), result["se"] ** 2, rtol=1e-8, atol=1e-12):
        raise AssertionError("Independent year sandwich mismatch")
    indep = None
    if independent:
        import statsmodels.api as sm

        fm = sm.GLM(np.column_stack([succ, fail]), np.array(gx), family=sm.families.Binomial()).fit(
            maxiter=200, tol=1e-10
        )
        diff = float(np.max(np.abs(fm.params - fitres["beta"])))
        indep = {
            "converged": bool(fm.converged),
            "max_all_parameter_difference": diff,
            "log_effect": float(fm.params[j]),
        }
        if not fm.converged or diff > 1e-6:
            raise AssertionError(("Independent binomial mismatch", indep))
    result.update(
        {
            "band": band,
            "reference": reference,
            "endpoint": endpoint,
            "geo": geo,
            "width": width,
            "corrected": corrected,
            "raw_strata": raw_n,
            "informative_strata": k,
            "zero_information_strata": table_zero_information,
            "years": sorted(int(z) for z in unique),
            "year_influence": byyear,
            "year_CR1_factor": factor,
            "audit": fitres["audit"],
            "coefficients": dict(zip(names, map(float, fitres["beta"]), strict=False)),
            "independent": indep,
            "input_events_across_reused_tables": int(y.sum()),
            "not_independent_event_count": True,
        }
    )
    return result


def contrast(a, b, name):
    ca = a["year_influence"]
    cb = b["year_influence"]
    shared = set(ca) & set(cb)
    cov = sum(ca[y] * cb[y] for y in shared)
    v = a["se"] ** 2 + b["se"] ** 2 - 2 * cov
    if v < -1e-10:
        raise AssertionError("negative direct variance")
    se = float(np.sqrt(max(v, 0)))
    d = a["log_ratio"] - b["log_ratio"]
    df = min(a["df"], b["df"])
    q = stats.t.ppf(0.975, df)
    return {
        "contrast": name,
        "log_ratio": d,
        "ratio": float(np.exp(d)),
        "low": float(np.exp(d - q * se)),
        "high": float(np.exp(d + q * se)),
        "p_value": float(2 * stats.t.sf(abs(d / se), df)),
        "se": se,
        "df": df,
        "shared_years": len(shared),
        "cross_covariance": cov,
    }
