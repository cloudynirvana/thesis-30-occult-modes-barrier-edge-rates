#!/usr/bin/env python3
"""Joint hybrid-mode and barrier-conductance toy for Thesis #30.

One three-node anatomical graph. Two directed edges. Each edge carries an
active-mode shedding rate, a desmoplastic conductance, and a quiescent
multiplier. The field switches once, at a known time, from the active mode
to the quiescent mode.

Site burdens see only the series flux. A mode-aware barrier schedule records
the stalled fraction once in each mode. A mode-blind schedule records one
duration-weighted pool per edge and no mode tag. A mode-partial schedule
resolves one edge and pools the other.

Fisher ranks and profiles are deterministic. Seed 20260921 draws the single
noisy cloud. Research only. Not a medical device.

The graph is not the four-node laboratory of Thesis #21, and the modes are
not the three-coordinate observer of Thesis #20. No number is read from
either deposit's results file.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parent
FIG = ROOT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

SEED = 20260921
DT = 0.05
TAU = 18.0
T_END = 48.0
SAMPLE_TIMES = np.arange(0.0, T_END + 1e-9, 2.0)
# Last active sample is the node at t = TAU. Quiescence starts after that node.

R_SOIL = np.array([0.12, 0.06, 0.045])
K_SOIL = np.array([6.0, 4.5, 3.5])
X0 = np.array([1.0, 0.0, 0.0])

# Active shedding, conductance, quiescent multiplier. Edge 1 is primary→filter.
# Edge 2 is filter→liver and is transport-limited in the active mode.
TRUTH = np.array([0.08, 0.24, 0.35, 0.05, 0.022, 0.50])
NAMES = ["lam1", "kap1", "alpha1", "lam2", "kap2", "alpha2"]
EDGE_OF = {
    "lam1": 1,
    "kap1": 1,
    "alpha1": 1,
    "lam2": 2,
    "kap2": 2,
    "alpha2": 2,
}
THETA0 = np.log(TRUTH)

SIGMA_NODE = 0.05
SIGMA_STALL = 0.02
FD_EPS = 1e-4
STRUCT_FLOOR = 1e-6
PRACTICAL_REL = 0.5
PROFILE_CHI = 3.841  # chi-square, 1 df, 95 percent
PROFILE_MULTS = np.array([0.45, 0.65, 0.80, 1.00, 1.25, 1.60, 2.20])
# Finer grid on the transport-limited shedding rate, where the coarse cut is ambiguous.
PROFILE_FINE = np.array([0.70, 0.78, 0.86, 0.92, 1.00, 1.08, 1.16, 1.28, 1.40, 1.55])

# Absolute bounds used by profiles and the noisy cloud.
LOG_LO = np.log(np.array([1e-4, 1e-4, 0.02, 1e-4, 1e-4, 0.02]))
LOG_HI = np.log(np.array([2.0, 2.0, 3.0, 2.0, 2.0, 3.0]))

W_ACTIVE = TAU / T_END
W_QUIET = (T_END - TAU) / T_END

SCHEDULES = (
    "nodes",
    "mode_blind",
    "mode_partial_proximal",
    "mode_partial_distal",
    "mode_aware",
)


def rk4_step(x: np.ndarray, h: float, phi: np.ndarray) -> np.ndarray:
    def f(z: np.ndarray) -> np.ndarray:
        p, filt, liver = z
        dp = R_SOIL[0] * p * (1.0 - p / K_SOIL[0]) - phi[0] * p
        df = (
            R_SOIL[1] * filt * (1.0 - filt / K_SOIL[1])
            - phi[1] * filt
            + phi[0] * p
        )
        dl = R_SOIL[2] * liver * (1.0 - liver / K_SOIL[2]) + phi[1] * filt
        return np.array([dp, df, dl])

    k1 = f(x)
    k2 = f(x + 0.5 * h * k1)
    k3 = f(x + 0.5 * h * k2)
    k4 = f(x + h * k3)
    return x + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def advance(x: np.ndarray, t: float, t_target: float, phi: np.ndarray) -> np.ndarray:
    while t_target - t > 1e-12:
        h = min(DT, t_target - t)
        x = rk4_step(x, h, phi)
        t = t + h
    return x


def trajectory(phi_a: np.ndarray, phi_q: np.ndarray) -> np.ndarray:
    """Node path. phi_* are length-2 series fluxes, edge 1 then edge 2."""
    x = X0.copy()
    out = np.empty((SAMPLE_TIMES.size, 3))
    t = 0.0
    si = 0
    segments = ((0.0, TAU, phi_a), (TAU, T_END, phi_q))
    for t0, t1, phi in segments:
        if abs(t - t0) > 1e-8:
            x = advance(x, t, t0, phi)
            t = t0
        while si < SAMPLE_TIMES.size and SAMPLE_TIMES[si] <= t1 + 1e-10:
            target = float(SAMPLE_TIMES[si])
            x = advance(x, t, target, phi)
            t = target
            out[si] = x
            si += 1
        if t1 - t > 1e-12:
            x = advance(x, t, t1, phi)
            t = t1
    if si != SAMPLE_TIMES.size:
        raise RuntimeError(f"sampler stopped at {si} of {SAMPLE_TIMES.size}")
    return out


def edge_laws(lam: float, kap: float, alpha: float) -> dict:
    la = lam
    lq = alpha * lam
    return {
        "phi_a": la * kap / (la + kap),
        "phi_q": lq * kap / (lq + kap),
        "s_a": la / (la + kap),
        "s_q": lq / (lq + kap),
        "lam_q": lq,
    }


def unpack(theta: np.ndarray) -> tuple[float, ...]:
    v = np.exp(theta)
    return tuple(float(z) for z in v)


def laws_from_theta(theta: np.ndarray) -> tuple[dict, dict]:
    lam1, kap1, a1, lam2, kap2, a2 = unpack(theta)
    return edge_laws(lam1, kap1, a1), edge_laws(lam2, kap2, a2)


def pool(law: dict) -> float:
    return W_ACTIVE * law["s_a"] + W_QUIET * law["s_q"]


def predict(theta: np.ndarray, schedule: str) -> np.ndarray:
    law1, law2 = laws_from_theta(theta)
    phi_a = np.array([law1["phi_a"], law2["phi_a"]])
    phi_q = np.array([law1["phi_q"], law2["phi_q"]])
    nodes = trajectory(phi_a, phi_q).ravel()
    if schedule == "nodes":
        extra = np.zeros(0)
    elif schedule == "mode_blind":
        extra = np.array([pool(law1), pool(law2)])
    elif schedule == "mode_partial_proximal":
        extra = np.array([law1["s_a"], law1["s_q"], pool(law2)])
    elif schedule == "mode_partial_distal":
        extra = np.array([pool(law1), law2["s_a"], law2["s_q"]])
    elif schedule == "mode_aware":
        extra = np.array([law1["s_a"], law1["s_q"], law2["s_a"], law2["s_q"]])
    else:
        raise KeyError(schedule)
    return np.concatenate([nodes, extra])


def sigma_of(schedule: str, n: int) -> np.ndarray:
    n_node = SAMPLE_TIMES.size * 3
    sig = np.empty(n)
    sig[:n_node] = SIGMA_NODE
    sig[n_node:] = SIGMA_STALL
    if sig.size != n:
        raise RuntimeError("sigma length")
    return sig


def sensitivity(theta: np.ndarray, schedule: str) -> tuple[np.ndarray, np.ndarray]:
    base_n = predict(theta, schedule).size
    sig = sigma_of(schedule, base_n)
    s = np.empty((base_n, theta.size))
    for j in range(theta.size):
        tp = theta.copy()
        tm = theta.copy()
        tp[j] += FD_EPS
        tm[j] -= FD_EPS
        s[:, j] = (predict(tp, schedule) - predict(tm, schedule)) / (2.0 * FD_EPS)
    return s, sig


def fisher_report(theta: np.ndarray, schedule: str) -> dict:
    sens, sig = sensitivity(theta, schedule)
    z = sens / sig[:, None]
    _u, sv, vt = np.linalg.svd(z, full_matrices=False)
    sv_ratio = sv / sv[0]
    structural = int(np.sum(sv_ratio > STRUCT_FLOOR))
    with np.errstate(divide="ignore"):
        principal_rel = np.where(sv > 0.0, 1.0 / sv, np.inf)
    practical = int(np.sum(principal_rel < PRACTICAL_REL))
    keep = sv_ratio > 1e-8
    inv_eig = np.zeros_like(sv)
    inv_eig[keep] = 1.0 / np.square(sv[keep])
    v = vt.T
    cov = (v * inv_eig) @ v.T
    if np.any(~keep):
        mass = np.sum(np.square(v[:, ~keep]), axis=1)
    else:
        mass = np.zeros(theta.size)
    marginal = np.sqrt(np.clip(np.diag(cov), 0.0, None))
    marginal = np.where(mass > 1e-6, np.inf, marginal)
    # Quiescent shedding is log λ + log α. Infinite if either factor is.
    lq = []
    for i_lam, i_alpha, edge in ((0, 2, 1), (3, 5, 2)):
        if not np.isfinite(marginal[i_lam]) or not np.isfinite(marginal[i_alpha]):
            se = np.inf
        else:
            var = cov[i_lam, i_lam] + cov[i_alpha, i_alpha] + 2.0 * cov[i_lam, i_alpha]
            se = float(np.sqrt(max(var, 0.0)))
        lq.append({"edge": edge, "name": f"lam{edge}_quiescent", "marginal_rel_se": se})
    stays = []
    for i, name in enumerate(NAMES):
        if not name.startswith("lam"):
            continue
        se = float(marginal[i])
        stays.append(
            {
                "name": name,
                "edge": EDGE_OF[name],
                "marginal_rel_se": se,
                "stays": bool(np.isfinite(se) and se < PRACTICAL_REL),
            }
        )
    for row in lq:
        row["stays"] = bool(np.isfinite(row["marginal_rel_se"]) and row["marginal_rel_se"] < PRACTICAL_REL)
        stays.append(row)
    return {
        "schedule": schedule,
        "n_obs": int(sens.shape[0]),
        "n_param": int(theta.size),
        "singular_values": sv.tolist(),
        "singular_value_ratios": sv_ratio.tolist(),
        "structural_rank": structural,
        "practical_rank": practical,
        "principal_rel_se": principal_rel.tolist(),
        "condition_structural": float(sv[0] / sv[structural - 1]) if structural else None,
        "marginal_rel_se": {NAMES[i]: float(marginal[i]) for i in range(theta.size)},
        "edge_rates": stays,
    }


def algebraic_jacobian(theta: np.ndarray) -> dict:
    """Chain-rule map from log parameters to fluxes and stalls. No ODE."""

    def pack(th: np.ndarray, kind: str) -> np.ndarray:
        law1, law2 = laws_from_theta(th)
        fluxes = np.array([law1["phi_a"], law1["phi_q"], law2["phi_a"], law2["phi_q"]])
        if kind == "nodes":
            extra = np.zeros(0)
        elif kind == "mode_blind":
            extra = np.array([pool(law1), pool(law2)])
        elif kind == "mode_partial_proximal":
            extra = np.array([law1["s_a"], law1["s_q"], pool(law2)])
        elif kind == "mode_partial_distal":
            extra = np.array([pool(law1), law2["s_a"], law2["s_q"]])
        elif kind == "mode_aware":
            extra = np.array(
                [law1["s_a"], law1["s_q"], law2["s_a"], law2["s_q"]]
            )
        else:
            raise KeyError(kind)
        return np.concatenate([fluxes, extra])

    out = {}
    for kind in SCHEDULES:
        base = pack(theta, kind)
        jac = np.empty((base.size, theta.size))
        for j in range(theta.size):
            tp = theta.copy()
            tm = theta.copy()
            tp[j] += FD_EPS
            tm[j] -= FD_EPS
            jac[:, j] = (pack(tp, kind) - pack(tm, kind)) / (2.0 * FD_EPS)
        sv = np.linalg.svd(jac, compute_uv=False)
        ratio = sv / sv[0]
        out[kind] = {
            "singular_values": sv.tolist(),
            "structural_rank": int(np.sum(ratio > STRUCT_FLOOR)),
            "outputs": int(base.size),
        }
    return out


def hyperbola_twin(theta: np.ndarray) -> dict:
    """A second (λ, κ, α) with the same two fluxes and a different stall."""
    law1, law2 = laws_from_theta(theta)
    twins = []
    new_s = (0.42, 0.55)
    params = []
    for law, s_new in zip((law1, law2), new_s):
        phi_a = law["phi_a"]
        phi_q = law["phi_q"]
        lam = phi_a / (1.0 - s_new)
        kap = phi_a / s_new
        alpha = phi_q * kap / (lam * (kap - phi_q))
        check = edge_laws(lam, kap, alpha)
        params.extend([lam, kap, alpha])
        twins.append(
            {
                "s_a_new": s_new,
                "lam": lam,
                "kap": kap,
                "alpha": alpha,
                "phi_a_err": abs(check["phi_a"] - phi_a),
                "phi_q_err": abs(check["phi_q"] - phi_q),
                "s_a": check["s_a"],
                "s_q": check["s_q"],
            }
        )
    theta_twin = np.log(np.array(params))
    y0 = predict(theta, "nodes")
    y1 = predict(theta_twin, "nodes")
    s0 = predict(theta, "mode_aware")[-4:]
    s1 = predict(theta_twin, "mode_aware")[-4:]
    return {
        "edges": twins,
        "node_rmse": float(np.sqrt(np.mean(np.square(y0 - y1)))),
        "stall_rmse": float(np.sqrt(np.mean(np.square(s0 - s1)))),
        "max_abs_log_shift": float(np.max(np.abs(theta_twin - theta))),
    }


def _free_index(fixed: int) -> np.ndarray:
    return np.array([i for i in range(6) if i != fixed])


def profile_parameter(schedule: str, fixed: int, multipliers: np.ndarray) -> dict:
    """Noise-free profile. Other log-parameters re-fit. Continuation from the truth."""
    truth_y = predict(THETA0, schedule)
    sig = sigma_of(schedule, truth_y.size)
    free = _free_index(fixed)
    lo = LOG_LO[free]
    hi = LOG_HI[free]

    def residuals(free_theta: np.ndarray, fixed_value: float) -> np.ndarray:
        th = np.empty(6)
        th[fixed] = fixed_value
        th[free] = free_theta
        return (predict(th, schedule) - truth_y) / sig

    # Walk outward from the truth so each start is the previous optimum.
    center = THETA0[fixed]
    rows = []
    # Evaluate the truth first, then each multiplier using the nearest cached start.
    sequence = [1.0] + [float(m) for m in multipliers if abs(float(m) - 1.0) > 1e-12]
    # Re-order sequence by walking left and right.
    left = sorted([m for m in sequence if m < 1.0], reverse=True)
    right = sorted([m for m in sequence if m > 1.0])
    sequence = [1.0] + left + right
    # Restart the right walk from the truth, not from the far left point.
    starts = {1.0: THETA0[free].copy()}
    for label, group in (("left", left), ("right", right)):
        prev = THETA0[free].copy()
        prev_m = 1.0
        for m in group:
            fixed_value = center + np.log(m)
            res = least_squares(
                residuals,
                prev,
                args=(fixed_value,),
                bounds=(lo, hi),
                method="trf",
                xtol=1e-10,
                ftol=1e-10,
                gtol=1e-10,
                max_nfev=48,
            )
            cost = float(np.dot(res.fun, res.fun))
            th = np.empty(6)
            th[fixed] = fixed_value
            th[free] = res.x
            rows.append(
                {
                    "multiplier": m,
                    "value": float(np.exp(fixed_value)),
                    "weighted_rss": cost,
                    "inside_threshold": bool(cost <= PROFILE_CHI),
                    "success": bool(res.success),
                    "nfev": int(res.nfev),
                    "max_abs_log_shift_free": float(np.max(np.abs(res.x - THETA0[free]))),
                    "walk": label,
                    "seeded_from_multiplier": prev_m,
                }
            )
            prev = res.x.copy()
            prev_m = m
            starts[m] = res.x.copy()
    # Truth row.
    fun0 = residuals(THETA0[free], center)
    rows.append(
        {
            "multiplier": 1.0,
            "value": float(np.exp(center)),
            "weighted_rss": float(np.dot(fun0, fun0)),
            "inside_threshold": True,
            "success": True,
            "nfev": 0,
            "max_abs_log_shift_free": 0.0,
            "walk": "truth",
            "seeded_from_multiplier": 1.0,
        }
    )
    rows.sort(key=lambda r: r["multiplier"])
    inside = [r["multiplier"] for r in rows if r["inside_threshold"]]
    return {
        "schedule": schedule,
        "parameter": NAMES[fixed],
        "threshold": PROFILE_CHI,
        "grid": rows,
        "multipliers_inside": inside,
        "n_inside": len(inside),
        "n_grid": len(rows),
    }


def continuous_alias() -> dict:
    """One (λ, κ) per edge, fitted to the hybrid mode-blind record."""
    target = predict(THETA0, "mode_blind")
    sig = sigma_of("mode_blind", target.size)

    def predict_c(psi: np.ndarray) -> np.ndarray:
        lam1, kap1, lam2, kap2 = np.exp(psi)
        # α = 1 forces the two modes to share a flux. The stall is the active formula.
        law1 = edge_laws(lam1, kap1, 1.0)
        law2 = edge_laws(lam2, kap2, 1.0)
        phi = np.array([law1["phi_a"], law2["phi_a"]])
        nodes = trajectory(phi, phi).ravel()
        stalls = np.array([law1["s_a"], law2["s_a"]])
        return np.concatenate([nodes, stalls])

    def residuals(psi: np.ndarray) -> np.ndarray:
        return (predict_c(psi) - target) / sig

    lo = np.log(np.array([1e-4, 1e-4, 1e-4, 1e-4]))
    hi = np.log(np.array([2.0, 2.0, 2.0, 2.0]))
    starts = [
        np.log(np.array([0.08, 0.24, 0.05, 0.022])),
        np.log(np.array([0.08 * 0.35, 0.24, 0.05 * 0.50, 0.022])),
        np.log(np.array([0.04, 0.10, 0.03, 0.04])),
        np.log(np.array([0.12, 0.08, 0.08, 0.015])),
    ]
    fits = []
    for psi0 in starts:
        res = least_squares(
            residuals,
            psi0,
            bounds=(lo, hi),
            method="trf",
            xtol=1e-12,
            ftol=1e-12,
            gtol=1e-12,
            max_nfev=80,
        )
        fits.append(
            {
                "success": bool(res.success),
                "nfev": int(res.nfev),
                "weighted_rss": float(np.dot(res.fun, res.fun)),
                "psi": res.x.tolist(),
                "lam1": float(np.exp(res.x[0])),
                "kap1": float(np.exp(res.x[1])),
                "lam2": float(np.exp(res.x[2])),
                "kap2": float(np.exp(res.x[3])),
            }
        )
    best = min(fits, key=lambda r: r["weighted_rss"])
    psi = np.array(best["psi"])
    # Local Fisher of the continuous model at that fit, on its own mean.
    sens = np.empty((target.size, 4))
    for j in range(4):
        tp = psi.copy()
        tm = psi.copy()
        tp[j] += FD_EPS
        tm[j] -= FD_EPS
        sens[:, j] = (predict_c(tp) - predict_c(tm)) / (2.0 * FD_EPS)
    z = sens / sig[:, None]
    sv = np.linalg.svd(z, compute_uv=False)
    ratio = sv / sv[0]
    principal = 1.0 / sv
    law_a_1 = edge_laws(TRUTH[0], TRUTH[1], TRUTH[2])
    law_a_2 = edge_laws(TRUTH[3], TRUTH[4], TRUTH[5])
    return {
        "starts": fits,
        "best": best,
        "structural_rank": int(np.sum(ratio > STRUCT_FLOOR)),
        "practical_rank": int(np.sum(principal < PRACTICAL_REL)),
        "singular_values": sv.tolist(),
        "principal_rel_se": principal.tolist(),
        "active_lam": [float(TRUTH[0]), float(TRUTH[3])],
        "quiescent_lam": [float(law_a_1["lam_q"]), float(law_a_2["lam_q"])],
        "log_abs_error_vs_active": [
            abs(np.log(best["lam1"] / TRUTH[0])),
            abs(np.log(best["lam2"] / TRUTH[3])),
        ],
        "log_abs_error_vs_quiescent": [
            abs(np.log(best["lam1"] / law_a_1["lam_q"])),
            abs(np.log(best["lam2"] / law_a_2["lam_q"])),
        ],
    }


def noisy_cloud() -> dict:
    rng = np.random.default_rng(SEED)
    out = {}
    for schedule in ("mode_blind", "mode_aware"):
        truth = predict(THETA0, schedule)
        sig = sigma_of(schedule, truth.size)
        draw = truth + rng.normal(0.0, sig)
        free_lo = LOG_LO
        free_hi = LOG_HI

        def residuals(th: np.ndarray, draw=draw, sig=sig, schedule=schedule) -> np.ndarray:
            return (predict(th, schedule) - draw) / sig

        starts = [THETA0.copy()]
        for _ in range(5):
            starts.append(np.clip(THETA0 + rng.normal(0.0, 0.35, size=6), free_lo, free_hi))
        fits = []
        for psi0 in starts:
            res = least_squares(
                residuals,
                psi0,
                bounds=(free_lo, free_hi),
                method="trf",
                xtol=1e-10,
                ftol=1e-10,
                gtol=1e-10,
                max_nfev=40,
            )
            cost = float(np.dot(res.fun, res.fun))
            est = np.exp(res.x)
            fits.append(
                {
                    "success": bool(res.success),
                    "nfev": int(res.nfev),
                    "weighted_rss": cost,
                    "estimate": est.tolist(),
                    "max_abs_log_error": float(np.max(np.abs(res.x - THETA0))),
                }
            )
        costs = np.array([f["weighted_rss"] for f in fits])
        best = float(costs.min())
        # Retain fits within 3.841 of the best, the same one-degree cut.
        retain = [f for f in fits if f["weighted_rss"] <= best + PROFILE_CHI]
        est = np.array([f["estimate"] for f in retain])
        out[schedule] = {
            "n_starts": len(fits),
            "n_retained": len(retain),
            "best_weighted_rss": best,
            "noise_sum_squares": float(np.sum(((draw - truth) / sig) ** 2)),
            "n_obs": int(truth.size),
            "fits": fits,
            "retained_min": est.min(axis=0).tolist() if len(retain) else None,
            "retained_max": est.max(axis=0).tolist() if len(retain) else None,
            "names": NAMES,
        }
    return out


def sanity(truth_laws: tuple[dict, dict], twin: dict, reports: dict) -> None:
    for law in truth_laws:
        for key in ("phi_a", "phi_q", "s_a", "s_q"):
            if not np.isfinite(law[key]) or law[key] <= 0.0:
                raise AssertionError(key)
        if not (law["phi_a"] < min(TRUTH[0], TRUTH[1]) + 1.0):
            raise AssertionError("phi bound")
        if not (0.0 < law["s_a"] < 1.0 and 0.0 < law["s_q"] < 1.0):
            raise AssertionError("stall")
    if twin["node_rmse"] > 1e-8:
        raise AssertionError(f"twin leaked into the nodes: {twin['node_rmse']}")
    if twin["stall_rmse"] < 1e-3:
        raise AssertionError("twin did not move the stalls")
    # Nodes cannot see six parameters. Structural rank at most 4.
    if reports["nodes"]["structural_rank"] > 4:
        raise AssertionError("node rank broke the flux ceiling")
    y = predict(THETA0, "mode_aware")
    if not np.all(np.isfinite(y)):
        raise AssertionError("non-finite prediction")


def _finite(x: float) -> float | None:
    if x is None or not np.isfinite(x):
        return None
    return float(x)


def json_ready(obj):
    if isinstance(obj, dict):
        return {k: json_ready(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [json_ready(v) for v in obj]
    if isinstance(obj, float):
        if not np.isfinite(obj):
            return None
        return obj
    if isinstance(obj, (np.floating,)):
        v = float(obj)
        return None if not np.isfinite(v) else v
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, np.ndarray):
        return json_ready(obj.tolist())
    return obj


def make_figures(path: np.ndarray, reports: dict, profiles: list[dict]) -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Serif",
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.dpi": 140,
        }
    )
    # Figure 4-1. Node paths.
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    labels = ("primary", "filter", "liver")
    styles = ("-", "--", ":")
    for j, (lab, style) in enumerate(zip(labels, styles)):
        ax.plot(SAMPLE_TIMES, path[:, j], linestyle=style, color="black", label=lab)
    ax.axvline(TAU, color="0.45", lw=0.8)
    ax.set_xlabel("toy time")
    ax.set_ylabel("node burden")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "fig_4_1_trajectories.png")
    plt.close(fig)

    # Figure 4-2. Spectra.
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    markers = {"nodes": "o", "mode_blind": "s", "mode_partial_proximal": "D", "mode_aware": "^"}
    show = ("nodes", "mode_blind", "mode_partial_proximal", "mode_aware")
    for name in show:
        sv = np.array(reports[name]["singular_values"])
        ax.semilogy(
            np.arange(1, sv.size + 1),
            sv,
            marker=markers[name],
            color="black",
            lw=0.8,
            label=name.replace("_", " "),
        )
    ax.axhline(2.0, color="0.45", lw=0.7, ls="--")
    ax.set_xlabel("principal component")
    ax.set_ylabel("weighted singular value")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "fig_4_2_spectra.png")
    plt.close(fig)

    # Figure 4-3. Marginal relative SE of the four shedding rates.
    fig, ax = plt.subplots(figsize=(6.6, 3.8))
    rate_names = ["lam1", "lam1_quiescent", "lam2", "lam2_quiescent"]
    sched_show = ("mode_blind", "mode_partial_proximal", "mode_partial_distal", "mode_aware")
    x = np.arange(len(rate_names))
    width = 0.18
    for i, sch in enumerate(sched_show):
        vals = []
        for row in reports[sch]["edge_rates"]:
            if row["name"] in rate_names:
                vals.append(row)
        vals = sorted(vals, key=lambda r: rate_names.index(r["name"]))
        heights = []
        for r in vals:
            se = r["marginal_rel_se"]
            heights.append(2.0 if se is None or not np.isfinite(se) else min(se, 2.0))
        ax.bar(x + (i - 1.5) * width, heights, width=width, color=str(0.15 + 0.2 * i), label=sch.replace("_", " "))
    ax.axhline(PRACTICAL_REL, color="0.2", lw=0.7, ls="--")
    ax.set_xticks(x)
    ax.set_xticklabels(["λ1 active", "λ1 quiescent", "λ2 active", "λ2 quiescent"])
    ax.set_ylabel("marginal relative SE (capped at 2)")
    ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(FIG / "fig_4_3_marginal.png")
    plt.close(fig)

    # Figure 4-4. Profiles of λ2.
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    styles = {
        "mode_blind": ("--", "o"),
        "mode_partial_proximal": ("-.", "s"),
        "mode_partial_distal": (":", "D"),
        "mode_aware": ("-", "^"),
    }
    for prof in profiles:
        if prof["parameter"] != "lam2":
            continue
        if prof["schedule"] not in styles:
            continue
        ms = [r["multiplier"] for r in prof["grid"]]
        cs = [max(r["weighted_rss"], 1e-6) for r in prof["grid"]]
        style, marker = styles[prof["schedule"]]
        ax.semilogy(
            ms,
            cs,
            linestyle=style,
            color="black",
            marker=marker,
            ms=3.5,
            label=prof["schedule"].replace("_", " "),
        )
    ax.axhline(PROFILE_CHI, color="0.45", lw=0.7, ls="--")
    ax.set_xlabel("multiplier of active λ on filter→liver")
    ax.set_ylabel("profiled weighted residual")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "fig_4_4_profile.png")
    plt.close(fig)


def generating_table() -> list[dict]:
    rows = []
    for edge, offset in ((1, 0), (2, 3)):
        lam, kap, alpha = (float(TRUTH[offset]), float(TRUTH[offset + 1]), float(TRUTH[offset + 2]))
        law = edge_laws(lam, kap, alpha)
        rows.append(
            {
                "edge": edge,
                "name": "primary→filter" if edge == 1 else "filter→liver",
                "lam_active": lam,
                "kap": kap,
                "alpha": alpha,
                "lam_quiescent": law["lam_q"],
                "phi_active": law["phi_a"],
                "phi_quiescent": law["phi_q"],
                "s_active": law["s_a"],
                "s_quiescent": law["s_q"],
                "s_pool": pool(law),
                "regime_active": "transport-limited" if kap < lam else "shedding-limited",
            }
        )
    return rows


def main() -> None:
    laws = laws_from_theta(THETA0)
    path = trajectory(
        np.array([laws[0]["phi_a"], laws[1]["phi_a"]]),
        np.array([laws[0]["phi_q"], laws[1]["phi_q"]]),
    )
    print("integrating sensitivities", flush=True)
    reports = {name: fisher_report(THETA0, name) for name in SCHEDULES}
    twin = hyperbola_twin(THETA0)
    algebraic = algebraic_jacobian(THETA0)
    sanity(laws, twin, reports)
    print("profiling", flush=True)
    profiles = []
    for schedule in (
        "mode_blind",
        "mode_partial_proximal",
        "mode_partial_distal",
        "mode_aware",
    ):
        print(f"  profile {schedule} lam1", flush=True)
        profiles.append(profile_parameter(schedule, 0, PROFILE_MULTS))
        print(f"  profile {schedule} lam2", flush=True)
        profiles.append(profile_parameter(schedule, 3, PROFILE_FINE))
    print("continuous alias", flush=True)
    alias = continuous_alias()
    print("noisy cloud", flush=True)
    cloud = noisy_cloud()
    make_figures(path, reports, profiles)

    horizon = {
        "primary_at_tau": float(path[SAMPLE_TIMES == TAU][0, 0]),
        "filter_at_end": float(path[-1, 1]),
        "liver_at_end": float(path[-1, 2]),
        "liver_at_tau": float(path[SAMPLE_TIMES == TAU][0, 2]),
        "min_node": float(path.min()),
    }
    payload = {
        "seed": SEED,
        "tau": TAU,
        "horizon": T_END,
        "sample_times": SAMPLE_TIMES.tolist(),
        "n_samples": int(SAMPLE_TIMES.size),
        "soil_r": R_SOIL.tolist(),
        "soil_K": K_SOIL.tolist(),
        "x0": X0.tolist(),
        "sigma_node": SIGMA_NODE,
        "sigma_stall": SIGMA_STALL,
        "fd_eps": FD_EPS,
        "structural_floor": STRUCT_FLOOR,
        "practical_rel_se_cutoff": PRACTICAL_REL,
        "profile_chi": PROFILE_CHI,
        "pool_weights": {"active": W_ACTIVE, "quiescent": W_QUIET},
        "truth": TRUTH.tolist(),
        "names": NAMES,
        "edges": generating_table(),
        "horizon_nodes": horizon,
        "twin": twin,
        "algebraic_rank": algebraic,
        "schedules": reports,
        "profiles": profiles,
        "continuous_alias": alias,
        "noisy_cloud": cloud,
        "calls": {
            sch: {
                "structural_rank": reports[sch]["structural_rank"],
                "practical_rank": reports[sch]["practical_rank"],
                "edge_rates": reports[sch]["edge_rates"],
            }
            for sch in SCHEDULES
        },
    }
    text = json.dumps(json_ready(payload), indent=2)
    (ROOT / "results.json").write_text(text + "\n", encoding="utf-8")
    print("wrote", ROOT / "results.json")
    for sch in SCHEDULES:
        rep = reports[sch]
        flags = ", ".join(
            f"{r['name']}={'stay' if r['stays'] else 'free'}({r['marginal_rel_se']})"
            for r in rep["edge_rates"]
        )
        print(
            f"{sch}: structural {rep['structural_rank']} practical {rep['practical_rank']} | {flags}"
        )


if __name__ == "__main__":
    main()
