"""Fine coupling sweep of row-normalised Kuramoto dynamics on a *connected* 3D periodic RGG.

This is the corrected follow-up to ``exact_spec_rK.py``: the network is well above the
percolation threshold (default N = 2,000, r0 = 0.1, mean degree ~8.5) and the dynamics are
run on the giant component only (any isolated nodes or small fragments are dropped and
reported). Seeds run in parallel worker processes.

For each seed (independent graph, frequencies and initial phases) and each coupling K the
script integrates

    dθ_i/dt = ω_i + K (1/k_i) Σ_j A_ij sin(θ_j − θ_i),   ω_i ~ N(0, 1),

with Heun's method (all K values are integrated together as columns of one matrix), and
records the time-averaged order parameter ⟨r⟩ and the susceptibility χ = N·var(r). The
critical coupling is estimated from the χ peak of the seed-averaged curve.

Avalanches. A 'spike' is a crossing of phase 0 (mod 2π) in *either* direction, so
oscillators with negative ω are counted. Crossing times are located inside each step by
linear interpolation, then binned with a width equal to the population mean inter-event
interval (Beggs & Plenz 2003). An avalanche is a run of non-empty bins bounded by empty
bins; runs touching the ends of the record are discarded. Size distributions are fitted by
discrete power-law maximum likelihood with xmin chosen by minimum KS distance (Clauset,
Shalizi & Newman 2009) and compared with a discrete exponential (geometric) tail by the
normalised log-likelihood ratio R and its Vuong p-value (R > 0 favours the power law).

Usage:  python connected_sweep.py [--n 2000] [--r0 0.1] [--seeds 5] [--dt 0.02]
                                  [--kmin 0] [--kmax 8] [--nk 41] [--outdir results]
"""
from __future__ import annotations

import argparse
import math
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.sparse.csgraph import connected_components
from scipy.special import zeta
from scipy.stats import norm

from rgg import periodic_rgg_pairs, row_normalise, symmetric_adjacency

T_TRANSIENT = 100.0
T_MEASURE = 200.0
TWO_PI = 2.0 * np.pi


def velocity(theta: np.ndarray, omega: np.ndarray, ks: np.ndarray, adj_norm) -> np.ndarray:
    """Kuramoto velocity for a (N, nK) phase matrix; column c uses coupling ks[c]."""
    z = np.exp(1j * theta)
    return omega[:, None] + ks[None, :] * np.imag(np.conj(z) * (adj_norm @ z))


def simulate(n: int, r0: float, ks: np.ndarray, seed: int, dt: float,
             record_spikes: bool) -> dict:
    """Run one seed for all couplings. Returns r time-series stats and optional spike times."""
    rng = np.random.default_rng(seed)
    pairs = periodic_rgg_pairs(n, r0, rng)
    adj = symmetric_adjacency(pairs, n)
    n_comp, labels = connected_components(adj, directed=False)
    keep = labels == np.argmax(np.bincount(labels))
    adj = adj[keep][:, keep].tocsr()
    n_total, n = n, int(keep.sum())
    adj_norm, deg = row_normalise(adj)
    omega = rng.normal(0.0, 1.0, n)
    theta = rng.uniform(0.0, TWO_PI, (n, ks.size))

    def step(th: np.ndarray) -> np.ndarray:
        k1 = velocity(th, omega, ks, adj_norm)
        k2 = velocity(th + dt * k1, omega, ks, adj_norm)
        return th + 0.5 * dt * (k1 + k2)

    for _ in range(round(T_TRANSIENT / dt)):
        theta = step(theta)

    n_meas = round(T_MEASURE / dt)
    r_sum = np.zeros(ks.size)
    r2_sum = np.zeros(ks.size)
    spikes: list[list[np.ndarray]] = [[] for _ in ks] if record_spikes else []
    for t in range(n_meas):
        new = step(theta)
        z = np.exp(1j * new).mean(axis=0)
        r = np.abs(z)
        r_sum += r
        r2_sum += r * r
        if record_spikes:
            w_old = np.floor(theta / TWO_PI)
            w_new = np.floor(new / TWO_PI)
            rows, cols = np.nonzero(w_old != w_new)
            if rows.size:
                # Boundary crossed: 2π·max(w_old, w_new); interpolate the crossing time.
                b = TWO_PI * np.maximum(w_old[rows, cols], w_new[rows, cols])
                a0, a1 = theta[rows, cols], new[rows, cols]
                frac = np.clip((b - a0) / (a1 - a0), 0.0, 1.0)
                times = (t + frac) * dt
                for c in np.unique(cols):
                    spikes[c].append(times[cols == c])
        theta = new
    mean_r = r_sum / n_meas
    var_r = np.maximum(r2_sum / n_meas - mean_r**2, 0.0)
    out = {"r": mean_r, "chi": n * var_r, "mean_degree": float(deg.mean()),
           "n_components": int(n_comp), "n_giant": n, "n_total": n_total}
    if record_spikes:
        out["spikes"] = [np.sort(np.concatenate(s)) if s else np.empty(0) for s in spikes]
    return out


def avalanche_sizes(times: np.ndarray, t_total: float) -> tuple[np.ndarray, float]:
    """Bin spike times with width = mean inter-event interval; return interior avalanche sizes."""
    if times.size < 2:
        return np.empty(0, dtype=np.int64), float("nan")
    width = float(np.mean(np.diff(times)))
    counts = np.bincount((times / width).astype(np.int64),
                         minlength=math.ceil(t_total / width))
    nonempty = counts > 0
    edges = np.diff(nonempty.astype(np.int8))
    starts = np.flatnonzero(edges == 1) + 1
    ends = np.flatnonzero(edges == -1) + 1
    if nonempty[0]:
        ends = ends[1:]  # drop run touching the start
    starts = starts[: ends.size]  # drop run touching the end
    csum = np.concatenate([[0], np.cumsum(counts)])
    return (csum[ends] - csum[starts]).astype(np.int64), width


def fit_power_law(sizes: np.ndarray, min_tail: int = 50) -> dict | None:
    """Discrete power-law MLE with KS-selected xmin; Vuong test against a geometric tail."""
    sizes = sizes[sizes > 0]
    best = None
    for xmin in np.unique(sizes):
        tail = sizes[sizes >= xmin]
        if tail.size < min_tail:
            break
        ln_sum = np.log(tail).sum()
        grid = np.linspace(1.05, 5.0, 791)
        ll = -grid * ln_sum - tail.size * np.log(zeta(grid, xmin))
        alpha = float(grid[np.argmax(ll)])
        xs = np.arange(xmin, tail.max() + 1)
        pmf = xs ** (-alpha) / zeta(alpha, xmin)
        model_cdf = np.cumsum(pmf)
        emp_cdf = np.searchsorted(np.sort(tail), xs, side="right") / tail.size
        ks_d = float(np.max(np.abs(emp_cdf - model_cdf)))
        if best is None or ks_d < best["ks"]:
            best = {"xmin": int(xmin), "alpha": alpha, "ks": ks_d, "n_tail": int(tail.size)}
    if best is None:
        return None
    tail = sizes[sizes >= best["xmin"]].astype(float)
    a, xmin = best["alpha"], best["xmin"]
    ll_pl = -a * np.log(tail) - np.log(zeta(a, xmin))
    p = 1.0 / (1.0 + np.mean(tail - xmin))  # geometric MLE on tail - xmin
    ll_geo = np.log(p) + (tail - xmin) * np.log1p(-p)
    d = ll_pl - ll_geo
    sd = d.std(ddof=1)
    r_norm = float(d.sum() / (sd * math.sqrt(d.size))) if sd > 0 else 0.0
    best["R_vs_exp"] = float(d.sum())
    best["vuong_z"] = r_norm
    best["vuong_p"] = float(2 * norm.sf(abs(r_norm)))
    return best


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--n", type=int, default=2000)
    ap.add_argument("--r0", type=float, default=0.1)
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--seed0", type=int, default=42)
    ap.add_argument("--dt", type=float, default=0.02)
    ap.add_argument("--kmin", type=float, default=0.0)
    ap.add_argument("--kmax", type=float, default=8.0)
    ap.add_argument("--nk", type=int, default=41)
    ap.add_argument("--outdir", type=Path, default=Path("results"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    ks = np.linspace(args.kmin, args.kmax, args.nk)
    r_all, chi_all, spikes0 = [], [], None
    seeds = [args.seed0 + i for i in range(args.seeds)]
    with ProcessPoolExecutor(max_workers=min(len(seeds), os.cpu_count() or 1)) as pool:
        futures = [pool.submit(simulate, args.n, args.r0, ks, sd, args.dt, i == 0)
                   for i, sd in enumerate(seeds)]
        results = [f.result() for f in futures]
    for sd, res in zip(seeds, results):
        r_all.append(res["r"])
        chi_all.append(res["chi"])
        print(f"seed {sd}: {res['n_components']} component(s), giant component "
              f"{res['n_giant']}/{res['n_total']} nodes, mean degree {res['mean_degree']:.2f}")
    spikes0 = results[0]["spikes"]
    r_all, chi_all = np.array(r_all), np.array(chi_all)
    r_m, r_s = r_all.mean(0), r_all.std(0, ddof=1) if args.seeds > 1 else np.zeros(ks.size)
    chi_m, chi_s = chi_all.mean(0), chi_all.std(0, ddof=1) if args.seeds > 1 else np.zeros(ks.size)
    k_peak = float(ks[np.argmax(chi_m)])
    k_mf = 2.0 / (math.pi * norm.pdf(0.0))

    print(f"\nN={args.n} (before giant-component restriction) r0={args.r0} dt={args.dt} seeds={args.seeds} "
          f"T_transient={T_TRANSIENT} T_measure={T_MEASURE}")
    print(f"{'K':>6} {'<r>':>7} {'sd':>6} {'chi':>8} {'sd':>7}")
    for k, a, b, c, d in zip(ks, r_m, r_s, chi_m, chi_s):
        print(f"{k:6.2f} {a:7.3f} {b:6.3f} {c:8.3f} {d:7.3f}")
    print(f"\nchi peak (seed-averaged) at K = {k_peak:.2f}; "
          f"per-seed peaks: {[round(float(ks[j]), 2) for j in chi_all.argmax(1)]}")
    print(f"mean-field (all-to-all) Kc = 2/(pi g(0)) = {k_mf:.3f} for omega ~ N(0,1)")

    j_peak = int(np.argmax(chi_m))
    probe = sorted({0, j_peak // 2, j_peak, min(j_peak + (args.nk - j_peak) // 2, args.nk - 1),
                    args.nk - 1})
    print(f"\nAvalanches (seed {args.seed0}, bin = mean inter-event interval):")
    aval = {}
    for j in probe:
        sizes, width = avalanche_sizes(spikes0[j], T_MEASURE)
        fit = fit_power_law(sizes)
        aval[float(ks[j])] = sizes
        line = (f"K={ks[j]:5.2f} spikes={spikes0[j].size} bin={width:.5f} "
                f"n_aval={sizes.size} max={sizes.max() if sizes.size else 0}")
        if fit:
            cap = " (at fit-grid limit 5.0: tail too steep for a power law)" \
                if fit["alpha"] >= 4.99 else ""
            line += (f" | PL alpha={fit['alpha']:.2f}{cap} xmin={fit['xmin']} n_tail={fit['n_tail']}"
                     f" KS={fit['ks']:.3f} R_vs_exp={fit['R_vs_exp']:.1f}"
                     f" z={fit['vuong_z']:.2f} p={fit['vuong_p']:.3g}")
        else:
            line += " | too few avalanches for a tail fit"
        print(line)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].errorbar(ks, r_m, yerr=r_s, fmt="o-", ms=3)
    ax[0].set(xlabel="K", ylabel="<r>", title="Order parameter")
    ax[1].errorbar(ks, chi_m, yerr=chi_s, fmt="o-", ms=3)
    ax[1].axvline(k_peak, ls="--", c="gray")
    ax[1].set(xlabel="K", ylabel="chi = N var(r)", title="Susceptibility")
    fig.tight_layout()
    fig.savefig(args.outdir / "r_chi_vs_K.png", dpi=120)
    fig.savefig(args.outdir / "r_chi_vs_K.svg")
    fig2, ax2 = plt.subplots(figsize=(5, 4))
    for k, s in aval.items():
        if s.size:
            v, c = np.unique(s, return_counts=True)
            ccdf = 1.0 - np.concatenate([[0], np.cumsum(c)[:-1]]) / s.size
            ax2.loglog(v, ccdf, ".", ms=3, label=f"K={k:.1f}")
    ax2.set(xlabel="avalanche size s", ylabel="P(S >= s)", title="Avalanche sizes")
    ax2.legend(fontsize=7)
    fig2.tight_layout()
    fig2.savefig(args.outdir / "avalanche_ccdf.png", dpi=120)
    fig2.savefig(args.outdir / "avalanche_ccdf.svg")
    csv_path = args.outdir / "sweep.csv"
    seeds_hdr = [f"r_seed{sd}" for sd in seeds] + [f"chi_seed{sd}" for sd in seeds]
    header = ",".join(["K", "r_mean", "r_sd", "chi_mean", "chi_sd"] + seeds_hdr)
    rows = np.column_stack([ks, r_m, r_s, chi_m, chi_s, r_all.T, chi_all.T])
    np.savetxt(csv_path, rows, delimiter=",", header=header, comments="")
    np.savez_compressed(args.outdir / "sweep.npz", ks=ks, r=r_all, chi=chi_all)


if __name__ == "__main__":
    main()
