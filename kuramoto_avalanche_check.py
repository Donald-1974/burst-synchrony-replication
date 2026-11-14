"""Sanity check of the 'Kuramoto avalanche layer' specification (Dec 2025 conversations).

Model: N = 2,000 oscillators in a 3D periodic random geometric graph (r0 = 0.1, mean degree
≈ 8.5), row-normalised Kuramoto coupling, ω ~ N(0, 1), forward-Euler dt = 0.01.
'Spikes' are forward crossings of phase 0 (mod 2π); an avalanche is a run of consecutive
time bins each containing ≥ 1 spike, bounded by empty bins. A run still open at the end of
the record is censored and not counted.

The point of the check: with thousands of independently rotating oscillators, nearly every
bin contains a spike unless bins are very short, so avalanche statistics depend mainly on
the arbitrary bin width rather than on the dynamics. Oscillators with negative ω rotate
backwards and are only counted if coupling drags them forward — a limitation inherited
from the spec's spike definition.

Usage:  python kuramoto_avalanche_check.py [--n 2000] [--r0 0.1] [--seed 42]
"""
from __future__ import annotations

import argparse

import numpy as np

from rgg import (kuramoto_velocity, order_parameter, periodic_rgg_pairs, row_normalise,
                 symmetric_adjacency)

COUPLINGS = [0.5, 1.0, 1.5, 2.0, 3.0, 5.0]
DT = 0.01
T_TRANSIENT = 50.0
T_MEASURE = 100.0
BIN_STEPS = [1, 10]  # bin widths in integration steps (0.01 and 0.10 time units)


def avalanche_sizes(counts: np.ndarray) -> list[int]:
    """Sizes (total spikes) of runs of non-empty bins that are terminated by an empty bin."""
    sizes: list[int] = []
    current = 0
    for c in counts:
        if c > 0:
            current += int(c)
        elif current > 0:
            sizes.append(current)
            current = 0
    return sizes


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--n", type=int, default=2000)
    ap.add_argument("--r0", type=float, default=0.1)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    pairs = periodic_rgg_pairs(args.n, args.r0, rng)
    adj_norm, deg = row_normalise(symmetric_adjacency(pairs, args.n))
    omega = rng.normal(0.0, 1.0, args.n)
    print(f"N={args.n} mean degree={deg.mean():.1f}")

    two_pi = 2 * np.pi
    for coupling in COUPLINGS:
        theta = rng.uniform(0.0, two_pi, args.n)
        for _ in range(int(round(T_TRANSIENT / DT))):
            theta = theta + DT * kuramoto_velocity(theta, omega, coupling, adj_norm)

        n_meas = int(round(T_MEASURE / DT))
        spikes = np.empty(n_meas, dtype=np.int64)
        rs = np.empty(n_meas)
        for t in range(n_meas):
            new = theta + DT * kuramoto_velocity(theta, omega, coupling, adj_norm)
            spikes[t] = np.count_nonzero(np.floor(new / two_pi) > np.floor(theta / two_pi))
            theta = new
            rs[t] = order_parameter(theta)

        for b in BIN_STEPS:
            counts = spikes[: len(spikes) // b * b].reshape(-1, b).sum(axis=1)
            sizes = avalanche_sizes(counts)
            print(f"K={coupling:4.1f} r={rs.mean():.3f} bin={b * DT:.2f} "
                  f"empty_bins={np.mean(counts == 0):.3f} n_avalanches={len(sizes)} "
                  f"max={max(sizes) if sizes else 0}")


if __name__ == "__main__":
    main()
