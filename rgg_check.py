"""Connectivity of 3D periodic random geometric graphs (RGGs) at the radii used in the specs.

For N nodes uniform in the unit cube (periodic boundaries) with connection radius r0, the
expected mean degree is N·(4/3)π·r0³. This script reports the measured mean degree,
number of connected components, giant-component fraction and isolated-node fraction.
It shows that r0 = 0.03 with N = 10,000 (the Dec 2025 spec) is far below the
percolation threshold, so that network is fragmented.

Usage:  python rgg_check.py [--n 10000] [--radii 0.03 0.05 0.07] [--seed 42]
"""
from __future__ import annotations

import argparse

import numpy as np
from scipy.sparse.csgraph import connected_components

from rgg import periodic_rgg_pairs, symmetric_adjacency


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--n", type=int, default=10_000, help="number of nodes")
    ap.add_argument("--radii", type=float, nargs="+", default=[0.03, 0.05, 0.07])
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    n = args.n
    for r0 in args.radii:
        pairs = periodic_rgg_pairs(n, r0, rng)
        adj = symmetric_adjacency(pairs, n)
        n_comp, labels = connected_components(adj, directed=False)
        giant = np.bincount(labels).max() / n
        deg = np.asarray(adj.sum(axis=1)).ravel()
        expected = n * 4.0 / 3.0 * np.pi * r0**3
        print(f"N={n} r0={r0}: mean degree={deg.mean():.2f} (expected {expected:.2f}), "
              f"components={n_comp}, giant component={giant:.1%}, isolated={np.mean(deg == 0):.1%}")


if __name__ == "__main__":
    main()
