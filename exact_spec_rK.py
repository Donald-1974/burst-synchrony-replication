"""Re-run the Dec 20, 2025 'Hyperplastic Criticality' r(K) specification exactly as written.

Spec: N = 10,000 oscillators placed uniformly in a 3D periodic unit box, linked when closer
than r0 = 0.03; row-normalised Kuramoto coupling; natural frequencies ω ~ N(0, 1);
forward-Euler integration (dt = 0.02) to t = 150, averaging r over t > 100.
The original conversation printed an 'Expected output' of r(K=0.85) ≈ 0.378, i.e. a
partially synchronised state. This script checks that claim.

Usage:  python exact_spec_rK.py [--n 10000] [--r0 0.03] [--seed 42]
"""
from __future__ import annotations

import argparse

import numpy as np

from rgg import (kuramoto_velocity, order_parameter, periodic_rgg_pairs, row_normalise,
                 symmetric_adjacency)

COUPLINGS = [0.5, 0.6, 0.68, 0.75, 0.85, 2.0, 5.0]
DT = 0.02          # Euler step; |dθ/dt| ≤ |ω| + K, so dt·K ≤ 0.1 for K ≤ 5
T_TOTAL = 150.0    # total simulated time
T_TRANSIENT = 100.0  # discard t ≤ this before averaging r


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--n", type=int, default=10_000)
    ap.add_argument("--r0", type=float, default=0.03)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    pairs = periodic_rgg_pairs(args.n, args.r0, rng)
    adj_norm, deg = row_normalise(symmetric_adjacency(pairs, args.n))
    omega = rng.normal(0.0, 1.0, args.n)
    print(f"edges={len(pairs)} <k>={deg.mean():.2f}")

    n_steps = int(round(T_TOTAL / DT))
    for coupling in COUPLINGS:
        theta = rng.uniform(0.0, 2 * np.pi, args.n)
        rs = []
        for step in range(n_steps):
            theta = theta + DT * kuramoto_velocity(theta, omega, coupling, adj_norm)
            if step * DT > T_TRANSIENT:
                rs.append(order_parameter(theta))
        print(f"K={coupling:4.2f}  <r>={np.mean(rs):.3f}  std={np.std(rs):.3f}")


if __name__ == "__main__":
    main()
