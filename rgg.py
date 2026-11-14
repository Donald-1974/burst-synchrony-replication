"""Shared helpers: 3D periodic random geometric graphs and row-normalised Kuramoto dynamics.

All functions are deterministic given a ``numpy.random.Generator``.
"""
from __future__ import annotations

import numpy as np
from scipy.sparse import csr_matrix
from scipy.spatial import cKDTree


def periodic_rgg_pairs(n: int, r0: float, rng: np.random.Generator) -> np.ndarray:
    """Place ``n`` points uniformly in the unit cube with periodic boundaries and
    return all unordered pairs closer than ``r0`` as an ``(m, 2)`` int array (i < j)."""
    pos = rng.random((n, 3))
    pairs = cKDTree(pos, boxsize=1.0).query_pairs(r0, output_type="ndarray")
    return pairs.reshape(-1, 2).astype(np.int64)


def symmetric_adjacency(pairs: np.ndarray, n: int) -> csr_matrix:
    """Unweighted symmetric adjacency matrix (CSR) from an edge list of unordered pairs."""
    rows = np.concatenate([pairs[:, 0], pairs[:, 1]])
    cols = np.concatenate([pairs[:, 1], pairs[:, 0]])
    return csr_matrix((np.ones(rows.size), (rows, cols)), shape=(n, n))


def row_normalise(adj: csr_matrix) -> tuple[csr_matrix, np.ndarray]:
    """Return (D^-1 A, degree vector). Isolated nodes keep an all-zero row
    (their degree is treated as 1 only to avoid division by zero)."""
    deg = np.asarray(adj.sum(axis=1)).ravel()
    inv = 1.0 / np.maximum(deg, 1.0)
    return csr_matrix(adj.multiply(inv[:, None])), deg


def kuramoto_velocity(theta: np.ndarray, omega: np.ndarray, coupling: float,
                      adj_norm: csr_matrix) -> np.ndarray:
    """dθ_i/dt = ω_i + K · (1/k_i) Σ_j A_ij sin(θ_j − θ_i), evaluated with complex phasors."""
    z = np.exp(1j * theta)
    return omega + coupling * np.imag(np.conj(z) * (adj_norm @ z))


def order_parameter(theta: np.ndarray) -> float:
    """Global Kuramoto order parameter r = |⟨exp(iθ)⟩|."""
    return float(np.abs(np.exp(1j * theta).mean()))
