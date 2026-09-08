"""Reciprocal-lattice operations."""

from __future__ import annotations

from itertools import product
from typing import Any

import numpy as np


def _vectors(vectors: Any, label: str = "vectors") -> np.ndarray:
    array = np.asarray(vectors, dtype=float)
    if array.shape != (3, 3):
        raise ValueError(f"{label} must have shape (3, 3)")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{label} must contain only finite values")
    if abs(float(np.linalg.det(array))) <= np.finfo(float).eps * max(1.0, np.linalg.norm(array) ** 3):
        raise ValueError(f"{label} must be linearly independent")
    return array


def reciprocal_vectors(direct_vectors: Any) -> np.ndarray:
    """Compute primitive reciprocal vectors from direct vectors.

    If ``A`` contains direct vectors as rows, the returned ``B`` obeys
    ``A @ B.T == 2π * eye(3)``.  This is equivalent to the cross-product
    construction in the report, with the conventional ``2π`` normalization.
    """

    direct = _vectors(direct_vectors, "direct_vectors")
    # Solving the transposed system is both clearer and numerically preferable
    # to explicitly constructing an inverse.  For rows A, A B^T = 2*pi I.
    return 2.0 * np.pi * np.linalg.solve(direct, np.eye(3)).T


def direct_reciprocal_duality(direct_vectors: Any, reciprocal: Any | None = None) -> np.ndarray:
    """Return the duality matrix ``A @ B.T`` for a pair of bases."""

    direct = _vectors(direct_vectors, "direct_vectors")
    reciprocal_array = reciprocal_vectors(direct) if reciprocal is None else _vectors(reciprocal, "reciprocal")
    return direct @ reciprocal_array.T


def reciprocal_lattice_points(
    reciprocal_basis: Any,
    max_index: int = 2,
    *,
    include_zero: bool = False,
    sort_by_norm: bool = True,
) -> np.ndarray:
    """Enumerate reciprocal-lattice points around the origin.

    Integer triples ``n`` in ``[-max_index, max_index]^3`` generate points
    ``n[0] b1 + n[1] b2 + n[2] b3``.  The result is sorted by norm by default,
    with the origin first when requested.
    """

    basis = _vectors(reciprocal_basis, "reciprocal_basis")
    if isinstance(max_index, (bool, np.bool_)) or not isinstance(max_index, (int, np.integer)) or max_index < 0:
        raise ValueError("max_index must be a non-negative integer")
    coefficients = np.asarray(
        list(product(range(-int(max_index), int(max_index) + 1), repeat=3)), dtype=int
    )
    if not include_zero:
        coefficients = coefficients[np.any(coefficients != 0, axis=1)]
    points = coefficients @ basis
    if sort_by_norm and len(points):
        order = np.argsort(np.linalg.norm(points, axis=1), kind="stable")
        points = points[order]
    return points


def reciprocal_cell_volume(direct_vectors: Any) -> float:
    """Return the reciprocal primitive-cell volume ``(2π)^3 / V``."""

    direct = _vectors(direct_vectors, "direct_vectors")
    return (2.0 * np.pi) ** 3 / abs(float(np.linalg.det(direct)))
