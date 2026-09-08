"""Point-group operations and irreducible-point reduction."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from itertools import permutations, product
from typing import Any

import numpy as np

from .lattices import Lattice


def _signed_permutation_matrices() -> list[np.ndarray]:
    matrices: list[np.ndarray] = []
    for permutation in permutations(range(3)):
        for signs in product((-1.0, 1.0), repeat=3):
            matrix = np.zeros((3, 3), dtype=float)
            for row, column in enumerate(permutation):
                matrix[row, column] = signs[row]
            matrices.append(matrix)
    return matrices


def _hexagonal_matrices() -> list[np.ndarray]:
    matrices: list[np.ndarray] = []
    for step in range(6):
        angle = step * np.pi / 3.0
        rotation = np.array(
            [[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]], dtype=float
        )
        for reflection in (np.eye(2), np.diag([1.0, -1.0])):
            basal = rotation @ reflection
            for z_sign in (-1.0, 1.0):
                matrix = np.eye(3)
                matrix[:2, :2] = basal
                matrix[2, 2] = z_sign
                matrices.append(matrix)
    return matrices


def _is_hexagonal(lattice: Lattice) -> bool:
    vectors = lattice.direct_vectors
    lengths = np.linalg.norm(vectors, axis=1)
    scale = max(1.0, float(np.max(lengths)))
    return bool(
        np.isclose(lengths[0], lengths[1], rtol=1e-8, atol=1e-10 * scale)
        and np.isclose(np.dot(vectors[0], vectors[1]), -0.5 * lengths[0] ** 2, rtol=1e-8, atol=1e-10 * scale**2)
        and np.isclose(np.dot(vectors[0], vectors[2]), 0.0, rtol=1e-8, atol=1e-10 * scale**2)
        and np.isclose(np.dot(vectors[1], vectors[2]), 0.0, rtol=1e-8, atol=1e-10 * scale**2)
    )


def _preserves_lattice(lattice: Lattice, operation: np.ndarray, tolerance: float) -> bool:
    """Check that a Cartesian operation maps primitive vectors to the lattice."""

    direct = lattice.direct_vectors
    integer_coordinates = direct @ operation.T @ np.linalg.inv(direct)
    return bool(np.allclose(integer_coordinates, np.rint(integer_coordinates), atol=tolerance, rtol=0.0))


def _unique_matrices(matrices: Iterable[np.ndarray]) -> np.ndarray:
    unique: dict[tuple[float, ...], np.ndarray] = {}
    for matrix in matrices:
        key = tuple(np.round(np.asarray(matrix).ravel(), decimals=10))
        unique.setdefault(key, np.asarray(matrix, dtype=float))
    result = list(unique.values())
    result.sort(key=lambda matrix: tuple(matrix.ravel()))
    identity_index = next(
        (index for index, matrix in enumerate(result) if np.allclose(matrix, np.eye(3))), None
    )
    if identity_index is not None:
        result.insert(0, result.pop(identity_index))
    return np.asarray(result)


def lattice_symmetry_operations(
    lattice: Lattice, *, tolerance: float = 1e-8, include_improper: bool = True
) -> np.ndarray:
    """Return Cartesian point-group matrices preserving ``lattice``.

    Cubic primitive cells are tested against all 48 signed permutation
    matrices.  Hexagonal cells use the 24 operations of ``D6h``.  For a
    custom cell the signed-permutation candidates are filtered, which yields
    a safe subgroup instead of silently applying operations that do not map
    the lattice to itself.
    """

    if not isinstance(lattice, Lattice):
        raise TypeError("lattice_symmetry_operations expects a Lattice")
    if tolerance <= 0 or not np.isfinite(tolerance):
        raise ValueError("tolerance must be a positive finite number")
    candidates = _hexagonal_matrices() if _is_hexagonal(lattice) else _signed_permutation_matrices()
    valid = [matrix for matrix in candidates if _preserves_lattice(lattice, matrix, tolerance)]
    if not include_improper:
        valid = [matrix for matrix in valid if np.linalg.det(matrix) > 0]
    return _unique_matrices(valid)


@dataclass(frozen=True)
class SymmetryReduction:
    """Result of reducing a point cloud into symmetry orbits."""

    representatives: np.ndarray
    labels: np.ndarray
    multiplicities: np.ndarray

    @property
    def n_orbits(self) -> int:
        return len(self.representatives)


def _transformed(point: np.ndarray, operations: np.ndarray) -> np.ndarray:
    # Points are rows while operations act on Cartesian column vectors.
    return np.asarray(point) @ np.transpose(operations, (0, 2, 1))


def symmetry_orbit(point: Any, operations: Any) -> np.ndarray:
    """Return all transformed copies of one Cartesian point."""

    point_array = np.asarray(point, dtype=float)
    if point_array.shape != (3,):
        raise ValueError("point must have shape (3,)")
    matrices = np.asarray(operations, dtype=float)
    if matrices.ndim != 3 or matrices.shape[1:] != (3, 3):
        raise ValueError("operations must have shape (n, 3, 3)")
    return _transformed(point_array, matrices)


def reduce_points_with_orbits(
    points: Any, operations: Any, *, tolerance: float = 1e-8
) -> SymmetryReduction:
    """Reduce points by choosing the lexicographically smallest orbit member.

    The returned labels map every input row to its representative index.  A
    point set should normally be restricted to a symmetry-invariant region
    (for example the first Brillouin zone) before reduction.
    """

    array = np.asarray(points, dtype=float)
    if array.ndim == 1:
        array = array[None, :]
    if array.ndim != 2 or array.shape[1] != 3:
        raise ValueError("points must have shape (n, 3)")
    matrices = np.asarray(operations, dtype=float)
    if matrices.ndim != 3 or matrices.shape[1:] != (3, 3) or len(matrices) == 0:
        raise ValueError("operations must have shape (n, 3, 3) with n > 0")
    if tolerance <= 0 or not np.isfinite(tolerance):
        raise ValueError("tolerance must be a positive finite number")
    representatives: list[np.ndarray] = []
    labels = np.empty(len(array), dtype=int)
    multiplicities: list[int] = []
    index_by_key: dict[tuple[int, int, int], int] = {}
    scale = max(1.0, float(np.max(np.linalg.norm(array, axis=1))) if len(array) else 1.0)
    quantum = tolerance * scale
    for input_index, point in enumerate(array):
        orbit = _transformed(point, matrices)
        canonical = orbit[np.lexsort((orbit[:, 2], orbit[:, 1], orbit[:, 0]))[0]]
        key = tuple(np.rint(canonical / quantum).astype(np.int64))
        representative_index = index_by_key.get(key)
        if representative_index is None:
            representative_index = len(representatives)
            index_by_key[key] = representative_index
            representatives.append(canonical)
            multiplicities.append(0)
        labels[input_index] = representative_index
        multiplicities[representative_index] += 1
    return SymmetryReduction(
        representatives=np.asarray(representatives, dtype=float).reshape((-1, 3)),
        labels=labels,
        multiplicities=np.asarray(multiplicities, dtype=int),
    )


def reduce_points(
    points: Any, operations: Any, *, tolerance: float = 1e-8
) -> np.ndarray:
    """Return only one representative per symmetry orbit."""

    return reduce_points_with_orbits(points, operations, tolerance=tolerance).representatives


def symmetry_reduce(
    points: Any, operations: Any, *, tolerance: float = 1e-8
) -> np.ndarray:
    """Alias for :func:`reduce_points`, matching the report terminology."""

    return reduce_points(points, operations, tolerance=tolerance)


def irreducible_wedge(points: Any, lattice: Lattice, *, tolerance: float = 1e-8) -> np.ndarray:
    """Reduce a point cloud with the point group of ``lattice``."""

    return reduce_points(points, lattice_symmetry_operations(lattice, tolerance=tolerance), tolerance=tolerance)
