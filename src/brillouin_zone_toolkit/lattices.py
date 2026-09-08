"""Common three-dimensional Bravais lattices.

All vectors are returned as rows of a floating-point ``(3, 3)`` array.  The
FCC and BCC constructors intentionally return primitive cells, while ``a``
is the conventional cubic edge.  The hexagonal constructor also returns a
primitive cell and uses the crystallographic convention ``a_1 · a_2 < 0``.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np


def _validated_vectors(vectors: Any) -> np.ndarray:
    """Return a validated copy of three linearly independent vectors."""

    array = np.asarray(vectors, dtype=float)
    if array.shape != (3, 3):
        raise ValueError("primitive_vectors must have shape (3, 3)")
    if not np.all(np.isfinite(array)):
        raise ValueError("primitive_vectors must contain only finite values")
    volume = float(np.linalg.det(array))
    if abs(volume) <= np.finfo(float).eps * max(1.0, float(np.linalg.norm(array)) ** 3):
        raise ValueError("primitive_vectors must be linearly independent")
    return array


@dataclass(frozen=True, eq=False)
class Lattice:
    """A Bravais lattice represented by primitive vectors.

    Parameters
    ----------
    vectors:
        Primitive direct-lattice vectors as rows.
    name:
        Human-readable identifier used in reports and plots.
    """

    vectors: np.ndarray
    name: str = "custom"

    def __post_init__(self) -> None:
        # Keep ownership of the array: freezing a caller's array in place is
        # surprising and makes otherwise unrelated code fail on assignment.
        vectors = _validated_vectors(self.vectors).copy()
        vectors.setflags(write=False)
        object.__setattr__(self, "vectors", vectors)
        if not isinstance(self.name, str) or not self.name.strip():
            raise ValueError("name must be a non-empty string")

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Lattice):
            return NotImplemented
        return self.name == other.name and np.array_equal(self.vectors, other.vectors)

    @property
    def direct_vectors(self) -> np.ndarray:
        """Primitive direct vectors as rows."""

        return self.vectors

    @property
    def volume(self) -> float:
        """Volume of the primitive direct cell."""

        return abs(float(np.linalg.det(self.vectors)))

    @property
    def reciprocal_vectors(self) -> np.ndarray:
        """Primitive reciprocal vectors satisfying ``A Bᵀ = 2π I``."""

        from .reciprocal import reciprocal_vectors

        return reciprocal_vectors(self.vectors)

    @property
    def reciprocal_volume(self) -> float:
        """Volume of the primitive reciprocal cell."""

        return abs(float(np.linalg.det(self.reciprocal_vectors)))

    def copy(self, *, name: str | None = None) -> Lattice:
        """Return an independent lattice object, optionally renamed."""

        return Lattice(self.vectors.copy(), self.name if name is None else name)


def _positive_length(value: float, label: str) -> float:
    value = float(value)
    if not np.isfinite(value) or value <= 0:
        raise ValueError(f"{label} must be a positive finite number")
    return value


def simple_cubic(a: float = 1.0) -> Lattice:
    """Simple-cubic lattice with conventional edge ``a``."""

    a = _positive_length(a, "a")
    return Lattice(a * np.eye(3), "simple cubic (SC)")


def face_centered_cubic(a: float = 1.0) -> Lattice:
    """Face-centred cubic lattice in a primitive basis.

    ``a`` is the conventional cubic edge, so the primitive volume is
    ``a**3 / 4``.
    """

    a = _positive_length(a, "a")
    vectors = 0.5 * a * np.array(
        [[0.0, 1.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 0.0]]
    )
    return Lattice(vectors, "face-centred cubic (FCC)")


def fcc(a: float = 1.0) -> Lattice:
    """Alias for :func:`face_centered_cubic`."""

    return face_centered_cubic(a)


def body_centered_cubic(a: float = 1.0) -> Lattice:
    """Body-centred cubic lattice in a primitive basis.

    ``a`` is the conventional cubic edge, so the primitive volume is
    ``a**3 / 2``.
    """

    a = _positive_length(a, "a")
    vectors = 0.5 * a * np.array(
        [[-1.0, 1.0, 1.0], [1.0, -1.0, 1.0], [1.0, 1.0, -1.0]]
    )
    return Lattice(vectors, "body-centred cubic (BCC)")


def bcc(a: float = 1.0) -> Lattice:
    """Alias for :func:`body_centered_cubic`."""

    return body_centered_cubic(a)


def hexagonal(a: float = 1.0, c: float | None = None) -> Lattice:
    """Hexagonal lattice in a primitive basis.

    Parameters
    ----------
    a:
        Basal-plane lattice constant.
    c:
        Out-of-plane lattice constant.  The ideal close-packed ratio
        ``c/a = sqrt(8/3)`` is used by default.
    """

    a = _positive_length(a, "a")
    c = _positive_length(np.sqrt(8.0 / 3.0) * a if c is None else c, "c")
    vectors = np.array(
        [[a, 0.0, 0.0], [-0.5 * a, 0.5 * np.sqrt(3.0) * a, 0.0], [0.0, 0.0, c]]
    )
    return Lattice(vectors, "hexagonal")


def make_lattice(kind: str, **parameters: Any) -> Lattice:
    """Construct a named lattice.

    Accepted names include ``sc``, ``fcc``, ``bcc`` and ``hex`` (plus their
    long English forms).  Names are case-insensitive and separators are
    ignored.
    """

    if not isinstance(kind, str):
        raise TypeError("kind must be a string")
    normalized = "".join(character for character in kind.lower() if character.isalnum())
    constructors: Mapping[str, Any] = {
        "sc": simple_cubic,
        "simplecubic": simple_cubic,
        "fcc": face_centered_cubic,
        "facecenteredcubic": face_centered_cubic,
        "bcc": body_centered_cubic,
        "bodycenteredcubic": body_centered_cubic,
        "hex": hexagonal,
        "hexagonal": hexagonal,
    }
    try:
        constructor = constructors[normalized]
    except KeyError as error:
        choices = ", ".join(sorted({"sc", "fcc", "bcc", "hex"}))
        raise ValueError(f"unknown lattice {kind!r}; choose one of {choices}") from error
    return constructor(**parameters)
