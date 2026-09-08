"""Numerical tools for reciprocal lattices and first Brillouin zones.

The package uses row vectors throughout: a ``(3, 3)`` array contains the three
primitive vectors as its rows.  Reciprocal vectors follow the convention
``a_i · b_j = 2π δ_ij``.
"""

from .brillouin import BrillouinZone, first_brillouin_zone
from .lattices import (
    Lattice,
    bcc,
    body_centered_cubic,
    face_centered_cubic,
    fcc,
    hexagonal,
    make_lattice,
    simple_cubic,
)
from .reciprocal import (
    direct_reciprocal_duality,
    reciprocal_cell_volume,
    reciprocal_lattice_points,
    reciprocal_vectors,
)
from .symmetry import (
    SymmetryReduction,
    irreducible_wedge,
    lattice_symmetry_operations,
    reduce_points,
    reduce_points_with_orbits,
    symmetry_orbit,
    symmetry_reduce,
)

__version__ = "0.1.0"

__all__ = [
    "BrillouinZone",
    "Lattice",
    "SymmetryReduction",
    "__version__",
    "bcc",
    "body_centered_cubic",
    "direct_reciprocal_duality",
    "face_centered_cubic",
    "fcc",
    "first_brillouin_zone",
    "hexagonal",
    "irreducible_wedge",
    "lattice_symmetry_operations",
    "make_lattice",
    "reciprocal_cell_volume",
    "reciprocal_lattice_points",
    "reciprocal_vectors",
    "reduce_points",
    "reduce_points_with_orbits",
    "simple_cubic",
    "symmetry_orbit",
    "symmetry_reduce",
]
