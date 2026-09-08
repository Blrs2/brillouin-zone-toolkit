import unittest

import numpy as np

from brillouin_zone_toolkit import (
    bcc,
    fcc,
    first_brillouin_zone,
    hexagonal,
    lattice_symmetry_operations,
    reduce_points_with_orbits,
    simple_cubic,
)


class GeometryTests(unittest.TestCase):
    def test_reciprocal_duality(self) -> None:
        for lattice in (simple_cubic(), fcc(), bcc(), hexagonal()):
            np.testing.assert_allclose(
                lattice.direct_vectors @ lattice.reciprocal_vectors.T,
                2.0 * np.pi * np.eye(3),
                atol=1e-11,
            )

    def test_builtin_polyhedra(self) -> None:
        expected = {
            "simple cubic (SC)": (8, 6),
            "face-centred cubic (FCC)": (24, 14),
            "body-centred cubic (BCC)": (14, 12),
            "hexagonal": (12, 8),
        }
        for lattice in (simple_cubic(), fcc(), bcc(), hexagonal()):
            zone = first_brillouin_zone(lattice, search=2)
            self.assertEqual((zone.n_vertices, zone.n_faces), expected[lattice.name])
            self.assertAlmostEqual(zone.volume, lattice.reciprocal_volume, places=8)
            self.assertTrue(np.all(zone.contains(zone.vertices)))

    def test_half_spaces_are_satisfied(self) -> None:
        zone = first_brillouin_zone(fcc(), search=2)
        residual = zone.vertices @ zone.halfspaces[:, :3].T + zone.halfspaces[:, 3]
        self.assertLessEqual(float(np.max(residual)), 1e-8)
        self.assertTrue(zone.contains([0.0, 0.0, 0.0]))
        self.assertFalse(zone.contains([20.0, 0.0, 0.0]))

    def test_symmetry_operations_and_orbits(self) -> None:
        lattice = fcc()
        operations = lattice_symmetry_operations(lattice)
        self.assertEqual(len(operations), 48)
        points = np.array([[1.0, 2.0, 3.0], [-1.0, 2.0, 3.0], [3.0, 2.0, 1.0]])
        reduced = reduce_points_with_orbits(points, operations)
        self.assertEqual(reduced.n_orbits, 1)
        self.assertEqual(int(reduced.multiplicities[0]), len(points))


if __name__ == "__main__":
    unittest.main()
