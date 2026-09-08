import unittest

import numpy as np

from brillouin_zone_toolkit import (
    bcc,
    fcc,
    hexagonal,
    lattice_symmetry_operations,
    reduce_points_with_orbits,
    simple_cubic,
    symmetry_orbit,
)


class SymmetryTests(unittest.TestCase):
    def test_point_group_orders(self):
        for lattice, order in (
            (simple_cubic(), 48),
            (fcc(), 48),
            (bcc(), 48),
            (hexagonal(), 24),
        ):
            operations = lattice_symmetry_operations(lattice)
            self.assertEqual(operations.shape, (order, 3, 3))
            for operation in operations:
                np.testing.assert_allclose(operation @ operation.T, np.eye(3), atol=1e-10)
            proper = lattice_symmetry_operations(lattice, include_improper=False)
            self.assertEqual(len(proper), order // 2)

    def test_reduction_returns_labels_and_multiplicities(self):
        operations = lattice_symmetry_operations(fcc())
        points = np.array([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        result = reduce_points_with_orbits(points, operations)
        self.assertEqual(result.n_orbits, 1)
        self.assertEqual(result.multiplicities.tolist(), [3])
        self.assertEqual(result.labels.tolist(), [0, 0, 0])

    def test_orbit(self):
        operations = lattice_symmetry_operations(simple_cubic())
        orbit = symmetry_orbit([1.0, 2.0, 3.0], operations)
        self.assertEqual(orbit.shape, (48, 3))
        self.assertEqual(len(np.unique(orbit, axis=0)), 48)


if __name__ == "__main__":
    unittest.main()
