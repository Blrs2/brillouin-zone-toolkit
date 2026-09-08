import unittest

import numpy as np

from brillouin_zone_toolkit import (
    bcc,
    direct_reciprocal_duality,
    fcc,
    hexagonal,
    make_lattice,
    reciprocal_cell_volume,
    reciprocal_vectors,
    simple_cubic,
)


class LatticeTests(unittest.TestCase):
    def test_reciprocal_duality(self):
        for lattice in (simple_cubic(), fcc(), bcc(), hexagonal()):
            np.testing.assert_allclose(
                direct_reciprocal_duality(lattice.direct_vectors),
                2.0 * np.pi * np.eye(3),
                rtol=1e-12,
                atol=1e-12,
            )
            np.testing.assert_allclose(
                reciprocal_cell_volume(lattice.direct_vectors), lattice.reciprocal_volume
            )

    def test_known_primitive_cell_volumes(self):
        self.assertAlmostEqual(simple_cubic(2).volume, 8.0)
        self.assertAlmostEqual(fcc(2).volume, 2.0)
        self.assertAlmostEqual(bcc(2).volume, 4.0)
        self.assertAlmostEqual(hexagonal(2, 3).volume, 2.0 * np.sqrt(3.0) * 3.0)

    def test_factory_aliases(self):
        self.assertEqual(make_lattice("SC"), simple_cubic())
        self.assertEqual(make_lattice("face-centered-cubic"), fcc())
        self.assertEqual(make_lattice("body_centered_cubic"), bcc())
        self.assertEqual(make_lattice("hexagonal"), hexagonal())

    def test_reciprocal_vectors_are_not_a_copy_of_input(self):
        direct = np.eye(3)
        lattice = simple_cubic()
        self.assertTrue(direct.flags.writeable)
        self.assertFalse(lattice.vectors.flags.writeable)
        direct[0, 0] = 3.0
        self.assertAlmostEqual(lattice.vectors[0, 0], 1.0)
        np.testing.assert_allclose(reciprocal_vectors(np.eye(3)), 2.0 * np.pi * np.eye(3))


if __name__ == "__main__":
    unittest.main()
