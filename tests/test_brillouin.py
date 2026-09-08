import unittest

import numpy as np

from brillouin_zone_toolkit import bcc, fcc, first_brillouin_zone, hexagonal, simple_cubic


class BrillouinZoneTests(unittest.TestCase):
    # (vertices, faces) of the Wigner--Seitz cells of reciprocal SC, BCC, FCC,
    # and hexagonal lattices, respectively.
    expected = {
        "sc": (8, 6),
        "fcc": (24, 14),
        "bcc": (14, 12),
        "hex": (12, 8),
    }

    def test_standard_zone_topology_and_volume(self):
        lattices = {"sc": simple_cubic(), "fcc": fcc(), "bcc": bcc(), "hex": hexagonal()}
        for name, lattice in lattices.items():
            with self.subTest(lattice=name):
                zone = first_brillouin_zone(lattice, search=2)
                self.assertEqual((zone.n_vertices, zone.n_faces), self.expected[name])
                self.assertEqual(len(zone.edges), zone.n_vertices + zone.n_faces - 2)
                self.assertAlmostEqual(zone.volume, lattice.reciprocal_volume, places=8)
                self.assertTrue(zone.contains(np.zeros(3)))
                self.assertTrue(np.all(zone.contains(zone.vertices)))

    def test_search_one_is_sufficient_for_standard_cells(self):
        for lattice in (simple_cubic(), fcc(), bcc(), hexagonal()):
            with self.subTest(lattice=lattice.name):
                one = first_brillouin_zone(lattice, search=1)
                many = first_brillouin_zone(lattice, search=2)
                self.assertEqual((one.n_vertices, one.n_faces), (many.n_vertices, many.n_faces))
                self.assertAlmostEqual(one.volume, many.volume, places=8)

    def test_array_input_is_a_reciprocal_basis(self):
        lattice = simple_cubic(2.0)
        zone = first_brillouin_zone(lattice.reciprocal_vectors)
        self.assertAlmostEqual(zone.volume, lattice.reciprocal_volume, places=8)

    def test_contains_shape_and_tolerance_validation(self):
        zone = first_brillouin_zone(simple_cubic())
        with self.assertRaises(ValueError):
            zone.contains(np.zeros((2, 2)))
        with self.assertRaises(ValueError):
            first_brillouin_zone(simple_cubic(), search=0)
        with self.assertRaises(ValueError):
            first_brillouin_zone(simple_cubic(), tolerance=0)


if __name__ == "__main__":
    unittest.main()
