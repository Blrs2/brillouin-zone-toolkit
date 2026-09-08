import unittest

from brillouin_zone_toolkit import fcc, first_brillouin_zone


class PlottingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            import matplotlib

            matplotlib.use("Agg")
            import matplotlib.pyplot as plt

            cls.plt = plt
        except ImportError:
            cls.plt = None

    def test_plots_can_be_rendered_headlessly(self):
        if self.plt is None:
            self.skipTest("matplotlib is not installed")
        from brillouin_zone_toolkit.plotting import (
            plot_brillouin_zone,
            plot_reciprocal_lattice,
            plot_symmetry_reduction,
            sample_cartesian_grid,
        )
        from brillouin_zone_toolkit.symmetry import lattice_symmetry_operations

        lattice = fcc()
        zone = first_brillouin_zone(lattice)
        points = sample_cartesian_grid(zone, points_per_axis=5)
        operations = lattice_symmetry_operations(lattice)
        figures = []
        for plot in (
            lambda ax: plot_brillouin_zone(zone, ax=ax),
            lambda ax: plot_reciprocal_lattice(lattice.reciprocal_vectors, ax=ax),
            lambda ax: plot_symmetry_reduction(zone, points, operations=operations, ax=ax),
        ):
            figure = self.plt.figure()
            axis = figure.add_subplot(111, projection="3d")
            self.assertIs(plot(axis), axis)
            figures.append(figure)
        for figure in figures:
            figure.canvas.draw()
            self.plt.close(figure)


if __name__ == "__main__":
    unittest.main()
