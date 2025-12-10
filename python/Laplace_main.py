import numpy as np
from GP_main import TetrahedronPotentialCalculatorOOP

class TetrahedronLaplaceTest:
    def __init__(self, A, B, C, D, G=1, sigma=1):
        self.A, self.B, self.C, self.D = A, B, C, D
        self.G, self.sigma = G, sigma

    def laplacian(self, func, x, y, z, h=1e-5):
        """Numerical Laplacian using central differences."""
        return (
            (func(x + h, y, z) - 2*func(x, y, z) + func(x - h, y, z)) / h**2 +
            (func(x, y + h, z) - 2*func(x, y, z) + func(x, y - h, z)) / h**2 +
            (func(x, y, z + h) - 2*func(x, y, z) + func(x, y, z - h)) / h**2
        )

    def compute_laplacian_grid(self, z_plane, x_range, y_range, resolution=100):
        """Compute Laplacian on a grid at z = z_plane."""
        x_vals = np.linspace(*x_range, resolution)
        y_vals = np.linspace(*y_range, resolution)
        X, Y = np.meshgrid(x_vals, y_vals)
        Z = np.full_like(X, z_plane)

        # Define φ using the OOP potential
        def phi(x, y, z):
            calc = TetrahedronPotentialCalculatorOOP(self.A, self.B, self.C, self.D, np.array([x, y, z]), self.G, self.sigma)
            return calc.compute()

        lap_grid = np.zeros_like(X)
        for i in range(resolution):
            for j in range(resolution):
                lap_grid[i, j] = self.laplacian(phi, X[i, j], Y[i, j], Z[i, j])

        return X, Y, lap_grid
