from abc import ABC, abstractmethod
from typing import Optional

import numpy as np

# constants
_ANG2BOHR = 1.8897259886

# default settings
_DEFAULT_DELTA_X = 0.01  # Angstrom


class BasePES(ABC):
    """Abstract PES class"""

    @staticmethod
    def _check_coords(
        num_atoms: float,
        coords: np.array,
    ) -> None:
        if coords.shape != (num_atoms, 3):
            raise ValueError(
                f"The shape of coords is {coords.shape}, which is not "
                f"equal to ({num_atoms}, 3). Please check your input."
            )

    @abstractmethod
    def calc_energy(
        self,
        coords: np.array,
    ) -> float:
        """Calculate enenrgy at given coordinates.

        Args:
            coords: Numpy array of coordinates in Angstrom, shape (N, 3)

        Returns:
            energy: Energy in Hartree
        """

        pass

    def calc_gradients(
        self,
        coords: np.array,
        delta_x: Optional[float] = _DEFAULT_DELTA_X,
    ) -> np.array:
        """Calculcate gradients at given coordinates.

        Args:
            coords: Numpy array of coordinates in Angstrom, shape (N, 3)
            delta_x: step size of finite difference, default 0.01 Angstrom,
                recommanded by Gaussian

        Returns:
            gradients: Numpy array of gradients in Hartree / Bohr, shape (N, 3)
        """

        num_atoms = len(coords)
        gradients = np.zeros_like(coords)

        for i in range(num_atoms):
            for j in range(3):
                coords[i, j] += delta_x
                energy1 = self.calc_energy(coords)
                coords[i, j] -= 2 * delta_x
                energy2 = self.calc_energy(coords)
                coords[i, j] += delta_x
                gradients[i, j] = (energy1 - energy2) / \
                    (2 * delta_x * _ANG2BOHR)

        return gradients

    def calc_force_constants(
        self,
        coords: np.array,
        delta_x: Optional[float] = _DEFAULT_DELTA_X,
    ) -> np.array:
        """Calculate force constants at given coordinates

        Args:
            coords: Numpy array of coordinates in Angstrom
            delta_x: step size of finite difference, default 0.01 Angstrom
                recommaded by Gaussian

        Returns:
            force_constants: Numpy array of lower triangular of Hessian, in the
                unit of Hartree^2 / Bohr^2, shape (3 * N * (3 * N + 1)) / 2
        """

        def _calc_gradient_component(X, i):
            return self.calc_gradients(
                X.reshape((-1, 3)),
                delta_x=delta_x,
            ).reshape((-1,))[i]

        X = coords.reshape((-1, 1))
        force_constants = []

        for i in range(len(X)):
            for j in range(i + 1):
                # forward
                X[i] += delta_x
                g1 = _calc_gradient_component(X, j)
                X[i] -= delta_x

                # backward
                X[i] -= delta_x
                g2 = _calc_gradient_component(X, j)
                X[i] += delta_x

                res = (g1 - g2) / (2.0 * delta_x * _ANG2BOHR)
                force_constants.append(res)

        return np.array(force_constants)
