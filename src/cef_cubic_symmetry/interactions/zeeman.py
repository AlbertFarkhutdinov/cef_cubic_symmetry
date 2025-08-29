"""The module contains CEF class."""
import numpy as np
from numpy import sqrt
from scipy.constants import physical_constants

from cef_cubic_symmetry.auxiliary import utils
from cef_cubic_symmetry.core.custom_datatypes import MagnetField
from cef_cubic_symmetry.interactions.base_interaction import BaseInteraction


class Zeeman(BaseInteraction):
    """Class for interaction with Zeeman field."""

    def __init__(self, magnet_field: MagnetField = None, **kwargs) -> None:
        """Initialize the CEF object or read it from a file."""
        super().__init__(**kwargs)
        self.magnet_field = magnet_field or MagnetField()

    def get_hamiltonian(self) -> np.ndarray:
        """Determine the Zeeman terms to the Hamiltonian."""
        size = self.sample.rare_earth.matrix_size
        hamiltonian = utils.get_empty_matrix(size)
        momentum = self.sample.rare_earth.info.total_momentum_ground
        squared_momentum = self.sample.rare_earth.squared_momentum
        factor = (
            float(self.sample.rare_earth.info.lande_factor)
            * physical_constants['Bohr magneton in eV/T'][0] * 1000
        )
        for row in range(size):
            # mqn1 =  m = -J...J
            mqn1 = row - momentum
            hamiltonian[row, row] -= factor * mqn1 * self.magnet_field.z_
            if row < (size - 1):
                column = row + 1
                mqn2 = mqn1 + 1
                hamiltonian[row, column] -= (
                    0.5 * factor * self.magnet_field.x_
                    * sqrt(squared_momentum - mqn1 * mqn2)
                )
                hamiltonian[column, row] = hamiltonian[row, column]
        return hamiltonian

    def __str__(self) -> str:
        output = []
        for key, value in self.magnet_field.__dict__.items():
            if value:
                output.append(f'H{key[0]} = {value:.4f};')
        return '\n'.join(output)
