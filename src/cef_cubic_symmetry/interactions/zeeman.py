"""The module contains CEF class."""


from numpy import sqrt
from scipy.constants import physical_constants

from common import utils
from core.custom_datatypes import MagnetField
from interactions.base_interaction import BaseInteraction


class Zeeman(BaseInteraction):
    """
    Class defining the trivalent rare earth compound,
    its crystal field parameters and the eigenvalues and eigenfunctions
    of the CEF Hamiltonian, if it is already diagonalized.

    """

    def __init__(self, magnet_field: MagnetField = None, **kwargs) -> None:
        """Initializes the CEF object or read it from a file."""
        super().__init__(**kwargs)
        self.magnet_field = magnet_field or MagnetField()

    def get_hamiltonian(self):
        """Determines the Zeeman terms to the Hamiltonian."""
        size = self.sample.rare_earth.matrix_size
        hamiltonian = utils.get_empty_matrix(size)
        momentum = self.sample.rare_earth.info.total_momentum_ground
        squared_momentum = self.sample.rare_earth.squared_momentum
        factor = (
            float(self.sample.rare_earth.info.lande_factor)
            * physical_constants['Bohr magneton in eV/T'][0] * 1000
        )
        for row in range(size):
            # mqn_1 =  m = -J...J
            mqn_1 = row - momentum
            hamiltonian[row, row] -= factor * mqn_1 * self.magnet_field.z_
            if row < (size - 1):
                column = row + 1
                mqn_2 = mqn_1 + 1
                hamiltonian[row, column] -= (
                    0.5 * factor * self.magnet_field.x_
                    * sqrt((squared_momentum - mqn_1 * mqn_2))
                )
                hamiltonian[column, row] = hamiltonian[row, column]
        return hamiltonian

    def __str__(self):
        output = []
        for key, value in self.magnet_field.__dict__.items():
            if value:
                output.append(f'H{key[0]} = {value:.4f};')
        return '\n'.join(output)
