"""The module contains CEF class."""


import numpy as np
from base_interaction import BaseInteraction

from common import physics, utils


class Thermostat(BaseInteraction):
    """
    Class defining the trivalent rare earth compound,
    its crystal field parameters and the eigenvalues and eigenfunctions
    of the CEF Hamiltonian, if it is already diagonalized.

    """

    def __init__(self, temperature: float = 0, **kwargs) -> None:
        """Initializes the CEF object or read it from a file."""
        super().__init__(**kwargs)
        self.temperature = temperature

    def get_boltzmann_factor(self, eigenvalues: np.ndarray) -> np.ndarray:
        """Determines boltzmann_factor at specified temperature."""
        thermal = physics.thermodynamics(self.temperature, eigenvalues)
        boltzmann_factors = utils.get_empty_matrix(
            self.sample.rare_earth.matrix_size,
            dimension=1,
        )
        if thermal['temperature'] <= 0:
            boltzmann_factors[0] = 1
        else:
            boltzmann_factors = (
                    thermal['boltzmann'] / sum(thermal['boltzmann'])
            )
        return boltzmann_factors

    def __str__(self):
        return f'Temperature: {self.temperature} K'
