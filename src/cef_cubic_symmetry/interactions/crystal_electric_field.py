"""The module contains CEF class."""


from common import physics, utils
from core.cef_parameters import BParameters
from interactions.base_interaction import BaseInteraction


class CrystalElectricField(BaseInteraction):
    """
    Class defining the trivalent rare earth compound,
    its crystal field parameters and the eigenvalues and eigenfunctions
    of the CEF Hamiltonian, if it is already diagonalized.

    """

    def __init__(self, parameters: BParameters, **kwargs) -> None:
        """Initializes the CEF object or read it from a file."""
        super().__init__(**kwargs)
        self.parameters = parameters

    def get_hamiltonian(self):
        """Determines the CEF Hamiltonian based on the input parameters."""
        momentum = self.sample.rare_earth.info.total_momentum_ground
        size = self.sample.rare_earth.matrix_size
        squared_momentum = self.sample.rare_earth.squared_momentum
        hamiltonian = utils.get_empty_matrix(size)
        for row in range(size):
            # row = 0...2J
            # mqn_1[1] = m = -J...J
            mqn_1 = [(row - momentum) ** i for i in range(5)]
            for key in ('20', '40', '60'):
                hamiltonian[row, row] += (
                        self.parameters.__getattribute__(f'b{key}') *
                        physics.steven_operators(
                            f'o{key}',
                            squared_momentum,
                            mqn_1,
                        )
                )
            for degree in range(2, size - row):
                mqn_2 = [(row - momentum + degree) ** i for i in range(5)]
                for key in ('22', '42', '62', '43', '63', '44', '64', '66'):
                    if key[-1] == str(degree):
                        hamiltonian[row, row + degree] += (
                                self.parameters.__getattribute__(f'b{key}') *
                                physics.steven_operators(
                                    f'o{key}',
                                    squared_momentum,
                                    mqn_1,
                                    mqn_2,
                                )
                        )
                hamiltonian[row + degree, row] = hamiltonian[row, row + degree]
        return hamiltonian

    def __str__(self):
        output = []
        for key, value in self.parameters.__dict__.items():
            if value:
                output.append(f'{key} = {value:.4f};')
        return '\n'.join(output)
