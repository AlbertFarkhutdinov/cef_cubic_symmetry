"""The module contains CEF class."""

import json

import numpy as np
from numpy import linspace, sqrt
from scipy.constants import physical_constants
from scipy.linalg import eigh

from cef_cubic_symmetry.auxiliary import physics, utils
from cef_cubic_symmetry.auxiliary.path_utils import get_paths
from cef_cubic_symmetry.auxiliary.utils import UTF8File, get_repr
from cef_cubic_symmetry.core.sample import Sample


class CEF:
    """
    Class for cristal electric field.

    Class defining the trivalent rare earth compound,
    its crystal field parameters and the eigenvalues and eigenfunctions
    of the CEF Hamiltonian, if it is already diagonalized.

    """

    resolution = 1e-2
    threshold = 1e-4

    def __init__(self, material: Sample) -> None:
        """Initialize the CEF object or read it from a file."""
        self.material = material
        self.file_name = get_paths(
            data_name='parameters',
            format_name='.json',
            material=self.material,
        )
        self.magnet_field = {'z': 0, 'x': 0}
        self.temperature = 0

    @property
    def parameters(self) -> dict[str, float]:
        """CEF parameters."""
        return dict.fromkeys(
            (
                'B20',
                'B40',
                'B60',
                'B22',
                'B42',
                'B62',
                'B43',
                'B63',
                'B44',
                'B64',
                'B66',
            ),
            0,
        )

    def load_data(self) -> None:
        """Load CEF object from file."""
        with UTF8File(self.file_name) as file:
            self.__dict__.update(json.load(file))

    def save_to_file(self) -> None:
        """Save parameters of the current object to file."""
        saved_object = {
            'crystal': self.material.crystal,
            'rare_earth': self.material.rare_earth.name,
            'parameters': self.parameters,
            'magnet_field': self.magnet_field,
        }
        with UTF8File(self.file_name, mode='w') as file:
            json.dump(saved_object, file, indent=4, sort_keys=True)

    def get_cef_hamiltonian(
        self,
        size: int,
        j: float,
        squared_j: float,
    ) -> np.ndarray:
        """Determine the CEF Hamiltonian based on the input parameters."""
        hamiltonian = utils.get_empty_matrix(size)
        parameters = self.parameters
        for row in range(size):
            # row = 0...2J
            # mqn_1[1] = m = -J...J
            mqn_1 = [(row - j) ** i for i in range(5)]
            for key in ('20', '40', '60'):
                hamiltonian[row, row] += (
                    parameters[f'B{key}'] *
                    physics.steven_operators(
                        f'o{key}',
                        squared_j,
                        mqn_1,
                    )
                )
            for degree in range(2, size - row):
                mqn_2 = [(row - j + degree) ** i for i in range(5)]
                for key in ('22', '42', '62', '43', '63', '44', '64', '66'):
                    if key[-1] == str(degree):
                        hamiltonian[row, row + degree] += (
                            parameters[f'B{key}'] *
                            physics.steven_operators(
                                f'o{key}',
                                squared_j,
                                mqn_1,
                                mqn_2,
                            )
                        )
                hamiltonian[row + degree, row] = hamiltonian[row, row + degree]
        return hamiltonian

    def get_zeeman_hamiltonian(self,
                               size: int,
                               j: float,
                               squared_j: float,
                               magnet_field: dict | None = None) -> np.ndarray:
        """Determine the Zeeman terms to the Hamiltonian."""
        if magnet_field is None:
            magnet_field = self.magnet_field
        hamiltonian = utils.get_empty_matrix(size)
        for row in range(size):
            # mqn_1 =  m = -J...J
            mqn_1 = row - j
            hamiltonian[row, row] -= (
                self.material.rare_earth.lande_factor *
                physical_constants['Bohr magneton in eV/T'][0] * 1000 *
                mqn_1 *
                magnet_field['z']
            )
            if row < (size - 1):
                column = row + 1
                mqn_2 = mqn_1 + 1
                hamiltonian[row, column] -= (
                    0.5 * self.material.rare_earth.lande_factor *
                    physical_constants['Bohr magneton in eV/T'][0] * 1000 *
                    sqrt(squared_j - mqn_1 * mqn_2) *
                    magnet_field['x']
                )
                hamiltonian[column, row] = hamiltonian[row, column]
        return hamiltonian

    def get_total_hamiltonian(
        self,
        magnet_field: dict | None = None,
    ) -> np.ndarray:
        """Return the total Hamiltonian including CEF and Zeeman terms."""
        size = self.material.rare_earth.matrix_size
        j = self.material.rare_earth.total_momentum_ground
        squared_j = j * (j + 1)
        return (
            self.get_cef_hamiltonian(size, j, squared_j)
            + self.get_zeeman_hamiltonian(size, j, squared_j, magnet_field)
        )

    def get_eigenvalues_and_eigenfunctions(
        self,
        total_hamiltonian: np.ndarray = None,
        *,
        ground_state_is_zero: bool = True,
    ) -> tuple:
        """Return eigenvalues and eigenfunctions of the total Hamiltonian."""
        if total_hamiltonian is None:
            total_hamiltonian = self.get_total_hamiltonian()
        eigenvalues, eigenfunctions = eigh(total_hamiltonian)
        if ground_state_is_zero:
            eigenvalues -= min(eigenvalues)
        return eigenvalues, eigenfunctions

    def get_transition_probabilities(
        self,
        eigenfunctions,
    ) -> tuple:
        """
        Return transition probabilities.

        Determine matrix elements for dipole transitions
        between eigenfunctions of the total Hamiltonian.

        """
        j = self.material.rare_earth.total_momentum_ground
        squared_j = j * (j + 1)
        size = int(2 * j + 1)
        j_ops = {
            'z': utils.get_empty_matrix(size),
            '+': utils.get_empty_matrix(size),
            '-': utils.get_empty_matrix(size),
        }
        transition_probability = utils.get_empty_matrix(size)
        for row in range(size):
            j_ops['z'][row, row] += (
                eigenfunctions[size - 1, row] ** 2 *
                (size - 1 - j)
            )
            for row_j in range(size - 1):
                j_ops['z'][row, row] += (
                    eigenfunctions[row_j, row] ** 2 *
                    (row_j - j)
                )
                j_ops['+'][row, row] += (
                    eigenfunctions[row_j + 1, row] *
                    eigenfunctions[row_j, row] *
                    sqrt(squared_j - (row_j - j) * (row_j - j + 1))
                )

            j_ops['-'][row, row] = j_ops['+'][row, row]
            for column in range(row + 1, size):
                mqn_1 = size - 1 - j
                j_ops['z'][row, column] += (
                    eigenfunctions[size - 1, row] *
                    eigenfunctions[size - 1, column] * mqn_1
                )
                for row_j in range(size - 1):
                    mqn_1 = row_j - j
                    j_ops['z'][row, column] += (
                        eigenfunctions[row_j, row] *
                        eigenfunctions[row_j, column] * mqn_1
                    )
                    column_j = row_j + 1
                    mqn_2 = column_j - j
                    common_root = sqrt(squared_j - mqn_1 * mqn_2)
                    j_ops['+'][row, column] += (
                        eigenfunctions[column_j, row] *
                        eigenfunctions[row_j, column] *
                        common_root
                    )
                    j_ops['-'][row, column] += (
                        eigenfunctions[row_j, row] *
                        eigenfunctions[column_j, column] *
                        common_root
                    )

                transition_probability[row, column] = (
                    (2 * j_ops['z'][row, column] ** 2 +
                     j_ops['+'][row, column] ** 2 +
                     j_ops['-'][row, column] ** 2) / 3
                )
                j_ops['z'][column, row] = j_ops['z'][row, column]
                j_ops['+'][column, row] = j_ops['-'][row, column]
                j_ops['-'][column, row] = j_ops['+'][row, column]
                transition_probability[
                    column,
                    row,
                ] = transition_probability[row, column]

        return j_ops, transition_probability

    def get_boltzmann_factor(
        self,
        size: int,
        eigenvalues,
        temperature=None,
    ) -> float:
        """Determine boltzmann_factor at specified temperature."""
        temperature = utils.get_default(temperature, self.temperature)
        thermal = physics.thermodynamics(temperature, eigenvalues)
        boltzmann_factor = utils.get_empty_matrix(size, dimension=1)
        if thermal['temperature'] <= 0:
            boltzmann_factor[0] = 1
        else:
            boltzmann_factor = thermal['boltzmann'] / sum(thermal['boltzmann'])
        return boltzmann_factor

    def get_all_peaks(
        self,
        temperature=None,
        magnet_field: dict | None = None,
    ) -> list:
        """Determine the peak properties from the total Hamiltonian."""
        size = self.material.rare_earth.matrix_size
        if magnet_field is None:
            magnet_field = self.magnet_field
        total_hamiltonian = self.get_total_hamiltonian(magnet_field)
        eigenvalues, eigenfunctions = self.get_eigenvalues_and_eigenfunctions(
            total_hamiltonian,
        )
        boltzmann_factor = self.get_boltzmann_factor(
            size, eigenvalues, temperature,
        )
        peaks = []
        _, transition_probabilities = self.get_transition_probabilities(
            eigenfunctions,
        )
        for level_1 in range(size):
            for level_2 in range(size):
                intensity_of_transition = (
                    transition_probabilities[level_2, level_1] *
                    boltzmann_factor[level_1]
                )
                if intensity_of_transition > 0:
                    peaks.append({
                        'energy': eigenvalues[level_2] - eigenvalues[level_1],
                        'intensity': intensity_of_transition,
                    })
        return peaks

    def get_peaks(
        self,
        temperature=None,
        magnet_field: dict | None = None,
    ) -> list:
        """Return peaks for non-degenerate levels."""
        result = []
        peaks = self.get_all_peaks(temperature, magnet_field)
        for peak in peaks:
            sum_peaks = peak['energy'] * peak['intensity']
            for other_peak in peaks:
                if (
                    peak['intensity'] > 0 and
                    other_peak is not peak and
                    (
                        abs(peak['energy'] - other_peak['energy'])
                        < self.__class__.resolution
                    )
                ):
                    peak['intensity'] += other_peak['intensity']
                    sum_peaks += other_peak['energy'] * other_peak['intensity']
                    other_peak['intensity'] = 0
            if peak['intensity'] > self.__class__.threshold:
                peak['energy'] = sum_peaks / peak['intensity']
                result.append((peak['energy'], peak['intensity']))
        result.sort()
        intensity_sum = 2 * (
            self.material.rare_earth.total_momentum_ground *
            (self.material.rare_earth.total_momentum_ground + 1)
        ) / 3
        intensities = [item[1] for item in result]
        if sum(intensities) != intensity_sum:
            result[0] = (result[0][0], intensity_sum - sum(intensities[1:]))

        return result

    def get_energies(self, peaks=None) -> list:
        """Return transition energies."""
        if peaks is None:
            peaks = self.get_peaks()
        return [peak[0] for peak in peaks]

    def get_intensities(self, peaks=None) -> list:
        """Return transition intensities."""
        if peaks is None:
            peaks = self.get_peaks()
        return [peak[1] for peak in peaks]

    def get_spectrum(
        self,
        energies=None,
        temperature=None,
        width_dict: dict | None = None,
        magnet_field: dict | None = None,
    ) -> np.ndarray:
        """Calculate the neutron scattering cross-section."""
        temperature = utils.get_default(temperature, self.temperature)
        peaks = self.get_peaks(temperature, magnet_field)
        eigenvalues, _ = self.get_eigenvalues_and_eigenfunctions()

        if energies is None:
            # 501 numbers in range from -1.1*E_max to 1.1*E_max
            energies = linspace(
                -1.1 * eigenvalues[-1],
                1.1 * eigenvalues[-1],
                501,
            )
        if width_dict is None:
            width_dict = {'sigma': 0.01 * (max(energies) - min(energies))}

        spectrum = utils.get_empty_matrix(energies.size, dimension=1)

        sigma = width_dict.get('sigma', None)
        gamma = width_dict.get('gamma', None)
        for peak in peaks:
            if sigma and not gamma:
                spectrum += peak[1] * physics.gaussian_normalized(
                    energies,
                    peak[0],
                    sigma,
                )
            elif gamma and not sigma:
                spectrum += peak[1] * physics.lorentzian_normalized(
                    energies,
                    peak[0],
                    gamma,
                )
            elif sigma and gamma:
                spectrum += peak[1] * physics.pseudo_voigt_normalized(
                    energies,
                    peak[0],
                    sigma,
                    gamma,
                )

        spectrum *= 72.65 * self.material.rare_earth.lande_factor ** 2

        return spectrum

    def get_moments(
        self,
        temperature=None,
        eigenvalues=None,
        eigenfunctions=None,
    ) -> tuple:
        """Calculate the magnetic moments of the CEF model."""
        if eigenvalues is None and eigenfunctions is None:
            eigenvalues, eigenfunctions = (
                self.get_eigenvalues_and_eigenfunctions()
            )
        j_ops, _ = self.get_transition_probabilities(eigenfunctions)
        temperature = utils.get_default(temperature, self.temperature)
        thermal = physics.thermodynamics(temperature, eigenvalues)
        if thermal['temperature'] > 0:
            j_average = {'z': 0, 'x': 0}
            statistic_sum = sum(thermal['boltzmann'])
            for index in range(eigenvalues.size):
                j_average['z'] += (j_ops['z'][index, index] *
                                   thermal['boltzmann'][index])
                j_average['x'] += (0.5 * (j_ops['+'][index, index] +
                                          j_ops['-'][index, index]) *
                                   thermal['boltzmann'][index])
            for key, value in j_average.items():
                j_average[key] = value / statistic_sum
        else:
            j_average = {
                'z': (sum(j_ops['z'][eigenvalues == 0, eigenvalues == 0]) /
                      eigenvalues[eigenvalues == 0].size),
                'x': (sum(0.5 * (j_ops['+'][eigenvalues == 0,
                                            eigenvalues == 0] +
                                 j_ops['-'][eigenvalues == 0,
                                            eigenvalues == 0])) /
                      eigenvalues[eigenvalues == 0].size),
            }
        magnetic_moment = {}
        for key, value in j_average.items():
            magnetic_moment[key] = (
                self.material.rare_earth.lande_factor
                * value
            )
            # magnetic moments are given in units of Bohr magneton
        return j_average, magnetic_moment

    def get_chi(
        self,
        temperature=None,
        eigenvalues=None,
        eigenfunctions=None,
    ) -> dict:
        """Calculate the susceptibility at a specified temperature."""
        if eigenvalues is None and eigenfunctions is None:
            eigenvalues, eigenfunctions = (
                self.get_eigenvalues_and_eigenfunctions()
            )
        j_ops, _ = self.get_transition_probabilities(eigenfunctions)
        thermal = physics.thermodynamics(utils.get_default(temperature,
                                                           self.temperature),
                                         eigenvalues)
        chi = {
            'curie': {'z': 0, 'x': 0},
            'van_vleck': {'z': 0, 'x': 0},
        }
        for row in range(eigenvalues.size):
            for column in range(eigenvalues.size):
                j_ops_square = {
                    key: val[row, column] ** 2
                    for key, val in j_ops.items()
                }
                row_value = eigenvalues[row]
                column_value = eigenvalues[column]
                if (
                    abs(column_value - row_value)
                    < 0.00001 * thermal['temperature']
                ):
                    chi['curie']['z'] += (
                        j_ops_square['z']
                        * thermal['boltzmann'][row]
                    )
                    chi['curie']['x'] += (
                        0.25 * (j_ops_square['+'] + j_ops_square['-']) *
                        thermal['boltzmann'][row]
                    )
                else:
                    chi['van_vleck']['z'] += (2 * j_ops_square['z'] *
                                              thermal['boltzmann'][row] /
                                              (column_value - row_value))
                    chi['van_vleck']['x'] += (0.5 * (j_ops_square['+'] +
                                                     j_ops_square['-']) *
                                              thermal['boltzmann'][row] /
                                              (column_value - row_value))
        coefficient = self.material.rare_earth.lande_factor ** 2
        if thermal['temperature'] > 0:
            coefficient /= sum(thermal['boltzmann'])
        for key in ('z', 'x'):
            chi['curie'][key] = (
                coefficient / thermal['temperature'] * chi['curie'][key]
            )
            chi['van_vleck'][key] = coefficient * chi['van_vleck'][key]
        return chi

    def get_chi_dependence(
        self,
        temperatures=None,
        eigenvalues=None,
        eigenfunctions=None,
    ) -> tuple:
        """Calculate the susceptibility at specified temperatures."""
        temperatures = utils.get_default(
            temperatures,
            linspace(1, 300, 300, dtype='float64'),
        )
        if eigenvalues is None and eigenfunctions is None:
            eigenvalues, eigenfunctions = (
                self.get_eigenvalues_and_eigenfunctions()
            )
        temperatures = utils.get_default(
            temperatures,
            linspace(1, 300, 300, dtype='float64'),
        )
        chi_curie = {
            'z': None,
            'x': None,
        }
        chi_van_vleck = {
            'z': None,
            'x': None,
        }
        chi = {
            'z': None,
            'x': None,
            'total': None,
            'inverse': None,
        }
        for key in ('z', 'x'):
            chi_curie[key] = utils.get_empty_matrix(temperatures.shape)
            chi_van_vleck[key] = utils.get_empty_matrix(temperatures.shape)
            chi[key] = utils.get_empty_matrix(temperatures.shape)

        chi = {
            'total': utils.get_empty_matrix(temperatures.shape),
            'inverse': utils.get_empty_matrix(temperatures.shape),
        }
        for temperature in temperatures:
            current_chi = self.get_chi(
                temperature,
                eigenvalues,
                eigenfunctions,
            )
            for key in ('z', 'x'):
                chi_curie[key] = current_chi['curie'][key]
                chi_van_vleck[key] = current_chi['van_vleck'][key]
                chi[key] = chi_curie[key] + chi_van_vleck[key]
            chi['total'] = (chi['z'] + 2 * chi['x']) / 3
            chi['inverse'] = 1 / chi['total']

        return chi_curie, chi_van_vleck, chi

    def __repr__(self) -> str:
        """Return string representation of the CEF object."""
        return get_repr(self, 'material')

    def __str__(self) -> str:
        """
        Return a summary of the model parameters.

        This includes the rare earth, the CEF parameters, and,
        if diagonalized, the eigenvalues and eigenfunctions.

        """
        output = [
            self.material.crystal,
            f'Rare-earth ion: {self.material.rare_earth.name};',
            f'Number of 4f-electrons = '
            f'{self.material.rare_earth.number_of_f_electrons};',
            f'J = {self.material.rare_earth.total_momentum_ground};',
        ]

        j_average, magnetic_moment = self.get_moments()

        for key, value in self.parameters.items():
            if value:
                output.append(f'{key} = {value:.4f};')
        threshold = 1e-9
        for key, value in self.magnet_field.items():
            if value:
                output.append(f'H{key} = {value:.4f};')
            if magnetic_moment[key] > threshold:
                line_to_append = '; '.join(
                    [
                        f'<J{key}> = {j_average[key]:7.3f}',
                        f'<mu_{key}> = {magnetic_moment[key]:7.3f} mu_Bohr',
                    ],
                )
                output.append(line_to_append)

        eigenvalues, eigenfunctions = self.get_eigenvalues_and_eigenfunctions()
        if eigenvalues.any():
            output.append('Crystal Field Eigenvalues and Eigenfunctions:')
            threshold = 1e-4
            for column in range(eigenvalues.size):
                line = [f'{eigenvalues[column]:8.3f}: ']
                for row in range(eigenvalues.size):
                    if abs(eigenfunctions[row, column]) > threshold:
                        tmg = self.material.rare_earth.total_momentum_ground
                        j_z = row - tmg
                        eigen_function = eigenfunctions[row, column]
                        line_to_append = ''.join(
                            [
                                f'{utils.get_sign(eigen_function)}',
                                f'{abs(eigen_function):7.4f}',
                                f'|{utils.get_sign(j_z)}{abs(j_z)}>',
                            ],
                        )
                        line.append(line_to_append)
                output.append(' '.join(line))

        peaks = self.get_peaks()
        if peaks:
            output.extend(
                (
                    'Crystal Field Transitions:',
                    f'Temperature: {self.temperature} K',
                ),
            )
            output.extend(
                f'Energy: {peak[0]:8.3f} meV  Intensity: {peak[1]:8.4f}'
                for peak in peaks
            )

        return '\n'.join(output)
