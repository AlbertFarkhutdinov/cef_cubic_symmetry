"""The module contains CEF class."""


from json import dump, load

from numpy import linspace, sqrt
from scipy.linalg import eigh

from common import utils, physics
from common.tabular_information import BOHR_MAGNETON
from common.path_utils import get_paths
from common.utils import OpenedFile, get_repr
from common.constants import Material


class CEF:
    """
    Class defining the trivalent rare earth compound,
    its crystal field parameters and the eigenvalues and eigenfunctions
    of the CEF Hamiltonian, if it is already diagonalized.

    """

    resolution = 1e-2
    threshold = 1e-4

    def __init__(self, material: Material):
        """Initializes the CEF object or read it from a file."""
        self.material = material
        self.file_name = get_paths(
            data_name='parameters',
            format_name='.json',
            material=self.material,
        )
        self.magnet_field = {'z': 0, 'x': 0}
        self.temperature = 0

    @property
    def parameters(self):
        """CEF parameters"""
        return {
            param: 0 for param in (
                'B20', 'B40', 'B60',
                'B22', 'B42', 'B62',
                'B43', 'B63',
                'B44', 'B64',
                'B66',
            )
        }

    def load_data(self):
        """Loads CEF object from file"""
        with OpenedFile(self.file_name) as file:
            self.__dict__.update(load(file))

    def save_to_file(self):
        """Saves parameters of the current object to file."""
        saved_object = {
            'crystal': self.material.crystal,
            'rare_earth': self.material.rare_earth.name,
            'parameters': self.parameters,
            'magnet_field': self.magnet_field,
        }
        with OpenedFile(self.file_name, mode='w') as file:
            dump(saved_object, file, indent=4, sort_keys=True)

    def get_cef_hamiltonian(self,
                            size: int,
                            j: float,
                            squared_j: float):
        """Determines the CEF Hamiltonian based on the input parameters."""
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
                               magnet_field: dict = None):
        """Determines the Zeeman terms to the Hamiltonian."""
        if magnet_field is None:
            magnet_field = self.magnet_field
        hamiltonian = utils.get_empty_matrix(size)
        for row in range(size):
            # mqn_1 =  m = -J...J
            mqn_1 = row - j
            hamiltonian[row, row] -= (
                self.material.rare_earth.lande_factor *
                BOHR_MAGNETON *
                mqn_1 *
                magnet_field['z']
            )
            if row < (size - 1):
                column = row + 1
                mqn_2 = mqn_1 + 1
                hamiltonian[row, column] -= (
                    0.5 * self.material.rare_earth.lande_factor *
                    BOHR_MAGNETON *
                    sqrt((squared_j - mqn_1 * mqn_2)) *
                    magnet_field['x']
                )
                hamiltonian[column, row] = hamiltonian[row, column]
        return hamiltonian

    def get_total_hamiltonian(self, magnet_field: dict = None):
        """Returns the total Hamiltonian including CEF and Zeeman terms."""
        size = self.material.rare_earth.matrix_size
        j = self.material.rare_earth.total_momentum_ground
        squared_j = j * (j + 1)
        return (self.get_cef_hamiltonian(size, j, squared_j) +
                self.get_zeeman_hamiltonian(size, j, squared_j, magnet_field))

    def get_eigenvalues_and_eigenfunctions(
            self,
            total_hamiltonian=None,
            ground_state_is_zero=True,
    ):
        """
        Calculates eigenvalues and eigenfunctions of the total Hamiltonian.

        """
        if total_hamiltonian is None:
            total_hamiltonian = self.get_total_hamiltonian()
        eigenvalues, eigenfunctions = eigh(total_hamiltonian)
        if ground_state_is_zero:
            eigenvalues = eigenvalues - min(eigenvalues)
        return eigenvalues, eigenfunctions

    def get_transition_probabilities(
            self,
            eigenfunctions,
    ):
        """
        Determines matrix elements for dipole transitions
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
                    row
                ] = transition_probability[row, column]

        return j_ops, transition_probability

    def get_bolzmann_factor(
            self,
            size: int,
            eigenvalues,
            temperature=None,
    ):
        """Determines bolzmann_factor at specified temperature."""
        temperature = utils.get_default(temperature or self.temperature)
        thermal = physics.thermodynamics(temperature, eigenvalues)
        bolzmann_factor = utils.get_empty_matrix(size, dimension=1)
        if thermal['temperature'] <= 0:
            bolzmann_factor[0] = 1
        else:
            bolzmann_factor = thermal['bolzmann'] / sum(thermal['bolzmann'])
        return bolzmann_factor

    def get_all_peaks(
            self,
            temperature=None,
            magnet_field: dict = None,
    ):
        """
        Determines the peak energies and intensities
        from the total Hamiltonian.

        """
        size = self.material.rare_earth.matrix_size
        if magnet_field is None:
            magnet_field = self.magnet_field
        total_hamiltonian = self.get_total_hamiltonian(magnet_field)
        eigenvalues, eigenfunctions = self.get_eigenvalues_and_eigenfunctions(
            total_hamiltonian
        )
        bolzmann_factor = self.get_bolzmann_factor(
            size, eigenvalues, temperature
        )
        peaks = []
        _, transition_probabilities = self.get_transition_probabilities(
            eigenfunctions
        )
        for level_1 in range(size):
            for level_2 in range(size):
                intensity_of_transition = (
                    transition_probabilities[level_2, level_1] *
                    bolzmann_factor[level_1]
                )
                if intensity_of_transition > 0:
                    peaks.append({
                        'energy': eigenvalues[level_2] - eigenvalues[level_1],
                        'intensity': intensity_of_transition
                    })
        return peaks

    def get_peaks(self,
                  temperature=None,
                  magnet_field: dict = None):
        """Returns peaks for non-degenerate levels."""
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

    def get_energies(self, peaks=None):
        """Returns transition energies"""
        if peaks is None:
            peaks = self.get_peaks()
        return [peak[0] for peak in peaks]

    def get_intensities(self, peaks=None):
        """Returns transition intensities"""
        if peaks is None:
            peaks = self.get_peaks()
        return [peak[1] for peak in peaks]

    def get_spectrum(self,
                     energies=None,
                     temperature=None,
                     width_dict: dict = None,
                     magnet_field: dict = None):
        """Calculates the neutron scattering cross section."""
        temperature = utils.get_default(temperature, self.temperature)
        peaks = self.get_peaks(temperature, magnet_field)
        eigenvalues, _ = self.get_eigenvalues_and_eigenfunctions()

        if energies is None:
            # 501 numbers in range from -1.1*E_max to 1.1*E_max
            energies = linspace(
                -1.1 * eigenvalues[-1],
                1.1 * eigenvalues[-1],
                501
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

    def get_moments(self,
                    temperature=None,
                    eigenvalues=None,
                    eigenfunctions=None):
        """Calculates the magnetic moments of the CEF model."""
        if eigenvalues is None and eigenfunctions is None:
            eigenvalues, eigenfunctions = self.get_eigenvalues_and_eigenfunctions()
        j_ops, _ = self.get_transition_probabilities(eigenfunctions)
        temperature = utils.get_default(temperature, self.temperature)
        thermal = physics.thermodynamics(temperature, eigenvalues)
        if thermal['temperature'] > 0:
            j_average = {'z': 0, 'x': 0}
            statistic_sum = sum(thermal['bolzmann'])
            for index in range(eigenvalues.size):
                j_average['z'] += (j_ops['z'][index, index] *
                                   thermal['bolzmann'][index])
                j_average['x'] += (0.5 * (j_ops['+'][index, index] +
                                          j_ops['-'][index, index]) *
                                   thermal['bolzmann'][index])
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
                      eigenvalues[eigenvalues == 0].size)
            }
        magnetic_moment = {}
        for key, value in j_average.items():
            magnetic_moment[key] = (
                    self.material.rare_earth.lande_factor
                    * value
            )
            # magnetic moments are given in units of Bohr magneton
        return j_average, magnetic_moment

    def get_chi(self,
                temperature=None,
                eigenvalues=None,
                eigenfunctions=None):
        """Calculates the susceptibility at a specified temperature."""
        if eigenvalues is None and eigenfunctions is None:
            eigenvalues, eigenfunctions = self.get_eigenvalues_and_eigenfunctions()
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
                            * thermal['bolzmann'][row]
                    )
                    chi['curie']['x'] += (
                            0.25 * (j_ops_square['+'] + j_ops_square['-']) *
                            thermal['bolzmann'][row]
                    )
                else:
                    chi['van_vleck']['z'] += (2 * j_ops_square['z'] *
                                              thermal['bolzmann'][row] /
                                              (column_value - row_value))
                    chi['van_vleck']['x'] += (0.5 * (j_ops_square['+'] +
                                                     j_ops_square['-']) *
                                              thermal['bolzmann'][row] /
                                              (column_value - row_value))
        coefficient = self.material.rare_earth.lande_factor ** 2
        if thermal['temperature'] > 0:
            coefficient = coefficient / sum(thermal['bolzmann'])
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
    ):
        """
        Calculates the susceptibility at a specified range of temperatures.

        """
        temperatures = utils.get_default(
            temperatures,
            linspace(1, 300, 300, dtype='float64')
        )
        if eigenvalues is None and eigenfunctions is None:
            eigenvalues, eigenfunctions = self.get_eigenvalues_and_eigenfunctions()
        temperatures = utils.get_default(
            temperatures,
            linspace(1, 300, 300, dtype='float64')
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
        for _, temperature in enumerate(temperatures):
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

    def __repr__(self):
        """Method returns string representation of the CEF object."""
        return get_repr(self, 'material')

    def __str__(self):
        """Return a summary of the model parameters.
        This includes the rare earth, the CEF parameters, and,
        if diagonalized, the eigenvalues and eigenfunctions."""
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
        for key, value in self.magnet_field.items():
            if value:
                output.append(f'H{key} = {value:.4f};')
            if magnetic_moment[key] > 1e-9:
                line_to_append = (
                        f'<J{key}> = {j_average[key]:7.3f}; ' +
                        f'<mu_{key}> = {magnetic_moment[key]:7.3f} mu_Bohr;'
                )
                output.append(line_to_append)

        eigenvalues, eigenfunctions = self.get_eigenvalues_and_eigenfunctions()
        if eigenvalues.any():
            output.append('Crystal Field Eigenvalues and Eigenfunctions:')
            for column in range(eigenvalues.size):
                line = [f'{eigenvalues[column]:8.3f}: ']
                for row in range(eigenvalues.size):
                    if abs(eigenfunctions[row, column]) > 0.0001:
                        j_z = (
                                row - self.material.rare_earth.total_momentum_ground
                        )
                        line_to_append = (
                                f'{utils.get_sign(eigenfunctions[row, column])}' +
                                f'{abs(eigenfunctions[row, column]):7.4f}' +
                                f'|{utils.get_sign(j_z)}{abs(j_z)}>')
                        line.append(line_to_append)
                output.append(' '.join(line))

        peaks = self.get_peaks()
        if peaks:
            output.append('Crystal Field Transitions:')
            output.append(f'Temperature: {self.temperature} K')
            for peak in peaks:
                output.append(
                    f'Energy: {peak[0]:8.3f} meV  Intensity: {peak[1]:8.4f}'
                )

        return '\n'.join(output)
