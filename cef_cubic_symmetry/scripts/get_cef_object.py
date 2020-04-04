"""The module contains CF class."""
from json import dump, load
from numpy import linspace, sqrt
from scripts.common import utils
from scripts.common.physics import thermodynamics, gauss, lorentz, pseudo_voigt, steven_operators
from scripts.common.tabular_information import BOHR_MAGNETON, RARE_EARTHS, F4
from scripts.common.path_utils import get_paths
from scripts.common.utils import OpenedFile
from scripts.common.constants import PATH_TO_SAVED_OBJECTS
from scripts.common.named_tuples import Material
from scipy.linalg import eigh

RARE_EARTHS_NAMES = [element.name for element in RARE_EARTHS]


class CF:
    """Class defining the trivalent rare earth compound,
    its crystal field parameters and the eigenvalues and eigenfunctions
    of the CF Hamiltonian, if it is already diagonalized."""

    resolution = 0.01
    threshold = 0.0001

    def __init__(self, material: Material, file_name=None):
        """Initializes the CF object or read it from a file."""
        if not file_name:
            self.rare_earth = material.rare_earth
            self.crystal = material.crystal
            self.parameters = {
                'B20': 0, 'B40': 0, 'B60': 0,
                'B22': 0, 'B42': 0, 'B62': 0,
                'B43': 0, 'B63': 0,
                'B44': 0, 'B64': 0,
                'B66': 0,
            }
            self.magnet_field = {'z': 0, 'x': 0}
        else:
            with OpenedFile(file_name) as file:
                self.__dict__.update(load(file))

        if 'rare_earth' in self.__dict__:
            self.number_of_f_electrons = RARE_EARTHS_NAMES.index(self.rare_earth) + 1
            self.rare_earth = RARE_EARTHS[self.number_of_f_electrons - 1]
        else:
            self.number_of_f_electrons = 1
            self.rare_earth = RARE_EARTHS[0]
        self.temperature = 0
        self.temperature_used = self.temperature

    def save_to_file(self):
        """Saves parameters of the current object to file."""
        saved_object = {
            'crystal': self.crystal,
            'rare_earth': self.rare_earth.name,
            'number_of_f_electrons': RARE_EARTHS.index(self.rare_earth) + 1,
            'parameters': self.parameters,
            'magnet_field': self.magnet_field,
        }
        file_name = get_paths(directory=PATH_TO_SAVED_OBJECTS,
                              data_name='parameters',
                              format_name='.json',
                              material=Material(crystal=self.crystal,
                                                rare_earth=self.rare_earth.name))
        with OpenedFile(file_name, mode='w') as file:
            dump(saved_object, file, indent=4, sort_keys=True)

    def get_llw_cubic_object(self, parameters: dict):
        """Returns CF object for crystals with cubic symmetry"""
        if self.rare_earth.name in ['Ce', 'Sm', 'Eu']:
            print('This element is not supported. The default object will be returned.')
        else:
            self.parameters['B40'] = parameters['w'] * parameters['x'] / F4
            self.parameters['B44'] = 5 * self.parameters['B40']
            self.parameters['B60'] = (parameters['w'] * (1 - abs(parameters['x'])) /
                                      self.rare_earth.f_6)
            self.parameters['B64'] = -21 * self.parameters['B60']
        return self

    def get_cef_hamiltonian(self):
        """Determines the CF Hamiltonian based on the input parameters."""
        size = self.rare_earth.matrix_size
        hamiltonian = utils.get_empty_matrix(size)
        j = self.rare_earth.total_momentum_ground
        squared_j = j * (j + 1)

        for row in range(size):
            # row = 0...2J
            # mqn_1[1] = m = -J...J
            mqn_1 = [(row - j) ** i for i in range(5)]
            for key in ('20', '40', '60'):
                hamiltonian[row, row] += (self.parameters[f'B{key}'] *
                                          steven_operators(f'o{key}', squared_j, mqn_1))
            for degree in range(2, size - row):
                mqn_2 = [(row - j + degree) ** i for i in range(5)]
                for key in ('22', '42', '62', '43', '63', '44', '64', '66'):
                    if key[-1] == str(degree):
                        hamiltonian[row, row + degree] += (self.parameters[f'B{key}'] *
                                                           steven_operators(f'o{key}',
                                                                            squared_j,
                                                                            mqn_1,
                                                                            mqn_2))
                hamiltonian[row + degree, row] = hamiltonian[row, row + degree]
        return hamiltonian

    def get_zeeman_hamiltonian(self):
        """Determines the Zeeman terms to the Hamiltonian."""
        size = self.rare_earth.matrix_size
        j = self.rare_earth.total_momentum_ground
        squared_j = j * (j + 1)
        hamiltonian = utils.get_empty_matrix(size)
        for row in range(size):
            # mqn_1 =  m = -J...J
            mqn_1 = row - j
            hamiltonian[row, row] -= (self.rare_earth.lande_factor *
                                      BOHR_MAGNETON *
                                      mqn_1 *
                                      self.magnet_field['z'])
            if row < (size - 1):
                column = row + 1
                mqn_2 = mqn_1 + 1
                hamiltonian[row, column] -= (0.5 * self.rare_earth.lande_factor *
                                             BOHR_MAGNETON *
                                             sqrt((squared_j - mqn_1 * mqn_2)) *
                                             self.magnet_field['x'])
                hamiltonian[column, row] = hamiltonian[row, column]
        return hamiltonian

    def get_total_hamiltonian(self):
        """Returns the total Hamiltonian including CF and Zeeman terms."""
        return self.get_cef_hamiltonian() + self.get_zeeman_hamiltonian()

    def get_eigenvalues_and_eigenfunctions(self):
        """Calculates eigenvalues and eigenfunctions of the total Hamiltonian."""
        eigenvalues, eigenfunctions = eigh(self.get_total_hamiltonian())
        # E = 0 - minimum of energy.
        eigenvalues = eigenvalues - min(eigenvalues)
        return eigenvalues, eigenfunctions

    def get_transition_probabilities(self):
        """Determines matrix elements for dipole transitions
        between eigenfunctions of the total Hamiltonian."""
        j = self.rare_earth.total_momentum_ground
        squared_j = j * (j + 1)
        size = self.rare_earth.matrix_size
        _, eigenfunctions = self.get_eigenvalues_and_eigenfunctions()
        j_z = utils.get_empty_matrix(size)
        j_plus = utils.get_empty_matrix(size)
        j_minus = utils.get_empty_matrix(size)
        transition_probability = utils.get_empty_matrix(size)

        for row in range(size):
            j_z[row, row] += (eigenfunctions[size - 1, row] ** 2 * (size - 1 - j))
            for row_j in range(size - 1):
                j_z[row, row] += (eigenfunctions[row_j, row] ** 2) * (row_j - j)
                j_plus[row, row] += (eigenfunctions[row_j + 1, row] *
                                     eigenfunctions[row_j, row] *
                                     sqrt(squared_j - (row_j - j) * (row_j - j + 1)))

            j_minus[row, row] = j_plus[row, row]
            for column in range(row + 1, size):
                mqn_1 = size - 1 - j
                j_z[row, column] += (eigenfunctions[size - 1, row] *
                                     eigenfunctions[size - 1, column] * mqn_1)
                for row_j in range(size - 1):
                    mqn_1 = row_j - j
                    j_z[row, column] += (eigenfunctions[row_j, row] *
                                         eigenfunctions[row_j, column] * mqn_1)
                    column_j = row_j + 1
                    mqn_2 = column_j - j
                    common_root = sqrt(squared_j - mqn_1 * mqn_2)
                    j_plus[row, column] += (eigenfunctions[column_j, row] *
                                            eigenfunctions[row_j, column] *
                                            common_root)
                    j_minus[row, column] += (eigenfunctions[row_j, row] *
                                             eigenfunctions[column_j, column] *
                                             common_root)

                transition_probability[row, column] = ((2 * j_z[row, column] ** 2 +
                                                        j_plus[row, column] ** 2 +
                                                        j_minus[row, column] ** 2) / 3)
                j_z[column, row] = j_z[row, column]
                j_plus[column, row] = j_minus[row, column]
                j_minus[column, row] = j_plus[row, column]
                transition_probability[column, row] = transition_probability[row, column]

        return j_z, j_plus, j_minus, transition_probability

    def get_peaks(self, temperature=None, magnet_field: dict = None):
        """Determines the peak intensities from the total Hamiltonian."""
        size = self.rare_earth.matrix_size
        old_magnet_field = {
            'x': self.magnet_field['x'],
            'z': self.magnet_field['z'],
        }
        if magnet_field is None:
            magnet_field = {
                'x': 0,
                'z': 0,
            }
        for key, value in magnet_field.items():
            if value:
                self.magnet_field[key] = value
        eigenvalues, _ = self.get_eigenvalues_and_eigenfunctions()
        _, _, _, transition_probabilities = self.get_transition_probabilities()
        bolzmann_factor = utils.get_empty_matrix(size, dimension=1)
        temperature = utils.get_new_if_old_is_none(temperature, self.temperature)
        thermal = thermodynamics(temperature, eigenvalues)
        if thermal['temperature'] <= 0:
            bolzmann_factor[0] = 1
        else:
            statistic_sum = sum(thermal['bolzmann'])
            bolzmann_factor = thermal['bolzmann'] / statistic_sum

        peaks = []
        for level_1 in range(size):
            for level_2 in range(size):
                energy_of_transition = eigenvalues[level_2] - eigenvalues[level_1]
                intensity_of_transition = (transition_probabilities[level_2, level_1] *
                                           bolzmann_factor[level_1])
                if intensity_of_transition > 0:
                    peaks.append({'energy': energy_of_transition,
                                  'intensity': intensity_of_transition})

        result = []
        for peak in peaks:
            if peak['intensity'] > 0:
                sum_peaks = peak['energy'] * peak['intensity']
                for other_peak in peaks:
                    if other_peak is not peak:
                        if (other_peak['intensity'] > 0 and
                                abs(peak['energy'] - other_peak['energy'])
                                < self.__class__.resolution):
                            peak['intensity'] += other_peak['intensity']
                            sum_peaks += other_peak['energy'] * other_peak['intensity']
                            other_peak['intensity'] = 0

                if peak['intensity'] > self.__class__.threshold:
                    peak['energy'] = sum_peaks / peak['intensity']
                    result.append((peak['energy'], peak['intensity']))
        result.sort()
        self.temperature_used = temperature
        for key, value in magnet_field.items():
            if value:
                self.magnet_field[key] = old_magnet_field[key]

        return result

    def get_energies(self):
        """Returns transition energies"""
        return [peak[0] for peak in self.get_peaks()]

    def get_intensities(self):
        """Returns transition intensities"""
        return [peak[1] for peak in self.get_peaks()]

    def get_spectrum(self, energy=None, temperature=None, width_dict: dict = None,
                     magnet_field: dict = None):
        """Calculates the neutron scattering cross section."""
        temperature = utils.get_new_if_old_is_none(temperature, self.temperature)
        peaks = self.get_peaks(temperature, magnet_field)
        eigenvalues, _ = self.get_eigenvalues_and_eigenfunctions()

        if energy is None:
            # 501 numbers in range from -1.1*E_max to 1.1*E_max
            energy = linspace(-1.1 * eigenvalues[-1], 1.1 * eigenvalues[-1], 501)
        if width_dict is None:
            width_dict = {'sigma': 0.01 * (max(energy) - min(energy))}

        spectrum = utils.get_empty_matrix(energy.size, dimension=1)

        sigma = width_dict.get('sigma', None)
        gamma = width_dict.get('gamma', None)
        for peak in peaks:
            if sigma and not gamma:
                spectrum += peak[1] * gauss(energy, peak[0],
                                            sigma)
            elif gamma and not sigma:
                spectrum += peak[1] * lorentz(energy, peak[0],
                                              gamma)
            elif sigma and gamma:
                spectrum += peak[1] * pseudo_voigt(energy, peak[0],
                                                   sigma,
                                                   gamma)

        spectrum *= 72.65 * self.rare_earth.lande_factor ** 2

        return spectrum

    def get_moments(self, temperature=None):
        """Calculates the magnetic moments of the CF model."""
        eigenvalues, _ = self.get_eigenvalues_and_eigenfunctions()
        j_z, j_plus, j_minus, _ = self.get_transition_probabilities()
        temperature = utils.get_new_if_old_is_none(temperature, self.temperature)
        if temperature != self.temperature_used:
            self.get_peaks(temperature)
        thermal = thermodynamics(temperature, eigenvalues)
        if thermal['temperature'] > 0:
            j_average = {'z': 0, 'x': 0}
            statistic_sum = sum(thermal['bolzmann'])
            for index in range(eigenvalues.size):
                j_average['z'] += (j_z[index, index] *
                                   thermal['bolzmann'][index])
                j_average['x'] += (0.5 * (j_plus[index, index] +
                                          j_minus[index, index]) *
                                   thermal['bolzmann'][index])
            for key, value in j_average.items():
                j_average[key] = value / statistic_sum
        else:
            j_average = {
                'z': (sum(j_z[eigenvalues == 0, eigenvalues == 0]) /
                      eigenvalues[eigenvalues == 0].size),
                'x': (sum(0.5 * (j_plus[eigenvalues == 0,
                                        eigenvalues == 0] +
                                 j_minus[eigenvalues == 0,
                                         eigenvalues == 0])) /
                      eigenvalues[eigenvalues == 0].size)
            }
        magnetic_moment = {}
        for key, value in j_average.items():
            magnetic_moment[key] = self.rare_earth.lande_factor * value
            # magnetic moments are given in units of Bohr magneton
        return j_average, magnetic_moment

    def get_chi(self, temperature=None):
        """Calculates the susceptibility at a specified temperature."""
        eigenvalues, _ = self.get_eigenvalues_and_eigenfunctions()
        j_z, j_plus, j_minus, _ = self.get_transition_probabilities()
        temperature = utils.get_new_if_old_is_none(temperature, self.temperature)
        thermal = thermodynamics(temperature, eigenvalues)
        temperature_as_energy = thermal['temperature']
        statistic_sum = 1
        if temperature_as_energy > 0:
            statistic_sum = sum(thermal['bolzmann'])

        coefficient = self.rare_earth.lande_factor ** 2 / statistic_sum
        chi_curie = {'z': 0, 'x': 0}
        chi_van_vleck = {'z': 0, 'x': 0}
        for row in range(eigenvalues.size):
            for column in range(eigenvalues.size):
                current_bolzmann = thermal['bolzmann'][row]
                j_z_square = j_z[row, column] ** 2
                j_plus_square = j_plus[row, column] ** 2
                j_minus_square = j_minus['-'][row, column] ** 2
                row_value = eigenvalues[row]
                column_value = eigenvalues[column]
                if abs(column_value - row_value) < 0.00001 * temperature_as_energy:
                    chi_curie['z'] += j_z_square * current_bolzmann
                    chi_curie['x'] += (0.25 * (j_plus_square + j_minus_square) * current_bolzmann)
                else:
                    chi_van_vleck['z'] += (2 * j_z_square *
                                           current_bolzmann / (column_value - row_value))
                    chi_van_vleck['x'] += (0.5 * (j_plus_square + j_minus_square) *
                                           current_bolzmann / (column_value - row_value))
        for key in ('z', 'x'):
            chi_curie[key] = coefficient / temperature_as_energy * chi_curie[key]
            chi_van_vleck[key] = coefficient * chi_van_vleck[key]
        return chi_curie, chi_van_vleck

    def get_chi_dependence(self, temperatures=None):
        """Calculates the susceptibility at a specified range of temperatures."""
        temperatures = utils.get_new_if_old_is_none(temperatures,
                                                    linspace(1, 300, 300, dtype='float64'))
        self.get_transition_probabilities()
        temperatures = utils.get_new_if_old_is_none(temperatures,
                                                    linspace(1, 300, 300, dtype='float64'))
        chi_curie = {'z': None, 'x': None}
        chi_van_vleck = {'z': None, 'x': None}
        chi = {'z': None, 'x': None, 'total': None, 'inverse': None}
        for key in ('z', 'x'):
            chi_curie[key] = utils.get_empty_matrix(temperatures.shape)
            chi_van_vleck[key] = utils.get_empty_matrix(temperatures.shape)
            chi[key] = utils.get_empty_matrix(temperatures.shape)

        chi = {
            'total': utils.get_empty_matrix(temperatures.shape),
            'inverse': utils.get_empty_matrix(temperatures.shape),
        }
        for _, temperature in enumerate(temperatures):
            current_chi = self.get_chi(temperature)
            for key in ('z', 'x'):
                chi_curie[key] = current_chi[0][key]
                chi_van_vleck[key] = current_chi[1][key]
                chi[key] = chi_curie[key] + chi_van_vleck[key]
            chi['total'] = (chi['z'] + 2 * chi['x']) / 3
            chi['inverse'] = 1 / chi['total']

        return chi_curie, chi_van_vleck, chi

    def __repr__(self):
        """Method returns string representation of the object."""
        return (f"{self.__class__.__name__}" +
                f"(material=Material(" +
                f"rare_earth='{self.rare_earth.name}', " +
                f"crystal='{self.crystal}'))")

    def __str__(self):
        """Return a summary of the model parameters.
        This includes the rare earth, the CF parameters, and,
        if diagonalized, the eigenvalues and eigenfunctions."""
        output = [self.crystal,
                  f'Rare-earth ion: {self.rare_earth.name};',
                  f'Number of 4f-electrons = {self.number_of_f_electrons};',
                  f'J = {self.rare_earth.total_momentum_ground};']

        for key, value in self.parameters.items():
            if value:
                output.append(f'{key} = {value:.4f};')
        for key, value in self.magnet_field.items():
            if value:
                output.append(f'H{key} = {value:.4f};')
        j_average, magnetic_moment = self.get_moments()
        for key in ('z', 'x'):
            if magnetic_moment[key] > 1e-9:
                line_to_append = (f'<J{key}> = {j_average[key]:7.3f}; ' +
                                  f'<mu_{key}> = {magnetic_moment[key]:7.3f} mu_Bohr;')
                output.append(line_to_append)

        eigenvalues, eigenfunctions = self.get_eigenvalues_and_eigenfunctions()
        if eigenvalues.any():
            output.append('Crystal Field Eigenvalues and Eigenfunctions:')
            for column in range(eigenvalues.size):
                line = [f'{eigenvalues[column]:8.3f}: ']
                for row in range(eigenvalues.size):
                    if abs(eigenfunctions[row, column]) > 0.0001:
                        j_z = row - self.rare_earth.total_momentum_ground
                        line_to_append = (f'{utils.get_sign(eigenfunctions[row, column])}' +
                                          f'{abs(eigenfunctions[row, column]):7.4f}' +
                                          f'|{utils.get_sign(j_z)}{abs(j_z)}>')
                        line.append(line_to_append)
                output.append(' '.join(line))

        peaks = self.get_peaks()
        if peaks:
            output.append('Crystal Field Transitions:')
            if self.temperature != self.temperature_used:
                peaks = self.get_peaks()
            output.append(f'Temperature: {self.temperature} K')
            for peak in peaks:
                output.append(f'Energy: {peak[0]:8.3f} meV  Intensity: {peak[1]:8.4f}')

        return '\n'.join(output)


if __name__ == '__main__':
    MATERIAL = Material(rare_earth='Tb', crystal='YNi2')
    PARAMETERS = {'w': 1, 'x': -1}
    CUBIC_SAMPLE = CF(MATERIAL).get_llw_cubic_object(PARAMETERS)
    print(repr(CUBIC_SAMPLE))
