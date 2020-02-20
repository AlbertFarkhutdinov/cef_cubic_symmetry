"""The module contains CF class."""
import os
import configparser
from numpy import sqrt, linspace
from scipy.linalg import eigh
from cef_object_scripts import common
from cef_object_scripts.tabular_information import BOHR_MAGNETON, RARE_EARTHS

RARE_EARTHS_NAMES = [element.name for element in RARE_EARTHS]
PARSER = configparser.ConfigParser()


class CF(object):
    """Class defining the trivalent rare earth compound and its crystal field parameters.
    There is only one basic object, defining the trivalent rare earth,
    its crystal field parameters, and, if already diagonalized,
    the eigenvalues and eigenfunctions of the CF Hamiltonian."""

    def __init__(self, rare_earth=None, name=None, par_file=None):
        """Initializes the CF object or read it from a file."""
        if not par_file:
            if rare_earth:
                if rare_earth in RARE_EARTHS_NAMES:
                    self.number_of_f_electrons = RARE_EARTHS_NAMES.index(rare_earth) + 1
                    self.rare_earth = RARE_EARTHS[self.number_of_f_electrons - 1]
            else:
                self.rare_earth = rare_earth[0]
            self.parameters = {
                'B20': 0, 'B40': 0, 'B60': 0,
                'B22': 0, 'B42': 0, 'B62': 0,
                'B43': 0, 'B63': 0,
                'B44': 0, 'B64': 0,
                'B66': 0,
            }
            self.magnet_field = {'z': 0, 'x': 0}
            self.name = name
        else:
            if not os.path.exists(par_file):
                raise OSError(f'"{os.path.realpath(par_file)}" does not exist')
            PARSER.read(par_file)
            self.name = PARSER.get('material', 'name')
            self.number_of_f_electrons = int(PARSER.get('material', 'number_of_f_electrons'))
            self.rare_earth = RARE_EARTHS[self.number_of_f_electrons - 1]
            for key in self.parameters.keys():
                self.parameters[key] = PARSER.getfloat('parameters', key)
            for key in self.magnet_field.keys():
                self.magnet_field[key] = PARSER.getfloat('parameters', f'H{key}')

        size = self.rare_earth.matrix_size
        self.H = common.get_empty_matrix(size)
        self.j_z = common.get_empty_matrix(size)
        self.j_plus = common.get_empty_matrix(size)
        self.j_minus = common.get_empty_matrix(size)
        self.transition_probability = common.get_empty_matrix(size)
        self.eigenfunctions = common.get_empty_matrix(size)
        self.eigenvalues = common.get_empty_matrix(size, dimension=1)
        self.peaks = []
        self.temperature = 0
        self.temperature_used = self.temperature
        self.resolution = 0.01
        self.threshold = 0.0001
        self.j_average = {'z': 0, 'x': 0}
        self.magnetic_moment = {'z': 0, 'x': 0}

    def save_to_par_file(self, par_file=None):
        """Stores the current object for later use."""
        PARSER.optionxform = str
        PARSER.add_section('material')
        PARSER.set('material', 'name', str(self.name))
        PARSER.set('material', 'rare_earth', str(self.rare_earth.name))
        PARSER.set('material', 'number_of_f_electrons', str(RARE_EARTHS.index(self.rare_earth) + 1))
        PARSER.add_section('parameters')
        for key, value in self.parameters.items():
            PARSER.set('parameters', key, str(value))
        for key, value in self.magnet_field.items():
            PARSER.set('parameters', f'H{key}', str(value))
        if par_file is None:
            par_file = f'{self.name}.cfg'
        common.check_path(par_file)
        with open(par_file, 'w', encoding='utf-8') as file:
            PARSER.write(file)

    def cef_hamiltonian(self):
        """Determines the CF Hamiltonian based on the input parameters."""
        size = self.rare_earth.matrix_size
        parameters = self.parameters
        cef_hamiltonian = common.get_empty_matrix(size)
        j = self.rare_earth.total_momentum_ground
        j_module_square = j * (j + 1)

        def down_operator(initial_number, degree_of_operator):
            result_of_lowering = 1
            for step in range(degree_of_operator):
                result_of_lowering *= (j_module_square -
                                       (initial_number - step) *
                                       (initial_number - step - 1))
            return sqrt(result_of_lowering)

        for row in range(size):  # row = 0...2J
            m = row - j  # m = -J...J
            m2 = m ** 2
            m4 = m ** 4
            o20 = 3 * m2 - j_module_square
            o40 = (35 * m4 -
                   30 * j_module_square * m2 +
                   25 * m2 -
                   6 * j_module_square +
                   3 * j_module_square ** 2)
            o60 = (231 * m ** 6 -
                   315 * j_module_square * m4 +
                   735 * m4 +
                   105 * j_module_square ** 2 * m2 -
                   525 * j_module_square * m2 +
                   294 * m2 -
                   5 * j_module_square ** 3 +
                   40 * j_module_square ** 2 -
                   60 * j_module_square)
            cef_hamiltonian[row, row] += (parameters['B20'] * o20 + parameters['B40'] * o40 +
                                          parameters['B60'] * o60)
            for degree in range(2, 7):
                if row < (size - degree):
                    column = row + degree
                    n = m + degree
                    if degree == 2:
                        o22 = 0.5 * down_operator(n, degree)
                        o42 = (3.5 * (m2 + n ** 2) - j_module_square - 5) * o22
                        o62 = (16.5 * (m4 + n ** 4) -
                               9 * (m2 + n ** 2) * j_module_square -
                               61.5 * (m2 + n ** 2) +
                               j_module_square ** 2 +
                               10 * j_module_square +
                               102) * o22
                        cef_hamiltonian[row, column] += (parameters['B22'] * o22 +
                                                         parameters['B42'] * o42 +
                                                         parameters['B62'] * o62)
                    if degree == 3:
                        o43 = 0.25 * down_operator(n, degree) * (m + n)
                        o63 = (0.25 * (11 * (m ** 3 + n ** 3) - 3 * (m + n) * j_module_square -
                                       59 * (m + n)) * down_operator(n, degree))
                        cef_hamiltonian[row, column] += (parameters['B43'] * o43 +
                                                         parameters['B63'] * o63)
                    if degree == 4:
                        o44 = 0.5 * down_operator(n, degree)
                        o64 = (5.5 * (m2 + n ** 2) - j_module_square - 38) * o44
                        cef_hamiltonian[row, column] += (parameters['B44'] * o44 +
                                                         parameters['B64'] * o64)
                    if degree == 6:
                        column = row + 6
                        n = m + 6
                        o66 = 0.5 * down_operator(n, degree)
                        cef_hamiltonian[row, column] += parameters['B66'] * o66
                    cef_hamiltonian[column, row] = cef_hamiltonian[row, column]

        return cef_hamiltonian

    def zeeman_hamiltonian(self):
        """Determines the Zeeman terms to the Hamiltonian."""
        size = self.rare_earth.matrix_size
        lande_factor = self.rare_earth.lande_factor
        j = self.rare_earth.total_momentum_ground
        magnet_field = self.magnet_field
        j_2 = j * (j + 1)
        hamiltonian = common.get_empty_matrix(size)
        for row in range(size):
            m = row - j
            hamiltonian[row, row] -= lande_factor * BOHR_MAGNETON * m * magnet_field['z']
            if row < (size - 1):
                column = row + 1
                n = m + 1
                hamiltonian[row, column] -= (0.5 * lande_factor * BOHR_MAGNETON *
                                             sqrt((j_2 - m * n)) * magnet_field['x'])
                hamiltonian[column, row] = hamiltonian[row, column]

        return hamiltonian

    def total_hamiltonian(self):
        """Returns the total Hamiltonian including CF and Zeeman terms."""
        return self.cef_hamiltonian() + self.zeeman_hamiltonian()

    def eigen_calculations(self):
        """Calculates eigenvalues and eigenfunctions of the total Hamiltonian."""
        hamiltonian = self.total_hamiltonian()
        self.eigenvalues, self.eigenfunctions = eigh(hamiltonian)
        self.eigenvalues = self.eigenvalues - min(self.eigenvalues)  # E = 0 - minimum of energy.

    def transition_probabilities(self):
        """Determines matrix elements for dipole transitions
        between eigenfunctions of the total Hamiltonian."""
        j = self.rare_earth.total_momentum_ground
        j_2 = j * (j + 1)
        size = self.rare_earth.matrix_size
        eigenfunctions = self.eigenfunctions
        j_z = common.get_empty_matrix(size)
        j_plus = common.get_empty_matrix(size)
        j_minus = common.get_empty_matrix(size)
        transition_probability = common.get_empty_matrix(size)

        for row in range(size):
            j_z[row, row] += (eigenfunctions[size - 1, row] ** 2) * (size - 1 - j)
            for row_j in range(size - 1):
                j_z[row, row] += (eigenfunctions[row_j, row] ** 2) * (row_j - j)
                j_plus[row, row] += (eigenfunctions[row_j + 1, row] * eigenfunctions[row_j, row] *
                                     sqrt(j_2 - (row_j - j) * (row_j - j + 1)))

            j_minus[row, row] = j_plus[row, row]
            for column in range(row + 1, size):
                m = size - 1 - j
                j_z[row, column] += (eigenfunctions[size - 1, row] *
                                     eigenfunctions[size - 1, column] * m)
                for row_j in range(size - 1):
                    m = row_j - j
                    j_z[row, column] += (eigenfunctions[row_j, row] *
                                         eigenfunctions[row_j, column] * m)
                    column_j = row_j + 1
                    n = column_j - j
                    common_root = sqrt(j_2 - m * n)
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

        self.j_z = j_z
        self.j_plus = j_plus
        self.j_minus = j_minus
        self.transition_probability = transition_probability

    def get_peaks(self, temperature=None, magnet_field_x=None, magnet_field_z=None):
        """Determines the peak intensities from the total Hamiltonian."""
        size = self.rare_earth.matrix_size
        old_magnet_field = {
            'x': self.magnet_field['x'],
            'z': self.magnet_field['z'],
        }
        magnet_field = {
            'x': magnet_field_x,
            'z': magnet_field_z,
        }
        for key, value in magnet_field.items():
            if value:
                self.magnet_field[key] = value
        self.eigen_calculations()
        self.transition_probabilities()

        bolzmann_factor = common.get_empty_matrix(size, dimension=1)
        temperature = common.get_temperature(temperature, self.temperature)
        thermal = common.thermodynamics(temperature, self.eigenvalues)
        if thermal['temperature'] <= 0:
            bolzmann_factor[0] = 1
        else:
            statistic_sum = sum(thermal['bolzmann'])
            bolzmann_factor = thermal['bolzmann'] / statistic_sum

        peaks = []
        for level_1 in range(size):
            for level_2 in range(size):
                energy_of_transition = self.eigenvalues[level_2] - self.eigenvalues[level_1]
                intensity_of_transition = (self.transition_probability[level_2, level_1] *
                                           bolzmann_factor[level_1])
                if intensity_of_transition > 0:
                    peaks.append({'energy': energy_of_transition,
                                  'intensity': intensity_of_transition})

        self.peaks = []
        for peak in peaks:
            if peak['intensity'] > 0:
                sum_peaks = peak['energy'] * peak['intensity']
                for other_peak in peaks:
                    if other_peak is not peak:
                        if (other_peak['intensity'] > 0 and
                                abs(peak['energy'] - other_peak['energy']) < self.resolution):
                            peak['intensity'] += other_peak['intensity']
                            sum_peaks += other_peak['energy'] * other_peak['intensity']
                            other_peak['intensity'] = 0

                if peak['intensity'] > self.threshold:
                    peak['energy'] = sum_peaks / peak['intensity']
                    self.peaks.append((peak['energy'], peak['intensity']))
        self.peaks.sort()
        self.temperature_used = temperature
        for key, value in magnet_field.items():
            if value:
                self.magnet_field[key] = old_magnet_field[key]

        return self.peaks

    def get_energies(self):
        """Returns transition energies"""
        return [peak[0] for peak in self.get_peaks()]

    def get_intensities(self):
        """Returns transition intensities"""
        return [peak[1] for peak in self.get_peaks()]

    def spectrum(self, energy=None, sigma=None, gamma=None, temperature=None,
                 magnet_field_x=None, magnet_field_z=None):
        """Calculates the neutron scattering cross section."""
        temperature = common.get_temperature(temperature, self.temperature)
        peaks = self.get_peaks(temperature, magnet_field_x, magnet_field_z)

        if energy is None:
            # 501 numbers in range from -1.1*E_max to 1.1*E_max
            energy = linspace(-1.1 * self.eigenvalues[-1], 1.1 * self.eigenvalues[-1], 501)
        if sigma is None and gamma is None:
            sigma = 0.01 * (max(energy) - min(energy))

        spectrum = common.get_empty_matrix(energy.size, dimension=1)

        for peak in peaks:
            if sigma and not gamma:
                spectrum += peak[1] * common.gauss(energy, peak[0], sigma)
            elif gamma and not sigma:
                spectrum += peak[1] * common.lorentz(energy, peak[0], gamma)
            elif sigma and gamma:
                spectrum += peak[1] * common.pseudo_voigt(energy, peak[0], sigma, gamma)

        spectrum *= 72.65 * self.rare_earth.lande_factor ** 2

        return spectrum

    def get_moments(self, temperature=None):
        """Calculates the magnetic moments of the CF model."""
        temperature = common.get_temperature(temperature, self.temperature)
        if temperature != self.temperature_used:
            self.get_peaks(temperature)
        thermal = common.thermodynamics(temperature, self.eigenvalues)
        if thermal['temperature'] > 0:
            j_average = {'z': 0, 'x': 0}
            statistic_sum = sum(thermal['bolzmann'])
            for index in range(self.eigenvalues.size):
                j_average['z'] += (self.j_z[index, index] *
                                   thermal['bolzmann'][index])
                j_average['x'] += (0.5 * (self.j_plus[index, index] + self.j_minus[index, index]) *
                                   thermal['bolzmann'][index])
            for key, value in j_average.items():
                j_average[key] = value / statistic_sum
        else:
            j_average = {
                'z': (sum(self.j_z[self.eigenvalues == 0, self.eigenvalues == 0]) /
                      self.eigenvalues[self.eigenvalues == 0].size),
                'x': (sum(0.5 * (self.j_plus[self.eigenvalues == 0, self.eigenvalues == 0] +
                                 self.j_minus[self.eigenvalues == 0, self.eigenvalues == 0])) /
                      self.eigenvalues[self.eigenvalues == 0].size)
            }
        for key, value in j_average.items():
            self.j_average[key] = value
            self.magnetic_moment[key] = self.rare_earth.lande_factor * value
            # magnetic moments are given in units of Bohr magneton
        return self.j_average, self.magnetic_moment

    def chi(self, temperature=None):
        """Calculates the susceptibility at a specified temperature."""
        temperature = common.get_temperature(temperature, self.temperature)
        thermal = common.thermodynamics(temperature, self.eigenvalues)
        temperature_as_energy = thermal['temperature']
        statistic_sum = 1
        if temperature_as_energy > 0:
            statistic_sum = sum(thermal['bolzmann'])

        coefficient = self.rare_earth.lande_factor ** 2 / statistic_sum
        chi_curie = {'z': 0, 'x': 0}
        chi_van_vleck = {'z': 0, 'x': 0}
        for row in range(self.eigenvalues.size):
            for column in range(self.eigenvalues.size):
                current_bolzmann = thermal['bolzmann'][row]
                j_z_square = self.j_z[row, column] ** 2
                j_plus_square = self.j_plus[row, column] ** 2
                j_minus_square = self.j_minus[row, column] ** 2
                row_value = self.eigenvalues[row]
                column_value = self.eigenvalues[column]
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

    def chi_s(self, temperatures=None):
        """Calculates the susceptibility at a specified range of temperatures."""
        temperatures = common.get_temperature(temperatures, linspace(1, 300, 300, dtype='float64'))
        self.eigen_calculations()
        self.transition_probabilities()
        temperatures = common.get_temperature(temperatures, linspace(1, 300, 300, dtype='float64'))
        chi_curie = {'z': None, 'x': None}
        chi_van_vleck = {'z': None, 'x': None}
        chi = {'z': None, 'x': None, 'total': None, 'inverse': None}
        for key in ('z', 'x'):
            chi_curie[key] = common.get_empty_matrix(temperatures.shape)
            chi_van_vleck[key] = common.get_empty_matrix(temperatures.shape)
            chi[key] = common.get_empty_matrix(temperatures.shape)

        chi = {
            'total': common.get_empty_matrix(temperatures.shape),
            'inverse': common.get_empty_matrix(temperatures.shape),
        }
        for index, temperature in enumerate(temperatures):
            current_chi = self.chi(temperatures[index])
            for key in ('z', 'x'):
                chi_curie[key] = current_chi[0][key]
                chi_van_vleck[key] = current_chi[1][key]
                chi[key] = chi_curie[key] + chi_van_vleck[key]
            chi['total'] = (chi['z'] + 2 * chi['x']) / 3
            chi['inverse'] = 1 / chi['total']

        return chi_curie, chi_van_vleck, chi

    def __str__(self):
        """Return a summary of the model parameters.
        This includes the rare earth, the CF parameters, and,
        if diagonalized, the eigenvalues and eigenvectors."""
        output = []
        if self.name:
            output.append(self.name)
        output.append(f'Rare-earth ion: {self.rare_earth.name};\n'
                      f'Number of 4f-electrons = {self.number_of_f_electrons};\n'
                      f'J = {self.rare_earth.total_momentum_ground};')

        for key, value in self.parameters.items():
            if value:
                output.append(f'{key} = {value:.4f};')
        for key, value in self.magnet_field.items():
            if value:
                output.append(f'H{key} = {value:.4f};')
        self.get_peaks()
        j_average, magnetic_moment = self.get_moments()
        for key in ('z', 'x'):
            if magnetic_moment[key] > 1e-9:
                line_to_append = (f'<J{key}> = {j_average[key]:7.3f}; ' +
                                  f'<mu_{key}> = {magnetic_moment[key]:7.3f} mu_Bohr;')
                output.append(line_to_append)

        if self.eigenvalues.any():
            output.append('Crystal Field Eigenvalues and Eigenfunctions:')
            for column in range(self.eigenvalues.size):
                line = [f'{self.eigenvalues[column]:8.3f}: ']
                for row in range(self.eigenvalues.size):
                    if abs(self.eigenfunctions[row, column]) > 0.0001:
                        j_z = row - self.rare_earth.total_momentum_ground
                        sign1 = common.get_sign(self.eigenfunctions[row, column])
                        sign2 = common.get_sign(j_z)
                        line_to_append = (f'{sign1}' +
                                          f'{abs(self.eigenfunctions[row, column]):7.4f}' +
                                          f'|{sign2}{abs(j_z)}>')
                        line.append(line_to_append)
                output.append(' '.join(line))

        if self.peaks:
            output.append('Crystal Field Transitions:')
            if self.temperature != self.temperature_used:
                self.get_peaks()
            output.append(f'Temperature: {self.temperature} K')
            for peak in self.peaks:
                output.append(f'Energy: {peak[0]:8.3f} meV  Intensity: {peak[1]:8.4f}')

        return '\n'.join(output)


if __name__ == '__main__':
    from cProfile import run
    from cef_object_scripts.get_results import get_object_with_parameters
    MATERIAL = {'crystal': 'YNi2', 'rare_earth': 'Tb'}
    PARAMETERS = {'w': 1, 'x': -1}
    CUBIC_SAMPLE = get_object_with_parameters(MATERIAL, PARAMETERS)
    CUBIC_SAMPLE.eigen_calculations()
    run('CUBIC_SAMPLE.transition_probabilities()')
