# Questions to think:
# What difference between numpy.linalg.eigh and scipy.linalg.eigh

import os
import configparser
import numpy
from scipy.linalg import eigh
from fractions import Fraction
from cef_object_scripts import common

parser = configparser.ConfigParser()
bohr_magneton = 5.788382e-2  # meV/T
# list of RE ions names
rare_earths = ['Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb']
# list of quantum numbers for total angular momentum of the ground state multiplet for RE ions
total_momentum_ground = [2.5, 4, 4.5, 4, 2.5, 0, 3.5, 6, 7.5, 8, 7.5, 6, 3.5]
# list of Lande g-factors of the ground state multiplet for RE ions
lande_factor = [Fraction(6, 7), Fraction(4, 5), Fraction(8, 11), Fraction(3, 5),
                Fraction(2, 7), 0, 2, Fraction(3, 2), Fraction(4, 3), Fraction(5, 4),
                Fraction(6, 5), Fraction(7, 6), Fraction(8, 7)]
# <r^2> radial integral for 3+ f-electrons
r2 = [0.3666, 0.3380, 0.3120, 0.2917, 0.2728, 0.2569, 0.2428, 0.2302, 0.2188, 0.2085, 0.1991, 0.1905, 0.1826]
# <r^4> radial integral for 3+ f-electrons
r4 = [0.3108, 0.2670, 0.2015, 0.1488, 0.1772, 0.1584, 0.1427, 0.1295, 0.1180, 0.1081, 0.0996, 0.0921, 0.0854]
# <r^6> radial integral for 3+ f-electrons
r6 = [0.5119, 0.4150, 0.3300, 0.2787, 0.2317, 0.1985, 0.1720, 0.1505, 0.1328, 0.1810, 0.1058, 0.0953, 0.0863]
# Stevens factors value from article of M.T.Hutchings (1964)
# alpha_j - second-degree Stevens factors for 3+ f-electrons
stevens_factor2 = [Fraction(-2, 35), Fraction(-52, 2475), Fraction(-7, 1089), Fraction(14, 1815),
                   Fraction(13, 315), 0, 0, Fraction(-1, 99), Fraction(-2, 315), Fraction(-1, 450),
                   Fraction(4, 1575), Fraction(1, 99), Fraction(2, 63)]
# beta_j - fourth-degree Stevens factors for 3+ f-electrons
stevens_factor4 = [Fraction(2, 315), Fraction(-4, 5445), Fraction(-136, 467181), Fraction(952, 2335905),
                   Fraction(26, 10395), 0, 0, Fraction(2, 16335), Fraction(-8, 135135), Fraction(-1, 30030),
                   Fraction(2, 45045), Fraction(8, 49005), Fraction(-2, 1155)]
# gamma_j - sixth-degree Stevens factors for 3+ f-electrons
stevens_factor6 = [0, Fraction(272, 4459455), Fraction(-1615, 42513471), Fraction(2584, 3864861),
                   0, 0, 0, Fraction(-1, 891891), Fraction(4, 3864861), Fraction(-5, 3864861),
                   Fraction(8, 3864861), Fraction(-5, 891891), Fraction(4, 27027)]


class CF(object):
    # Class defining the trivalent rare earth compound and its crystal field parameters. 
    # There is only one basic object, defining the trivalent rare earth, its crystal field parameters, and,
    # if already diagonalized, the eigenvalues and eigenfunctions of the CF Hamiltonian.

    def __init__(self, rare_earth=None, name=None, par_file=None):
        # Initializes the CF object or read it from a file.
        if not par_file:
            if rare_earth:
                self.rare_earth = rare_earth
            else:
                self.rare_earth = 'Ce'
            self.B20 = 0
            self.B40 = 0
            self.B60 = 0
            self.B22 = 0
            self.B42 = 0
            self.B62 = 0
            self.B43 = 0
            self.B63 = 0
            self.B44 = 0
            self.B64 = 0
            self.B66 = 0
            self.magnet_field_z = 0
            self.magnet_field_x = 0
            self.name = name
        else:
            self.load_from_par_file(par_file)

        self.index_of_rare_earth = rare_earths.index(self.rare_earth)
        self.number_of_f_electrons = self.index_of_rare_earth + 1
        self.j = total_momentum_ground[self.index_of_rare_earth]
        self.size = int(2 * self.j + 1)
        self.lande_factor = lande_factor[self.index_of_rare_earth]
        self.r2 = r2[self.index_of_rare_earth]
        self.r4 = r4[self.index_of_rare_earth]
        self.r6 = r6[self.index_of_rare_earth]
        self.stevens_factor2 = stevens_factor2[self.index_of_rare_earth]
        self.stevens_factor4 = stevens_factor4[self.index_of_rare_earth]
        self.stevens_factor6 = stevens_factor6[self.index_of_rare_earth]
        self.H = numpy.zeros((self.size, self.size), dtype='float64')
        self.j_z = numpy.zeros((self.size, self.size), dtype='float64')
        self.j_plus = numpy.zeros((self.size, self.size), dtype='float64')
        self.j_minus = numpy.zeros((self.size, self.size), dtype='float64')
        self.transition_probability = numpy.zeros((self.size, self.size), dtype='float64')
        self.eigenfunctions = numpy.zeros((self.size, self.size), dtype='float64')
        self.eigenvalues = numpy.zeros(self.size, dtype='float64')
        self.peaks = []
        self.temperature = 0
        self.temperature_used = self.temperature
        self.resolution = 0.01
        self.threshold = 0.0001
        self.j_z_average = 0
        self.j_x_average = 0
        self.magnetic_moment_z = 0
        self.magnetic_moment_x = 0

    def save_to_par_file(self, par_file=None):
        # Store the current object for later use.
        parser.optionxform = str
        parser.add_section('material')
        parser.set('material', 'name', str(self.name))
        parser.set('material', 'rare_earth', str(self.rare_earth))
        parser.add_section('parameters')
        parser.set('parameters', 'B20', str(self.B20))
        parser.set('parameters', 'B22', str(self.B22))
        parser.set('parameters', 'B40', str(self.B40))
        parser.set('parameters', 'B42', str(self.B42))
        parser.set('parameters', 'B43', str(self.B43))
        parser.set('parameters', 'B44', str(self.B44))
        parser.set('parameters', 'B60', str(self.B60))
        parser.set('parameters', 'B62', str(self.B62))
        parser.set('parameters', 'B63', str(self.B63))
        parser.set('parameters', 'B64', str(self.B64))
        parser.set('parameters', 'B66', str(self.B66))
        parser.set('parameters', 'Hz', str(self.magnet_field_z))
        parser.set('parameters', 'Hx', str(self.magnet_field_x))
        if par_file is None:
            par_file = f'{self.name}.cfg'
        common.check_path(par_file)
        with open(par_file, 'w', encoding='utf-8') as file:
            parser.write(file)

    def load_from_par_file(self, par_file):
        # Load a saved version of the current object.
        if not os.path.exists(par_file):
            raise OSError(f'"{os.path.realpath(par_file)}" does not exist')

        parser.read(par_file)
        self.name = parser.get('material', 'name')
        self.rare_earth = parser.get('material', 'rare_earth')
        self.B20 = parser.getfloat('parameters', 'B20')
        self.B22 = parser.getfloat('parameters', 'B22')
        self.B40 = parser.getfloat('parameters', 'B40')
        self.B42 = parser.getfloat('parameters', 'B42')
        self.B43 = parser.getfloat('parameters', 'B43')
        self.B44 = parser.getfloat('parameters', 'B44')
        self.B60 = parser.getfloat('parameters', 'B60')
        self.B62 = parser.getfloat('parameters', 'B62')
        self.B63 = parser.getfloat('parameters', 'B63')
        self.B64 = parser.getfloat('parameters', 'B64')
        self.B66 = parser.getfloat('parameters', 'B66')
        self.magnet_field_z = parser.getfloat('parameters', 'Hz')
        self.magnet_field_x = parser.getfloat('parameters', 'Hx')
        self.index_of_rare_earth = rare_earths.index(self.rare_earth)
        self.number_of_f_electrons = self.index_of_rare_earth + 1
        self.j = total_momentum_ground[self.index_of_rare_earth]
        self.size = int(2 * self.j + 1)
        self.lande_factor = lande_factor[self.index_of_rare_earth]
        self.r2 = r2[self.index_of_rare_earth]
        self.r4 = r4[self.index_of_rare_earth]
        self.r6 = r6[self.index_of_rare_earth]
        self.stevens_factor2 = stevens_factor2[self.index_of_rare_earth]
        self.stevens_factor4 = stevens_factor4[self.index_of_rare_earth]
        self.stevens_factor6 = stevens_factor6[self.index_of_rare_earth]
        self.H = numpy.zeros((self.size, self.size), dtype='float64')
        self.j_z = numpy.zeros((self.size, self.size), dtype='float64')
        self.j_plus = numpy.zeros((self.size, self.size), dtype='float64')
        self.j_minus = numpy.zeros((self.size, self.size), dtype='float64')
        self.transition_probability = numpy.zeros((self.size, self.size), dtype='float64')
        self.eigenfunctions = numpy.zeros((self.size, self.size), dtype='float64')
        self.eigenvalues = numpy.zeros(self.size, dtype='float64')

    def cef_hamiltonian(self):
        # Determine the CF Hamiltonian based on the input parameters.
        cef_hamiltonian = numpy.zeros((self.size, self.size), dtype='float64')
        j = self.j
        j_module_square = j * (j + 1)

        def down_operator(initial_number, degree_of_operator):
            result_of_lowering = 1
            for step in range(degree_of_operator):
                result_of_lowering *= j_module_square - (initial_number - step) * (initial_number - step - 1)
            return numpy.sqrt(result_of_lowering)

        for row in range(self.size):  # row = 0...2J
            m = row - j  # m = -J...J
            m2 = m ** 2
            m4 = m ** 4
            o20 = 3 * m2 - j_module_square
            o40 = 35 * m4 - 30 * j_module_square * m2 + 25 * m2 - 6 * j_module_square + 3 * j_module_square ** 2
            o60 = (231 * m ** 6 - 315 * j_module_square * m4 + 735 * m4 + 105 * j_module_square ** 2 * m2 -
                   525 * j_module_square * m2 + 294 * m2 - 5 * j_module_square ** 3 + 40 * j_module_square ** 2 -
                   60 * j_module_square)
            cef_hamiltonian[row, row] += self.B20 * o20 + self.B40 * o40 + self.B60 * o60
            for degree in range(2, 7):
                if row < (self.size - degree):
                    column = row + degree
                    n = m + degree
                    if degree == 2:
                        o22 = 0.5 * down_operator(n, degree)
                        o42 = (3.5 * (m2 + n ** 2) - j_module_square - 5) * o22
                        o62 = (16.5 * (m4 + n ** 4) - 9 * (m2 + n ** 2) * j_module_square - 61.5 * (m2 + n ** 2) +
                               j_module_square ** 2 + 10 * j_module_square + 102) * o22
                        cef_hamiltonian[row, column] += self.B22 * o22 + self.B42 * o42 + self.B62 * o62
                    if degree == 3:
                        o43 = 0.25 * down_operator(n, degree) * (m + n)
                        o63 = (0.25 * (11 * (m ** 3 + n ** 3) - 3 * (m + n) * j_module_square -
                                       59 * (m + n)) * down_operator(n, degree))
                        cef_hamiltonian[row, column] += self.B43 * o43 + self.B63 * o63
                    if degree == 4:
                        o44 = 0.5 * down_operator(n, degree)
                        o64 = (5.5 * (m2 + n ** 2) - j_module_square - 38) * o44
                        cef_hamiltonian[row, column] += self.B44 * o44 + self.B64 * o64
                    if degree == 6:
                        column = row + 6
                        n = m + 6
                        o66 = 0.5 * down_operator(n, degree)
                        cef_hamiltonian[row, column] += self.B66 * o66
                    cef_hamiltonian[column, row] = cef_hamiltonian[row, column]

        return cef_hamiltonian

    def zeeman_hamiltonian(self):
        # Determine the Zeeman terms to the Hamiltonian.
        j = self.j
        j_2 = j * (j + 1)
        size = self.size
        hamiltonian = numpy.zeros((size, size), dtype='float64')
        for row in range(size):
            m = row - j
            hamiltonian[row, row] -= self.lande_factor * bohr_magneton * m * self.magnet_field_z
            if row < (size - 1):
                column = row + 1
                n = m + 1
                hamiltonian[row, column] -= (0.5 * self.lande_factor * bohr_magneton *
                                             numpy.sqrt((j_2 - m * n)) * self.magnet_field_x)
                hamiltonian[column, row] = hamiltonian[row, column]

        return hamiltonian

    def total_hamiltonian(self):
        # Returns the total Hamiltonian including CF and Zeeman terms.
        return self.cef_hamiltonian() + self.zeeman_hamiltonian()

    def eigen_calculations(self):
        # Calculate eigenvalues and eigenfunctions of the total Hamiltonian.
        hamiltonian = self.total_hamiltonian()
        self.eigenvalues, self.eigenfunctions = eigh(hamiltonian)
        self.eigenvalues = self.eigenvalues - min(self.eigenvalues)  # E = 0 - minimum of energy.

    def transition_probabilities(self):
        # Determine matrix elements for dipole transitions between eigenfunctions of the total Hamiltonian.
        j = self.j
        j_2 = j * (j + 1)
        size = self.size
        eigenfunctions = self.eigenfunctions
        self.j_z = numpy.zeros((size, size), dtype='float64')
        self.j_plus = numpy.zeros((size, size), dtype='float64')
        self.j_minus = numpy.zeros((size, size), dtype='float64')
        for row in range(size):
            for row_j in range(size):
                self.j_z[row, row] += (eigenfunctions[row_j, row] ** 2) * (row_j - j)
                if row_j < (size - 1):
                    self.j_plus[row, row] += (eigenfunctions[row_j + 1, row] * eigenfunctions[row_j, row] *
                                              numpy.sqrt(j_2 - (row_j - j) * (row_j - j + 1)))
            self.j_minus[row, row] = self.j_plus[row, row]
            for column in range(row + 1, size):
                for row_j in range(size):
                    m = row_j - j
                    self.j_z[row, column] += (eigenfunctions[row_j, row] * eigenfunctions[row_j, column] * m)
                    if row_j < (size - 1):
                        column_j = row_j + 1
                        n = column_j - j
                        common_root = numpy.sqrt(j_2 - m * n)
                        self.j_plus[row, column] += (eigenfunctions[column_j, row] * eigenfunctions[row_j, column] *
                                                     common_root)
                        self.j_minus[row, column] += (eigenfunctions[row_j, row] * eigenfunctions[column_j, column] *
                                                      common_root)
                self.transition_probability[row, column] = ((2 * self.j_z[row, column] ** 2 +
                                                             self.j_plus[row, column] ** 2 +
                                                             self.j_minus[row, column] ** 2) / 3)

                self.j_z[column, row] = self.j_z[row, column]
                self.j_plus[column, row] = self.j_minus[row, column]
                self.j_minus[column, row] = self.j_plus[row, column]
                self.transition_probability[column, row] = self.transition_probability[row, column]

    def get_peaks(self, temperature=None, magnet_field_x=None, magnet_field_z=None):
        # Determine the peak intensities from the total Hamiltonian.
        old_magnet_field_x = self.magnet_field_x
        old_magnet_field_z = self.magnet_field_z
        if magnet_field_x:
            self.magnet_field_x = magnet_field_x
        if magnet_field_z:
            self.magnet_field_z = magnet_field_z
        self.eigen_calculations()
        self.transition_probabilities()

        bolzmann_factor = numpy.zeros(self.size, dtype='float64')
        temperature = common.get_temperature(temperature, self.temperature)
        thermal = common.thermodynamics(temperature, self.eigenvalues)
        if thermal['temperature'] <= 0:
            bolzmann_factor[0] = 1
        else:
            statistic_sum = sum(thermal['bolzmann'])
            bolzmann_factor = thermal['bolzmann'] / statistic_sum

        peaks = []
        for level_1 in range(self.size):
            for level_2 in range(self.size):
                energy_of_transition = self.eigenvalues[level_2] - self.eigenvalues[level_1]
                intensity_of_transition = self.transition_probability[level_2, level_1] * bolzmann_factor[level_1]
                if intensity_of_transition > 0:
                    peaks.append({'energy': energy_of_transition, 'intensity': intensity_of_transition})

        self.peaks = []
        for peak in peaks:
            if peak['intensity'] > 0:
                sum_peaks = peak['energy'] * peak['intensity']
                for other_peak in peaks:
                    if other_peak is not peak:
                        if other_peak['intensity'] > 0 and abs(peak['energy'] - other_peak['energy']) < self.resolution:
                            peak['intensity'] += other_peak['intensity']
                            sum_peaks += other_peak['energy'] * other_peak['intensity']
                            other_peak['intensity'] = 0

                if peak['intensity'] > self.threshold:
                    peak['energy'] = sum_peaks / peak['intensity']
                    self.peaks.append((peak['energy'], peak['intensity']))
        self.peaks.sort()
        self.temperature_used = temperature
        if magnet_field_x:
            self.magnet_field_x = old_magnet_field_x
        if magnet_field_z:
            self.magnet_field_z = old_magnet_field_z

        return self.peaks

    def get_energies(self):
        return [peak[0] for peak in self.get_peaks()]

    def get_intensities(self):
        return [peak[1] for peak in self.get_peaks()]

    def spectrum(self, energy=None, sigma=None, gamma=None, temperature=None, magnet_field_x=None, magnet_field_z=None):
        # Calculates the neutron scattering cross section.
        temperature = common.get_temperature(temperature, self.temperature)
        peaks = self.get_peaks(temperature, magnet_field_x, magnet_field_z)

        if energy is None:
            # 501 numbers in range from -1.1*E_max to 1.1*E_max
            energy = numpy.linspace(-1.1 * self.eigenvalues[-1], 1.1 * self.eigenvalues[-1], 501)
        if sigma is None and gamma is None:
            sigma = 0.01 * (max(energy) - min(energy))

        spectrum = numpy.zeros(energy.size, dtype='float64')

        for peak in peaks:
            if sigma and not gamma:
                spectrum += peak[1] * common.gauss(energy, peak[0], sigma)
            elif gamma and not sigma:
                spectrum += peak[1] * common.lorentz(energy, peak[0], gamma)
            elif sigma and gamma:
                spectrum += peak[1] * common.pseudo_voigt(energy, peak[0], sigma, gamma)

        spectrum *= 72.65 * self.lande_factor ** 2

        return spectrum

    def get_moments(self, temperature=None):
        # Calculate the magnetic moments of the CF model.
        temperature = common.get_temperature(temperature, self.temperature)
        if temperature != self.temperature_used:
            self.get_peaks(temperature)
        thermal = common.thermodynamics(temperature, self.eigenvalues)
        if thermal['temperature'] > 0:
            j_z_average = 0
            j_x_average = 0
            statistic_sum = sum(thermal['bolzmann'])
            for index in range(self.eigenvalues.size):
                j_z_average += (self.j_z[index, index] *
                                thermal['bolzmann'][index])
                j_x_average += (0.5 * (self.j_plus[index, index] + self.j_minus[index, index]) *
                                thermal['bolzmann'][index])
            j_z_average = j_z_average / statistic_sum
            j_x_average = j_x_average / statistic_sum
        else:
            j_z_average = (sum(self.j_z[self.eigenvalues == 0, self.eigenvalues == 0]) /
                           self.eigenvalues[self.eigenvalues == 0].size)
            j_x_average = (sum(0.5 * (self.j_plus[self.eigenvalues == 0, self.eigenvalues == 0] +
                                      self.j_minus[self.eigenvalues == 0, self.eigenvalues == 0])) /
                           self.eigenvalues[self.eigenvalues == 0].size)
        self.j_z_average = j_z_average
        self.j_x_average = j_x_average
        # magnetic moments are given in units of Bohr magneton
        self.magnetic_moment_z = self.lande_factor * self.j_z_average
        self.magnetic_moment_x = self.lande_factor * self.j_x_average
        moments = {
            'j_x_average': self.j_x_average,
            'j_z_average': self.j_z_average,
            'magnetic_moment_x': self.magnetic_moment_x,
            'magnetic_moment_z': self.magnetic_moment_z
        }
        return moments

    def chi(self, temperature=None):
        # Calculate the susceptibility at a specified temperature.
        temperature = common.get_temperature(temperature, self.temperature)
        thermal = common.thermodynamics(temperature, self.eigenvalues)
        temperature_as_energy = thermal['temperature']
        statistic_sum = 1
        if temperature_as_energy > 0:
            statistic_sum = sum(thermal['bolzmann'])

        coefficient = self.lande_factor ** 2 / statistic_sum
        chi_curie_z = 0
        chi_curie_x = 0
        chi_van_vleck_z = 0
        chi_van_vleck_x = 0
        for row in range(self.eigenvalues.size):
            for column in range(self.eigenvalues.size):
                current_bolzmann = thermal['bolzmann'][row]
                j_z_square = self.j_z[row, column] ** 2
                j_plus_square = self.j_plus[row, column] ** 2
                j_minus_square = self.j_minus[row, column] ** 2
                row_value = self.eigenvalues[row]
                column_value = self.eigenvalues[column]
                if abs(column_value - row_value) < 0.00001 * temperature_as_energy:
                    chi_curie_z += j_z_square * current_bolzmann
                    chi_curie_x += (0.25 * (j_plus_square + j_minus_square) * current_bolzmann)
                else:
                    chi_van_vleck_z += (2 * j_z_square * current_bolzmann / (column_value - row_value))
                    chi_van_vleck_x += (0.5 * (j_plus_square + j_minus_square) * current_bolzmann /
                                        (column_value - row_value))
        return {
            'chi_curie_z': coefficient / temperature_as_energy * chi_curie_z,
            'chi_curie_x': coefficient / temperature_as_energy * chi_curie_x,
            'chi_van_vleck_z': coefficient * chi_van_vleck_z,
            'chi_van_vleck_x': coefficient * chi_van_vleck_x,
        }

    def chi_s(self, temperatures=None):
        temperatures = common.get_temperature(temperatures, numpy.linspace(1, 300, 300, dtype=numpy.float32))
        self.eigen_calculations()
        self.transition_probabilities()
        temperatures = common.get_temperature(temperatures, numpy.linspace(1, 300, 300, dtype=numpy.float32))
        chi_curie_z = numpy.zeros(shape=temperatures.shape, dtype=numpy.float32)
        chi_van_vleck_z = numpy.zeros(shape=temperatures.shape, dtype=numpy.float32)
        chi_z = numpy.zeros(shape=temperatures.shape, dtype=numpy.float32)
        chi_curie_x = numpy.zeros(shape=temperatures.shape, dtype=numpy.float32)
        chi_van_vleck_x = numpy.zeros(shape=temperatures.shape, dtype=numpy.float32)
        chi_x = numpy.zeros(shape=temperatures.shape, dtype=numpy.float32)
        chi = numpy.zeros(shape=temperatures.shape, dtype=numpy.float32)
        inverse_chi = numpy.zeros(shape=temperatures.shape, dtype=numpy.float32)
        for index in range(len(temperatures)):
            current_chi = self.chi(temperatures[index])
            chi_curie_z[index] = current_chi['chi_curie_z']
            chi_van_vleck_z[index] = current_chi['chi_van_vleck_z']
            chi_z[index] = chi_curie_z[index] + chi_van_vleck_z[index]
            chi_curie_x[index] = current_chi['chi_curie_x']
            chi_van_vleck_x[index] = current_chi['chi_van_vleck_x']
            chi_x[index] = chi_curie_x[index] + chi_van_vleck_x[index]
            chi[index] = (chi_z[index] + 2 * chi_x[index]) / 3
            inverse_chi[index] = 1 / chi[index]

        return {
            'chi_curie_z': chi_curie_z,
            'chi_van_vleck_z': chi_van_vleck_z,
            'chi_z': chi_z,
            'chi_curie_x': chi_curie_x,
            'chi_van_vleck_x': chi_van_vleck_x,
            'chi_x': chi_x,
            'chi': chi,
            'inverse_chi': inverse_chi,
        }

    def __str__(self):
        # Return a summary of the model parameters.
        # This includes the rare earth, the CF parameters, and, if diagonalized, the eigenvalues and eigenvectors.
        if self.name:
            output = [self.name]
        else:
            output = []
        output.append(f'{self.rare_earth}: Number of 4f-electrons = {self.number_of_f_electrons}, J = {self.j}')
        line = []

        def append_if_not_zero(name, value):
            if value:
                line.append(f'{name} = {value:.4f},')

        append_if_not_zero('B20', self.B20)
        append_if_not_zero('B22', self.B22)
        append_if_not_zero('B40', self.B40)
        append_if_not_zero('B42', self.B42)
        append_if_not_zero('B43', self.B43)
        append_if_not_zero('B44', self.B44)
        append_if_not_zero('B60', self.B60)
        append_if_not_zero('B62', self.B62)
        append_if_not_zero('B63', self.B63)
        append_if_not_zero('B64', self.B64)
        append_if_not_zero('B66', self.B66)
        append_if_not_zero('Hz', self.magnet_field_z)
        append_if_not_zero('Hx', self.magnet_field_x)
        if line:
            output.append(' '.join(line))
        self.get_peaks()
        j_x_average = self.get_moments()['j_x_average']
        j_z_average = self.get_moments()['j_z_average']
        magnetic_moment_x = self.get_moments()['magnetic_moment_x']
        magnetic_moment_z = self.get_moments()['magnetic_moment_z']

        if magnetic_moment_z or magnetic_moment_x:
            output.append(f'<Jz> = {j_z_average:7.3f} <mu_z> = {magnetic_moment_z:7.3f} mu_Bohr')
            output.append(f'<Jx> = {j_x_average:7.3f} <mu_x> = {magnetic_moment_x:7.3f} mu_Bohr')

        if self.eigenvalues.any():
            output.append('Crystal Field Eigenvalues and Eigenfunctions')
            for column in range(self.eigenvalues.size):
                line = [f'{self.eigenvalues[column]:8.3f}: ']
                for row in range(self.eigenvalues.size):
                    if abs(self.eigenfunctions[row, column]) > 0.0001:
                        j_z = row - self.j
                        sign1 = common.get_sign(self.eigenfunctions[row, column])
                        sign2 = common.get_sign(j_z)
                        line.append(f'{sign1}{abs(self.eigenfunctions[row, column]):7.4f}|{sign2}{abs(j_z)}>')
                output.append(' '.join(line))

        if self.peaks:
            output.append('Crystal Field Transitions')
            if self.temperature != self.temperature_used:
                self.get_peaks()
            output.append(f'Temperature: {self.temperature} K')
            for peak in self.peaks:
                output.append(f'Energy: {peak[0]:8.3f} meV  Intensity: {peak[1]:8.4f}')

        return '\n'.join(output)


if __name__ == '__main__':
    f4 = 60
    f6 = {'Pr': 1260, 'Nd': 2520, 'Pm': 1260, 'Gd': 1260, 'Tb': 7560, 'Dy': 13860, 'Ho': 13860, 'Er': 13860, 'Tm': 7560,
          'Yb': 1260}  # add check for Ce, Sm, Eu
    rare = 'Tb'
    w = 1
    x = -1
    cubic_sample = CF(name='YNi2: ' + rare + '3+', rare_earth=rare)
    cubic_sample.B40 = w * x / f4
    cubic_sample.B44 = 5 * cubic_sample.B40
    cubic_sample.B60 = w * (1 - abs(x)) / f6[rare]
    cubic_sample.B64 = -21 * cubic_sample.B60

    ham = cubic_sample.total_hamiltonian()
    for i in range(ham.shape[0]):
        for k in range(ham.shape[1]):
            print(f'{ham[i][k]: 9.3f}', end='\t')
        print()

    a = eigh(ham)
    b = numpy.linalg.eigh(ham)
    print()
    for i in range(a[0].shape[0]):
        print(f'{a[0][i]: 9.3f}', end='\t')
    print('\n')
    for i in range(a[1].shape[0]):
        for k in range(a[1].shape[1]):
            print(f'{a[1][i][k]: 9.3f}', end='\t')
        print()

    print()
    for i in range(b[0].shape[0]):
        print(f'{b[0][i]: 9.3f}', end='\t')
    print('\n')
    for i in range(b[1].shape[0]):
        for k in range(b[1].shape[1]):
            print(f'{b[1][i][k]: 9.3f}', end='\t')
        print()