"""The module contains CEF class."""


import numpy as np
from pretty_repr import RepresentableObject
from scipy.linalg import eigh

from common import physics, utils
from common.constants import RESOLUTION, THRESHOLD
from core.sample import Sample


class Transitions(RepresentableObject):

    def __init__(self, sample: Sample, hamiltonian: np.ndarray):
        """Initializes the CEF object or read it from a file."""
        self.sample = sample
        self.hamiltonian = hamiltonian

    def get_eigenvalues_and_eigenfunctions(
            self,
            is_ground_state_zero: bool = True,
    ):
        """Return eigenvalues and eigenfunctions of the hamiltonian."""
        eigenvalues, eigenfunctions = eigh(self.hamiltonian)
        if is_ground_state_zero:
            eigenvalues = eigenvalues - min(eigenvalues)
        return eigenvalues, eigenfunctions

    def get_transition_probabilities(self, eigenfunctions: np.ndarray):
        """
        Determines matrix elements for dipole transitions
        between eigenfunctions of the total Hamiltonian.

        """
        momentum = self.sample.rare_earth.info.total_momentum_ground
        squared_momentum = self.sample.rare_earth.squared_momentum
        size = self.sample.rare_earth.matrix_size
        j_ops = {
            'z': utils.get_empty_matrix(size),
            '+': utils.get_empty_matrix(size),
            '-': utils.get_empty_matrix(size),
        }
        transition_probability = utils.get_empty_matrix(size)
        for row in range(size):
            j_ops['z'][row, row] += (
                eigenfunctions[size - 1, row] ** 2 *
                (size - 1 - momentum)
            )
            for row_j in range(size - 1):
                j_ops['z'][row, row] += (
                    eigenfunctions[row_j, row] ** 2 *
                    (row_j - momentum)
                )
                j_ops['+'][row, row] += (
                    eigenfunctions[row_j + 1, row] *
                    eigenfunctions[row_j, row] *
                    np.sqrt(
                        squared_momentum
                        - (row_j - momentum)
                        * (row_j - momentum + 1)
                    )
                )

            j_ops['-'][row, row] = j_ops['+'][row, row]
            for column in range(row + 1, size):
                mqn_1 = size - 1 - momentum
                j_ops['z'][row, column] += (
                    eigenfunctions[size - 1, row] *
                    eigenfunctions[size - 1, column] * mqn_1
                )
                for row_j in range(size - 1):
                    mqn_1 = row_j - momentum
                    j_ops['z'][row, column] += (
                        eigenfunctions[row_j, row] *
                        eigenfunctions[row_j, column] * mqn_1
                    )
                    column_j = row_j + 1
                    mqn_2 = column_j - momentum
                    common_root = np.sqrt(squared_momentum - mqn_1 * mqn_2)
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

    def get_all_peaks(self, boltzmann_factors: np.ndarray):
        """
        Determines the peak energies and intensities
        from the total Hamiltonian.

        """
        size = self.sample.rare_earth.matrix_size
        eigen_v, eigen_f = self.get_eigenvalues_and_eigenfunctions()
        peaks = []
        _, transition_probabilities = self.get_transition_probabilities(
            eigenfunctions=eigen_f
        )
        for level_1 in range(size):
            for level_2 in range(size):
                intensity_of_transition = (
                    transition_probabilities[level_2, level_1] *
                    boltzmann_factors[level_1]
                )
                if intensity_of_transition > 0:
                    peaks.append({
                        'energy': eigen_v[level_2] - eigen_v[level_1],
                        'intensity': intensity_of_transition
                    })
        return peaks

    def get_peaks(self, boltzmann_factors: np.ndarray):
        """Returns peaks for non-degenerate levels."""
        result = []
        peaks = self.get_all_peaks(boltzmann_factors)
        for peak in peaks:
            sum_peaks = peak['energy'] * peak['intensity']
            for other_peak in peaks:
                if (
                        peak['intensity'] > 0 and
                        other_peak is not peak and
                        abs(peak['energy'] - other_peak['energy']) < RESOLUTION
                ):
                    peak['intensity'] += other_peak['intensity']
                    sum_peaks += other_peak['energy'] * other_peak['intensity']
                    other_peak['intensity'] = 0
            if peak['intensity'] > THRESHOLD:
                peak['energy'] = sum_peaks / peak['intensity']
                result.append((peak['energy'], peak['intensity']))
        result.sort()
        intensity_sum = 2 * (
                self.sample.rare_earth.info.total_momentum_ground *
                (self.sample.rare_earth.info.total_momentum_ground + 1)
        ) / 3
        intensities = [item[1] for item in result]
        if sum(intensities) != intensity_sum:
            result[0] = (result[0][0], intensity_sum - sum(intensities[1:]))

        return result

    def get_energies(self, boltzmann_factors: np.ndarray):
        """Returns transition energies"""
        return [peak[0] for peak in self.get_peaks(boltzmann_factors)]

    def get_intensities(self, boltzmann_factors: np.ndarray):
        """Returns transition intensities"""
        return [peak[1] for peak in self.get_peaks(boltzmann_factors)]

    def get_spectrum(
            self,
            boltzmann_factors: np.ndarray,
            energies=None,
            width_dict: dict = None,
    ):
        """Calculates the neutron scattering cross-section."""
        peaks = self.get_peaks(boltzmann_factors)
        eigenvalues, _ = self.get_eigenvalues_and_eigenfunctions()

        if energies is None:
            # 501 numbers in range from -1.1*E_max to 1.1*E_max
            energies = np.linspace(
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

        spectrum *= 72.65 * self.sample.rare_earth.lande_factor ** 2

        return spectrum

    def __str__(self):
        output = []
        momentum = self.sample.rare_earth.info.total_momentum_ground
        eigen_v, eigen_f = self.get_eigenvalues_and_eigenfunctions()
        if eigen_v.any():
            output.append('Crystal Field Eigenvalues and Eigenfunctions:')
            for column in range(eigen_v.size):
                line = [f'{eigen_v[column]:8.3f}: ']
                for row in range(eigen_v.size):
                    if abs(eigen_f[row, column]) > 0.0001:
                        j_z = row - momentum
                        line_to_append = (
                                f'{utils.get_sign(eigen_f[row, column])}' +
                                f'{abs(eigen_f[row, column]):7.4f}' +
                                f'|{utils.get_sign(j_z)}{abs(j_z)}>'
                        )
                        line.append(line_to_append)
                output.append(' '.join(line))

        return '\n'.join(output)
