"""The module contains CEF class."""


import json

from numpy import linspace
from pretty_repr import RepresentableObject

from common import physics, utils
from common.path_utils import get_paths
from common.utils import UTF8File
from core.cef_parameters import BParameters
from core.custom_datatypes import MagnetField
from core.sample import Sample
from core.transitions import Transitions
from interactions import CrystalElectricField, Thermostat, Zeeman


class System(RepresentableObject):
    """
    Class defining the trivalent rare earth compound,
    its crystal field parameters and the eigenvalues and eigenfunctions
    of the CEF Hamiltonian, if it is already diagonalized.

    """

    def __init__(
            self,
            sample: Sample,
            parameters: BParameters = None,
            temperature: float = 0,
            magnet_field: MagnetField = None,
    ):
        """Initializes the CEF object or read it from a file."""
        self.sample = sample
        self.file_name = get_paths(
            data_name='parameters',
            format_name='.json',
            sample=self.sample,
        )
        self.magnet_field = magnet_field or MagnetField()
        self.parameters = parameters or BParameters()
        self.temperature = temperature
        self.interactions = {
            'CEF': CrystalElectricField(sample=sample, parameters=parameters),
            'Temperature': Thermostat(sample=sample, temperature=temperature),
            'MagnetField': Zeeman(sample=sample, magnet_field=magnet_field),
        }

    @property
    def excluded_attributes_for_repr(self) -> set[str]:
        return {'file_name', 'interactions'}

    def save_to_file(self):
        """Saves parameters of the current object to file."""
        saved_object = {
            'crystal': self.sample.crystal.name,
            'rare_earth': self.sample.rare_earth.info.symbol,
            'temperature': self.temperature.__dict__,
            'magnet_field': self.magnet_field.__dict__,
            'parameters': self.parameters.__dict__,
        }
        with UTF8File(self.file_name, mode='w') as file:
            json.dump(saved_object, file, indent=4, sort_keys=True)

    def load_data(self):
        """Loads CEF object from file"""
        with UTF8File(self.file_name) as file:
            properties = json.load(file)
        self.sample = Sample(
            crystal=properties.get('crystal'),
            rare_earth=properties.get('rare_earth'),
        )
        self.temperature = properties.get('temperature')
        self.magnet_field = BParameters(**properties.get('magnet_field'))
        self.parameters = BParameters(**properties.get('parameters'))

    def get_hamiltonian(self):
        cef = self.interactions['CEF']
        zeeman = self.interactions['MagnetField']
        return cef.get_hamiltonian() + zeeman.get_hamiltonian()

    def get_moments(self):
        """Calculates the magnetic moments of the CEF model."""
        transitions = Transitions(
            sample=self.sample,
            hamiltonian=self.get_hamiltonian(),
        )
        eigen_v, eigen_f = transitions.get_eigenvalues_and_eigenfunctions()
        j_ops, _ = transitions.get_transition_probabilities(eigen_f)
        thermal = physics.thermodynamics(self.temperature, eigen_v)
        if thermal['temperature'] > 0:
            j_average = {'z': 0, 'x': 0}
            statistic_sum = sum(thermal['boltzmann'])
            for index in range(eigen_v.size):
                j_average['z'] += (j_ops['z'][index, index] *
                                   thermal['boltzmann'][index])
                j_average['x'] += (0.5 * (j_ops['+'][index, index] +
                                          j_ops['-'][index, index]) *
                                   thermal['boltzmann'][index])
            for key, value in j_average.items():
                j_average[key] = value / statistic_sum
        else:
            j_average = {
                'z': (sum(j_ops['z'][eigen_v == 0, eigen_v == 0]) /
                      eigen_v[eigen_v == 0].size),
                'x': (sum(0.5 * (j_ops['+'][eigen_v == 0,
                                            eigen_v == 0] +
                                 j_ops['-'][eigen_v == 0,
                                            eigen_v == 0])) /
                      eigen_v[eigen_v == 0].size)
            }
        magnetic_moment = {}
        for key, value in j_average.items():
            magnetic_moment[key] = (
                    float(self.sample.rare_earth.info.lande_factor)
                    * value
            )
            # magnetic moments are given in units of Bohr magneton
        return j_average, magnetic_moment

    def get_chi(self):
        """Calculates the susceptibility at a specified temperature."""
        transitions = Transitions(
            sample=self.sample,
            hamiltonian=self.get_hamiltonian(),
        )
        eigen_v, eigen_f = transitions.get_eigenvalues_and_eigenfunctions()
        j_ops, _ = transitions.get_transition_probabilities(eigen_f)
        thermal = physics.thermodynamics(self.temperature, eigen_v)
        chi = {
            'curie': {'z': 0, 'x': 0},
            'van_vleck': {'z': 0, 'x': 0},
        }
        for row in range(eigen_v.size):
            for column in range(eigen_v.size):
                j_ops_square = {
                    key: val[row, column] ** 2
                    for key, val in j_ops.items()
                }
                row_value = eigen_v[row]
                column_value = eigen_v[column]
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
        coefficient = self.sample.rare_earth.lande_factor ** 2
        if thermal['temperature'] > 0:
            coefficient = coefficient / sum(thermal['boltzmann'])
        for key in ('z', 'x'):
            chi['curie'][key] = (
                    coefficient / thermal['temperature'] * chi['curie'][key]
            )
            chi['van_vleck'][key] = coefficient * chi['van_vleck'][key]
        return chi

    def get_chi_dependence(self, temperatures=None):
        """
        Calculates the susceptibility at a specified range of temperatures.

        """
        temperatures = utils.get_default(
            temperatures,
            linspace(1, 300, 300, dtype='float64')
        )
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
            self.temperature = temperature
            current_chi = self.get_chi()
            for key in ('z', 'x'):
                chi_curie[key] = current_chi['curie'][key]
                chi_van_vleck[key] = current_chi['van_vleck'][key]
                chi[key] = chi_curie[key] + chi_van_vleck[key]
            chi['total'] = (chi['z'] + 2 * chi['x']) / 3
            chi['inverse'] = 1 / chi['total']

        return chi_curie, chi_van_vleck, chi

    def __str__(self):
        """Return a summary of the model parameters.
        This includes the rare earth, the CEF parameters, and,
        if diagonalized, the eigenvalues and eigenfunctions."""
        output = [str(self.sample)]
        for interaction in self.interactions.values():
            output.append(str(interaction))
        transitions = Transitions(
            sample=self.sample,
            hamiltonian=self.get_hamiltonian()
        )
        output.append(str(transitions))
        eigen_v, _ = transitions.get_eigenvalues_and_eigenfunctions()
        peaks = transitions.get_peaks(
            self.interactions['Temperature'].get_boltzmann_factor(
                eigenvalues=eigen_v,
            )
        )
        if peaks:
            output.append('Crystal Field Transitions:')
            output.append(f'Temperature: {self.temperature} K')
            for peak in peaks:
                output.append(
                    f'Energy: {peak[0]:8.3f} meV  Intensity: {peak[1]:8.4f}'
                )

        return '\n'.join(output)
