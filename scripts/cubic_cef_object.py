"""
The module contains class for CEF with cubic symmetry.

"""


import sys

from numpy import linspace

from scripts.cef_object import CEF
from scripts.common.constants import CrossPoint, Material
from scripts.common.tabular_information import F4
from scripts.common.utils import (
    get_time_of_execution,
    OpenedFile,
    write_row,
    get_ratios_names,
)
from scripts.common.utils import get_repr
from scripts.common.path_utils import get_paths, remove_if_exists


class Cubic(CEF):
    """
    Class defining the trivalent rare earth compound in crystal
    with cubic symmetry,
    its crystal field parameters and the eigenvalues and eigenfunctions
    of the CEF Hamiltonian, if it is already diagonalized.

    """

    def __init__(self, material: Material, llw_parameters: dict):
        """Initializes the Cubic object or read it from a file."""
        super().__init__(material=material)
        if self.material.rare_earth.name in ['Ce', 'Sm', 'Eu']:
            print(
                f"The element '{self.material.rare_earth.name}' "
                f"is not supported."
            )
            sys.exit(1)
        else:
            self.llw_parameters = llw_parameters

    @property
    def parameters(self):
        """CEF parameters"""
        parameters = super().parameters
        try:
            parameters['B40'] = (
                    self.llw_parameters['w'] * self.llw_parameters['x'] / F4
            )
            parameters['B44'] = 5 * parameters['B40']
            parameters['B60'] = (self.llw_parameters['w'] *
                                 (1 - abs(self.llw_parameters['x'])) /
                                 self.material.rare_earth.f_6)
            parameters['B64'] = -21 * parameters['B60']
        except KeyError:
            pass
        return parameters

    def __repr__(self):
        """Method returns string representation of the Cubic object."""
        return get_repr(self, 'material', 'llw_parameters')

    def get_one_dot(self):
        """Prints information about RE ion with specified parameters"""
        for key in 'wx':
            print(f'{key}:\t\t{self.llw_parameters[key]: 9.3f}')
        for i, energy in enumerate(self.get_energies()):
            print(f'E[{i + 1}]:\t{energy: 9.3f} meV')

    def get_file_name(self, data_name: str, parameters=None):
        """Returns file_name for data saving"""
        parameters = self.llw_parameters if parameters is None else parameters
        return get_paths(
            data_name=data_name,
            material=self.material,
            parameters=parameters,
        )

    @get_time_of_execution
    def save_peak_dat(self, number_of_intervals: int, choice=0):
        """
        Saves the dependence of transition energies
        or intensities on parameter x to file.

        """
        file_name = self.get_file_name(
            data_name='energies' if choice == 0 else 'intensities'
        )
        remove_if_exists(file_name)
        print(f'Saving of {"energies" if choice == 0 else "intensities"} '
              f'datafiles will take some time...')
        with OpenedFile(file_name, mode='a') as file:
            for x_parameter in linspace(-1, 1, number_of_intervals + 1):
                self.llw_parameters['x'] = x_parameter
                row = (
                    self.get_energies()
                    if choice == 0
                    else self.get_intensities()
                )
                write_row(file, (x_parameter, *row))

    @get_time_of_execution
    def save_spectra_with_one_temperature(
            self,
            gamma: float,
            temperature: float,
    ):
        """
        Saves inelastic neutron scattering spectra
        at specified temperature to file.

        """
        energies = linspace(-5, 30, 10001)
        spectrum = self.get_spectrum(
            energies=energies,
            width_dict={'gamma': gamma},
            temperature=temperature,
        )
        file_name = self.get_file_name(
            data_name='spectra',
            parameters={
                **self.llw_parameters,
                'gamma': gamma,
                'T': temperature,
            }
        )
        remove_if_exists(file_name)
        with OpenedFile(file_name, mode='a') as file:
            for index, energy in enumerate(energies):
                write_row(file, (energy, spectrum[index]))

    def save_spectra_with_many_temperatures(
            self,
            gamma: float,
            temperatures,
    ):
        """
        Saves inelastic neutron scattering spectra
        at several specified temperatures to file.

        """
        lines = {}
        parameters = {
            **self.llw_parameters,
            'gamma': gamma,
        }
        data = {
            'energies': [],
        }
        for t_number, temperature in enumerate(temperatures):
            data[temperature] = []
            parameters['T'] = temperature
            file_name = self.get_file_name(
                data_name='spectra',
                parameters=parameters,
            )
            with OpenedFile(file_name) as file:
                lines[temperature] = list(file)

            for index, line in enumerate(lines[temperature]):
                row = line.split('\t')
                if t_number == 0:
                    data['energies'].append(float(row[0]))
                data[temperature].append(float(row[1]))

        del parameters['T']
        file_name = self.get_file_name(
            data_name='spectra',
            parameters=parameters,
        )
        remove_if_exists(file_name)
        with OpenedFile(file_name, mode='a') as file:
            for index, _ in enumerate(data['energies']):
                write_row(file, row=[data[key][index] for key in data])
        return data

    @get_time_of_execution
    def save_susceptibility(self):
        """
        Saves temperature dependence of magnetic susceptibilities to file.

        """
        temperatures = linspace(0.1, 100.0, 300)
        common_file_name = self.get_file_name(
            data_name='susceptibilities',
        )
        chi_curie, chi_van_vleck, chi = self.get_chi_dependence(temperatures)
        for axis in ('z', 'x', 'total'):
            file_name = common_file_name.replace('.dat', f'_chi_{axis}.dat')
            remove_if_exists(file_name)
            with OpenedFile(file_name, mode='a') as file:
                row = ['T(Kelvin)']
                if axis in ('z', 'x'):
                    row += [
                        f'chi_curie_{axis}',
                        f'chi_van_vleck_{axis}',
                        f'chi_{axis}',
                    ]
                else:
                    row += [
                        'chi_total',
                        'inverse_chi',
                    ]
                write_row(file, row)
                for i, temperature in enumerate(temperatures):
                    row = [temperature]
                    if axis in ('z', 'x'):
                        row += [
                            chi_curie[axis][i],
                            chi_van_vleck[axis][i],
                            chi[axis][i],
                        ]
                    else:
                        row += [
                            chi['total'][i],
                            chi['inverse'][i],
                        ]
                    write_row(file, row)

    @get_time_of_execution
    def get_ratios(self, choice=0):
        """
        Saves the dependence of transition energies ratio on parameter x to file.

        """
        peak_data = 'energies' if choice == 0 else 'intensities'
        levels_number = 7
        parameters = {
            'w': self.llw_parameters['w'],
        }
        ratio_file_name = self.get_file_name(
            data_name=f'ratios_{peak_data}',
            parameters=parameters
        )
        peak_file_name = self.get_file_name(
            data_name=peak_data,
            parameters=parameters
        )
        remove_if_exists(ratio_file_name)
        with OpenedFile(ratio_file_name, mode='a') as ratio_file:
            with OpenedFile(peak_file_name) as peak_file:
                for line in peak_file:
                    line = line.rstrip('\n')
                    peak_row = [float(energy) for energy in line.split('\t')]
                    for _ in range(len(peak_row), levels_number):
                        peak_row.append(0)
                    ratios = [peak_row[0]]
                    for low in range(1, levels_number):
                        for high in range(low + 1, levels_number):
                            if peak_row[low] == 0:
                                ratios.append(0)
                            else:
                                ratios.append(peak_row[high] / peak_row[low])
                    write_row(ratio_file, ratios)

    def check_ratios(
            self,
            numbers,
            experimental_value: float,
            points,
            accuracy: float,
    ):
        """
        Checks the array of ratios
        if one of them is approximately equal to the given value.

        """
        ratios = numbers[1:]
        for index, ratio in enumerate(ratios):
            if abs(experimental_value - ratio) < accuracy:
                current = CrossPoint(
                    rare_earth=self.material.rare_earth.name,
                    w=self.llw_parameters['w'],
                    x=numbers[0],
                    ratio_name=get_ratios_names(0)[index],
                    difference=experimental_value - ratio,
                )
                if not points:
                    points.append(current)
                else:
                    previous = points[-1]
                    if (
                            current.rare_earth == previous.rare_earth and
                            current.ratio_name == previous.ratio_name and
                            abs(current.w - previous.w) < accuracy and
                            abs(current.x - previous.x) < accuracy
                    ):
                        current_x = ((current.x * previous.difference -
                                      previous.x * current.difference) /
                                     (previous.difference - current.difference)
                                     )
                        points[-1] = CrossPoint(
                            rare_earth=self.material.rare_earth.name,
                            w=self.llw_parameters['w'],
                            x=current_x,
                            difference=0,
                            ratio_name=get_ratios_names(0)[index],
                        )
                    else:
                        points.append(current)
        return points

    def find_cross(
            self,
            experimental_value: float,
            experimental_energy: float,
            accuracy=0.005,
    ):
        """Returns points of cross experimental and calculated curves,
        recalculated with correct value of W."""
        points = []
        for w_parameter in (
                abs(self.llw_parameters['w']),
                -abs(self.llw_parameters['w']),
        ):
            self.llw_parameters['w'] = w_parameter
            ratio_file_name = self.get_file_name(
                data_name='ratios_energies',
                parameters={'w': w_parameter},
            )
            with OpenedFile(ratio_file_name) as ratio_file:
                for line in ratio_file:
                    line = line.rstrip('\n')
                    numbers = [float(number) for number in line.split('\t')]
                    if any(
                            abs(experimental_value - value) < accuracy
                            for value in numbers[1:]
                    ):
                        points = self.check_ratios(
                            numbers=numbers,
                            points=points,
                            experimental_value=experimental_value,
                            accuracy=accuracy
                        )
        for index, point in enumerate(points):
            self.llw_parameters['x'] = point.x
            self.llw_parameters['w'] = point.w
            level = int(point.ratio_name[-2]) - 1
            old_energy = self.get_energies()[level]
            new_w = experimental_energy / old_energy
            new_w = -new_w if self.llw_parameters['w'] < 0 else new_w
            points[index] = CrossPoint(
                w=new_w,
                rare_earth=point.rare_earth,
                x=point.x,
                difference=point.difference,
                ratio_name=point.ratio_name,
            )
            self.llw_parameters = {'w': points[index].w, 'x': points[index].x}
        return points

    @get_time_of_execution
    def save_intensities(self):
        """
        Saves temperature dependence of magnetic susceptibilities to file.

        """
        temperatures = linspace(0, 200, 1001)
        file_name = self.get_file_name(
            data_name='intensities_on_temperature',
        )
        remove_if_exists(file_name)
        with OpenedFile(file_name, mode='a') as file:
            for _, temperature in enumerate(temperatures):
                peaks = self.get_peaks(temperature=temperature)
                intensities = [peak[1] for peak in peaks if peak[0] >= 0]
                row = [temperature] + intensities
                write_row(file, row)


if __name__ == '__main__':
    print(
        Cubic(
            material=Material(rare_earth='Pr', crystal='YNi2'),
            llw_parameters={'w': -0.505, 'x': -0.107},
        )
    )
