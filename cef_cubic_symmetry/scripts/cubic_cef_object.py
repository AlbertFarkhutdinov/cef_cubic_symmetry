"""The module contains class for CEF with cubic symmetry."""
import sys
from numpy import linspace
from scripts.cef_object import CEF
from scripts.common.constants import CrossPoint, Material, RATIOS_NAMES
from scripts.common.tabular_information import F4
from scripts.common.utils import get_time_of_execution, OpenedFile, write_row
from scripts.common.path_utils import get_paths, remove_if_exists


class Cubic(CEF):
    """Class defining the trivalent rare earth compound in crystal with cubic symmetry,
    its crystal field parameters and the eigenvalues and eigenfunctions
    of the CEF Hamiltonian, if it is already diagonalized."""

    def __init__(self, material: Material, llw_parameters: dict):
        """Initializes the Cubic object or read it from a file."""
        super().__init__(material=material)
        if self.material.rare_earth.name in ['Ce', 'Sm', 'Eu']:
            print(f"The element '{self.material.rare_earth.name}' is not supported.")
            sys.exit(1)
        else:
            self.llw_parameters = llw_parameters

    @property
    def parameters(self):
        """CEF parameters"""
        parameters = super().parameters
        try:
            parameters['B40'] = self.llw_parameters['w'] * self.llw_parameters['x'] / F4
            parameters['B44'] = 5 * parameters['B40']
            parameters['B60'] = (self.llw_parameters['w'] *
                                 (1 - abs(self.llw_parameters['x'])) /
                                 self.material.rare_earth.f_6)
            parameters['B64'] = -21 * parameters['B60']
        except KeyError:
            pass
        return parameters

    def __repr__(self):
        """Method returns string representation of the object."""
        return ''.join(
            [
                f"{self.__class__.__name__}",
                f"(material={self.material!r}, ",
                f"llw_parameters={self.llw_parameters!r}')"
            ]
        )

    def get_one_dot(self):
        """Prints information about RE ion with specified parameters"""
        for key in 'wx':
            print(f'{key}:\t\t{self.llw_parameters[key]: 9.3f}')
        for i, energy in enumerate(self.get_energies()):
            print(f'E[{i + 1}]:\t{energy: 9.3f} meV')

    @get_time_of_execution
    def save_energy_dat(self, number_of_intervals):
        """Saves the dependence of transition energies on parameter x to file."""
        file_name = get_paths(
            data_name='energies',
            material=self.material,
            parameters=self.llw_parameters
        )
        remove_if_exists(file_name)
        print(f'Saving of energy datafiles will take some time...')
        with OpenedFile(file_name, mode='a') as file:
            for x_parameter in linspace(-1, 1, number_of_intervals + 1):
                self.llw_parameters['x'] = x_parameter
                energies = self.get_energies()
                write_row(file, (x_parameter, *energies))

    @get_time_of_execution
    def save_spectra_with_one_temperature(self, gamma, temperature):
        """Saves inelastic neutron scattering spectra at specified temperature to file."""
        energies = linspace(-5, 30, 10001)
        spectrum = self.get_spectrum(energies=energies,
                                     width_dict={'gamma': gamma},
                                     temperature=temperature)
        file_name = get_paths(
            data_name='spectra',
            material=self.material,
            parameters={
                **self.llw_parameters,
                'gamma': gamma,
                'T': temperature,
            }
        )
        with OpenedFile(file_name, mode='a') as file:
            for index, energy in enumerate(energies):
                write_row(file, (energy, spectrum[index]))

    def save_spectra_with_two_temperatures(self, gamma, temperature_1, temperature_2):
        """Saves inelastic neutron scattering spectra at two specified temperatures to file."""
        lines = {}
        parameters = {**self.llw_parameters, 'gamma': gamma}
        for temperature in (temperature_1, temperature_2):
            parameters['T'] = temperature
            file_name = get_paths(
                data_name='spectra',
                material=self.material,
                parameters=parameters,
            )
            with OpenedFile(file_name) as file:
                lines[temperature] = list(file)

        data = {
            'energies': [],
            temperature_1: [],
            temperature_2: [],
            'diff': [],
        }

        for index, line in enumerate(lines[temperature_1]):
            line_1 = line.split('\t')
            line_2 = lines[temperature_2][index].split('\t')
            data['energies'].append(float(line_1[0]))
            data[temperature_1].append(float(line_1[1]))
            data[temperature_2].append(float(line_2[1]))
            data['diff'].append(float(line_1[1]) - float(line_2[1]))

        parameters['T'] = f'{temperature_1}-{temperature_2}'
        file_name = get_paths(
            data_name='spectra',
            material=self.material,
            parameters=parameters,
        )
        with OpenedFile(file_name, mode='a') as file:
            for index, _ in enumerate(data['energies']):
                write_row(file, row=[data[key][index] for key in data])
        return data

    @get_time_of_execution
    def save_susceptibility(self):
        """Saves temperature dependence of magnetic susceptibilities to file."""
        temperatures = linspace(0.1, 100.0, 300)
        common_file_name = get_paths(
            data_name='susceptibilities',
            material=self.material,
            parameters=self.llw_parameters,
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
    def get_ratios(self):
        """Saves the dependence of transition energies ratio on parameter x to file."""
        levels_number = 7
        path_args = {
            'material': self.material,
            'parameters': {'w': self.llw_parameters['w']},
        }
        ratio_file_name = get_paths(data_name='ratios', **path_args)
        energy_file_name = get_paths(data_name='energies', **path_args)
        remove_if_exists(ratio_file_name)
        with OpenedFile(ratio_file_name, mode='a') as ratio_file:
            with OpenedFile(energy_file_name) as energy_file:
                for line in energy_file:
                    line = line.rstrip('\n')
                    energies = [float(energy) for energy in line.split('\t')]
                    for _ in range(len(energies), levels_number):
                        energies.append(0)
                    ratios = [energies[0]]
                    for low in range(1, levels_number):
                        for high in range(low + 1, levels_number):
                            if energies[low] == 0:
                                ratios.append(0)
                            else:
                                ratios.append(energies[high] / energies[low])
                    write_row(ratio_file, ratios)

    def check_ratios(self, numbers, experimental_value, points, accuracy):
        """Checks the array of ratios
        if one of them is approximately equal to the given value."""
        ratios = numbers[1:]
        for index, ratio in enumerate(ratios):
            if abs(experimental_value - ratio) < accuracy:
                current = CrossPoint(rare_earth=self.material.rare_earth.name,
                                     w=self.llw_parameters['w'],
                                     x=numbers[0],
                                     ratio_name=RATIOS_NAMES[index],
                                     difference=experimental_value - ratio)
                if not points:
                    points.append(current)
                else:
                    previous = points[-1]
                    if (
                            current.rare_earth == previous.rare_earth and
                            current.w == previous.w and
                            current.ratio_name == previous.ratio_name and
                            abs(current.x - previous.x) < 1e-3
                    ):
                        current_x = ((current.x * previous.difference -
                                      previous.x * current.difference) /
                                     (previous.difference - current.difference)
                                     )
                        points[-1] = CrossPoint(rare_earth=self.material.rare_earth.name,
                                                w=self.llw_parameters['w'],
                                                x=current_x,
                                                difference=0,
                                                ratio_name=RATIOS_NAMES[index])
                    else:
                        points.append(current)
        return points

    def find_cross(self, experimental_value, experimental_energy, accuracy=0.003):
        """Returns points of cross experimental and calculated curves,
        recalculated with correct value of W."""
        points = []
        for w_parameter in (abs(self.llw_parameters['w']), -abs(self.llw_parameters['w'])):
            self.llw_parameters['w'] = w_parameter
            ratio_file_name = get_paths(
                data_name='ratios',
                material=self.material,
                parameters={'w': w_parameter},
            )
            with OpenedFile(ratio_file_name) as ratio_file:
                for line in ratio_file:
                    line = line.rstrip('\n')
                    numbers = [float(number) for number in line.split('\t')]
                    if any(abs(experimental_value - value) < accuracy for value in numbers[1:]):
                        points = self.check_ratios(
                            numbers=numbers,
                            points=points,
                            experimental_value=experimental_value,
                            accuracy=accuracy
                        )
        for point in points:
            self.llw_parameters['x'] = point.x
            level = int(point.ratio_name[-2]) - 1
            old_energy = self.get_energies()[level]
            new_w = experimental_energy / old_energy
            new_w = -new_w if self.llw_parameters['w'] < 0 else new_w
            point = CrossPoint(w=new_w,
                               rare_earth=point.rare_earth,
                               x=point.x,
                               difference=point.difference,
                               ratio_name=point.ratio_name)
            self.llw_parameters = {'w': point.w, 'x': point.x}
        return points


if __name__ == '__main__':
    MATERIAL = Material(crystal='YNi2', rare_earth='Ce')
    LLW_PARAMETERS = {'w': 1, 'x': 0.1}
    CUBIC_SAMPLE = Cubic(
        material=MATERIAL,
        llw_parameters=LLW_PARAMETERS
    )
    print(dir(CUBIC_SAMPLE))
    print(repr(CUBIC_SAMPLE))
    print(CUBIC_SAMPLE)
    CUBIC_SAMPLE.get_one_dot()
    CUBIC_SAMPLE.save_to_file()
