"""The module contains functions for printing or saving results."""
from numpy import linspace
from scripts.common import constants
from scripts.common.path_utils import get_paths, remove_if_exists
from scripts.common.utils import get_time_of_execution, value_to_write, OpenedFile
from scripts.common.named_tuples import CrossPoint, Material
from scripts.get_cef_object import CF


def get_one_dot(material: Material, parameters: dict):
    """Prints information about RE ion with specified parameters"""
    cef_object = CF(material).get_llw_cubic_object(parameters)
    energies = cef_object.get_energies()
    print(f'x:\t\t{parameters["x"]: 9.3f}')
    for i, energy in enumerate(energies):
        print(f'E[{i + 1}]:\t{energy: 9.3f} meV')


def load_data(material: Material):
    """Loads CF object from file"""
    file_name = get_paths(constants.PATH_TO_SAVED_OBJECTS, 'parameters',
                          format_name='.json',
                          material=material)
    return CF(material=material, file_name=file_name)


@get_time_of_execution
def save_parameters(material: Material, parameters: dict):
    """Saves CEF parameters to file."""
    cef_object = CF(material).get_llw_cubic_object(parameters)
    print(f'\nSaving CEF parameters')
    cef_object.save_to_file()


@get_time_of_execution
def save_energy_dat(material: Material, parameters: dict, number_of_intervals):
    """Saves the dependence of transition energies on parameter x to file."""
    x_space = linspace(-1, 1, number_of_intervals + 1)
    file_name = get_paths(constants.PATH_TO_ENERGY_DATAFILES, 'energy',
                          material=material, parameters=parameters)
    remove_if_exists(file_name)
    print(f'\nSaving energy datafiles\nIt will take some time...')
    with OpenedFile(file_name, mode='a') as file:
        for x_parameter in x_space:
            parameters['x'] = x_parameter
            cef_object = CF(material).get_llw_cubic_object(parameters)
            energies = cef_object.get_energies()
            file.write(value_to_write(x_parameter, '\t'))
            for level, energy in enumerate(energies):
                if level == len(energies) - 1:
                    file.write(value_to_write(energy, '\n'))
                else:
                    file.write(value_to_write(energy, '\t'))


@get_time_of_execution
def save_spectra_with_one_temperature(material: Material, parameters: dict):
    """Saves inelastic neutron scattering spectra at specified temperature to file."""
    energies = linspace(-5, 30, 10001)
    print('\nSaving neutron inelastic scattering spectra')
    cef_object = CF(material).get_llw_cubic_object(parameters)
    spectrum = cef_object.get_spectrum(energy=energies,
                                       width_dict={'gamma': parameters['gamma']},
                                       temperature=parameters['T'])
    file_name = get_paths(constants.PATH_TO_SPECTRA_DATAFILES, 'spectrum',
                          material=material, parameters=parameters)
    with OpenedFile(file_name, mode='a') as file:
        for index, energy in enumerate(energies):
            file.write(value_to_write(energy, '\t'))
            file.write(value_to_write(spectrum[index], '\n'))


def save_spectra_with_two_temperatures(material: Material, parameters: dict,
                                       temperature_1, temperature_2):
    """Saves inelastic neutron scattering spectra at two specified temperatures to file."""
    print('\nSaving neutron inelastic scattering spectra')
    energies = []
    intensities = {'diff': []}
    lines = {}
    for temperature in (temperature_1, temperature_2):
        parameters['T'] = temperature
        intensities[temperature] = []
        file_name = get_paths(constants.PATH_TO_SPECTRA_DATAFILES, 'spectrum',
                              material=material, parameters=parameters)
        with OpenedFile(file_name) as file:
            lines[temperature] = list(file)

    for index in range(len(lines[temperature_1])):
        line_1 = lines[temperature_1][index].split('\t')
        line_2 = lines[temperature_2][index].split('\t')
        energies.append(float(line_1[0]))
        intensities[temperature_1].append(float(line_1[1]))
        intensities[temperature_2].append(float(line_2[1]))
        intensities['diff'].append(float(line_1[1]) - float(line_2[1]))

    parameters['T'] = f'{temperature_1}-{temperature_2}'
    file_name = get_paths(constants.PATH_TO_SPECTRA_DATAFILES, 'spectrum',
                          material=material, parameters=parameters)
    with OpenedFile(file_name, mode='a') as file:
        for index, energy in enumerate(energies):
            for value in (energy, intensities[temperature_1][index],
                          intensities[temperature_2][index]):
                file.write(value_to_write(value, '\t'))
            file.write(value_to_write(intensities['diff'][index], '\n'))
    return energies, intensities


@get_time_of_execution
def save_susceptibility(material: Material, parameters: dict):
    """Saves temperature dependence of magnetic susceptibilities to file."""
    temperatures = linspace(0.1, 100.0, 300)
    print('\nSaving magnetic susceptibilities')
    common_file_name = get_paths(constants.PATH_TO_SUSCEPTIBILITY_DATAFILES,
                                 'susceptibility',
                                 material=material, parameters=parameters)
    cef_object = CF(material).get_llw_cubic_object(parameters)
    chi_curie, chi_van_vleck, chi = cef_object.get_chi_dependence(temperatures)
    for axis in ('z', 'x', 'total'):
        file_name = common_file_name.replace('.dat', f'_chi_{axis}.dat')
        remove_if_exists(file_name)
        with OpenedFile(file_name, mode='a') as file:
            if axis in ('z', 'x'):
                file.write(f'T(Kelvin)\tchi_curie_{axis}\tchi_van_vleck_{axis}\tchi_{axis}\n')
            else:
                file.write('T(Kelvin)\tchi_total\tinverse_chi\n')
            for i, temperature in enumerate(temperatures):
                file.write(value_to_write(temperature, '\t'))
                if axis in ('z', 'x'):
                    file.write(value_to_write(chi_curie[axis][i], '\t'))
                    file.write(value_to_write(chi_van_vleck[axis][i], '\t'))
                    file.write(value_to_write(chi[axis][i], '\n'))
                else:
                    file.write(value_to_write(chi['total'][i], '\t'))
                    file.write(value_to_write(chi['inverse'][i], '\n'))


@get_time_of_execution
def get_ratios(material: Material, parameters: dict):
    """Saves the dependence of transition energies ratio on parameter x to file."""
    levels_number = 7
    ratio_file_name = get_paths(constants.PATH_TO_RATIO_DATAFILES,
                                'ratio',
                                material=material, parameters=parameters)
    energy_file_name = get_paths(constants.PATH_TO_ENERGY_DATAFILES,
                                 'energy',
                                 material=material, parameters=parameters)
    remove_if_exists(ratio_file_name)
    with OpenedFile(ratio_file_name, mode='a') as ratio_file:
        with OpenedFile(energy_file_name) as energy_file:
            print('\n', 'Saving ratio datafiles',
                  'It will take some time...', sep='\n')
            for line in energy_file:
                line = line.rstrip('\n')
                energies = [float(energy) for energy in line.split('\t')]
                if len(energies) < levels_number:
                    for _ in range(len(energies), levels_number):
                        energies.append(0)
                ratios = [energies[0]]
                for low in range(1, levels_number):
                    for high in range(low + 1, levels_number):
                        if energies[low] == 0:
                            ratios.append(0)
                        else:
                            ratios.append(energies[high] / energies[low])

                for level, ratio in enumerate(ratios):
                    if level == len(ratios) - 1:
                        ratio_file.write(value_to_write(ratio, '\n'))
                    else:
                        ratio_file.write(value_to_write(ratio, '\t'))


def check_ratios(numbers, points, material: Material, parameters: dict):
    """Checks the array of ratios
    if one of them is approximately equal to the given value."""
    ratios = numbers[1:]
    for index, ratio in enumerate(ratios):
        if abs(parameters['value'] - ratio) < parameters['accuracy']:
            current = CrossPoint(rare_earth=material.rare_earth,
                                 w=parameters['w'],
                                 x=numbers[0],
                                 ratio_name=constants.RATIOS_NAMES[index],
                                 difference=parameters['value'] - ratio)
            if not points:
                points.append(current)
            else:
                previous = points[len(points) - 1]
                if (current.rare_earth == previous.rare_earth and
                        current.w == previous.w and
                        current.ratio_name == previous.ratio_name and
                        abs(current.x - previous.x) < 1e-3):
                    current_x = ((current.x * previous.difference -
                                  previous.x * current.difference) /
                                 (previous.difference - current.difference))
                    points[-1] = CrossPoint(rare_earth=material.rare_earth,
                                            w=parameters['w'],
                                            x=current_x,
                                            difference=0,
                                            ratio_name=constants.RATIOS_NAMES[index])
                else:
                    points.append(current)
    return points


def find_cross(experimental_value, material: Material, parameters: dict, accuracy=0.003):
    """Returns points of cross experimental and calculated curves."""
    points = []
    for w_parameter in (abs(parameters['w']), -abs(parameters['w'])):
        parameters['w'] = w_parameter
        ratio_file_name = get_paths(constants.PATH_TO_RATIO_DATAFILES,
                                    'ratio',
                                    material=material, parameters=parameters)
        with OpenedFile(ratio_file_name) as ratio_file:
            for line in ratio_file:
                line = line.rstrip('\n')
                numbers = [float(number) for number in line.split('\t')]
                if any(abs(experimental_value - value) < accuracy for value in numbers[1:]):
                    points = check_ratios(numbers, points, material,
                                          parameters={'value': experimental_value,
                                                      'w': w_parameter,
                                                      'accuracy': accuracy})
    return points


def recalculation(points, experimental_energy, material: Material):
    """Returns points of cross experimental and calculated curves with correct value of W."""
    recalculated = []
    for point in points:
        x_parameter = point.x
        old_w = point.w
        level = int(point.ratio_name[-2]) - 1
        cef = CF(material).get_llw_cubic_object({'w': old_w, 'x': x_parameter})
        old_energy = cef.get_energies()[level]
        new_w = experimental_energy / old_energy
        new_w = -new_w if old_w < 0 else new_w
        recalculated.append(CrossPoint(rare_earth=material.rare_earth,
                                       w=new_w,
                                       x=x_parameter,
                                       difference=point.difference,
                                       ratio_name=point.ratio_name))
    return recalculated


if __name__ == '__main__':
    for rare_earth in ('Tb', 'Tm', 'Er', 'Ho'):
        for w_value in (1, -1):
            MATERIAL = Material(crystal='YNi2', rare_earth=rare_earth)
            PARAMETERS = {'w': w_value}
            get_ratios(material=MATERIAL, parameters=PARAMETERS)
