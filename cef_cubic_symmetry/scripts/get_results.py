"""The module contains functions for printing or saving results."""
from collections import namedtuple
from numpy import linspace, savetxt
from .common import constants
from .common.tabular_information import F4
from .common.path_utils import get_paths, check_path, remove_if_exists
from .common.utils import create_table, get_time_of_execution, value_to_write
from .get_cef_object import CF

CrossPoint = namedtuple('CrossPoint', ['rare_earth', 'w', 'x', 'ratio_name', 'difference'])


def get_object_with_parameters(material: dict, parameters: dict):
    """Returns CF object for crystals with cubic symmetry"""
    cef_object = CF(name=f'{material["crystal"]}: {material["rare_earth"]}3+',
                    rare_earth=material["rare_earth"])
    if material["rare_earth"] in ['Ce', 'Sm', 'Eu']:
        print('This element is not supported. The default object will be returned.')
    else:
        cef_object.parameters['B40'] = parameters['w'] * parameters['x'] / F4
        cef_object.parameters['B44'] = 5 * cef_object.parameters['B40']
        cef_object.parameters['B60'] = (parameters['w'] * (1 - abs(parameters['x'])) /
                                        cef_object.rare_earth.f_6)
        cef_object.parameters['B64'] = -21 * cef_object.parameters['B60']
    return cef_object


def get_one_dot(material: dict, parameters: dict):
    """Prints information about RE ion with specified parameters"""
    cef_object = get_object_with_parameters(material, parameters)
    energies = cef_object.get_energies()
    print(f'x:\t\t{parameters["x"]: 9.3f}')
    for i, energy in enumerate(energies):
        print(f'E[{i + 1}]:\t{energy: 9.3f} meV')


def load_data(material: dict, parameters: dict):
    """Loads CF object from file"""
    file_name = get_paths(constants.PATH_TO_SAVED_OBJECTS, 'parameters', 'cfg',
                          material=material, parameters=parameters)
    return CF(name=f'{material["crystal"]}: {material["rare_earth"]}3+',
              rare_earth=material["rare_earth"], par_file=file_name)


@get_time_of_execution
def save_energy_dat(material: dict, parameters: dict, number_of_intervals):
    """Saves the dependence of transition energies on parameter x to file."""
    x_space = linspace(-1, 1, number_of_intervals + 1)
    file_name = get_paths(constants.PATH_TO_ENERGY_DATAFILES, 'energy', 'dat',
                          material=material, parameters=parameters)
    remove_if_exists(file_name)
    my_file = open(file_name, 'a', encoding='utf-8')
    print(f'\nSaving energy datafiles\nSaving file "{file_name}"...\nIt will take some time...')
    for x_parameter in x_space:
        parameters['x'] = x_parameter
        cef_object = get_object_with_parameters(material, parameters)
        energies = cef_object.get_energies()
        my_file.write(value_to_write(x_parameter, '\t'))
        for level, energy in enumerate(energies):
            if level == len(energies) - 1:
                my_file.write(value_to_write(energy, '\n'))
            else:
                my_file.write(value_to_write(energy, '\t'))
    my_file.close()
    print(f'File "{file_name}" is saved.')


@get_time_of_execution
def save_parameters(material: dict, parameters: dict):
    """Saves CEF parameters to file."""
    cef_object = get_object_with_parameters(material, parameters)
    file_name = get_paths(constants.PATH_TO_SAVED_OBJECTS, 'parameters', 'cfg',
                          material=material, parameters=parameters)
    print(f'\nSaving CEF parameters\nSaving file "{file_name}"...')
    cef_object.save_to_par_file(file_name)
    print(f'File "{file_name}" is saved.')


@get_time_of_execution
def save_spectra_with_one_temperature(material: dict, parameters: dict):
    """Saves inelastic neutron scattering spectra at specified temperature to file."""
    energies = linspace(-5, 30, 10001)
    print('\nSaving neutron inelastic scattering spectra')
    cef_object = get_object_with_parameters(material, parameters)
    spectrum = cef_object.spectrum(energy=energies,
                                   gamma=parameters['gamma'],
                                   temperature=parameters['T'])
    file_name = get_paths(constants.PATH_TO_SPECTRA_DATAFILES, 'spectrum', 'dat',
                          material=material, parameters=parameters)
    check_path(file_name)
    print(f'Saving file {file_name}...')
    savetxt(file_name, create_table(energies, spectrum), delimiter='\t')
    print(f'File "{file_name}" is saved')


# @get_time_of_execution
def save_spectra_with_two_temperatures(material: dict, parameters: dict,
                                       temperature_1, temperature_2):
    """Saves inelastic neutron scattering spectra at two specified temperatures to file."""
    print('\nSaving neutron inelastic scattering spectra')
    energies = []
    intensities = {'diff': []}
    lines = {}
    for temperature in (temperature_1, temperature_2):
        parameters['T'] = temperature
        intensities[temperature] = []
        file_name = get_paths(constants.PATH_TO_SPECTRA_DATAFILES, 'spectrum', 'dat',
                              material=material, parameters=parameters)
        check_path(file_name)
        with open(file_name, 'r', encoding='utf-8') as file_name:
            lines[temperature] = list(file_name)

    for index in range(len(lines[temperature_1])):
        line_1 = lines[temperature_1][index].split('\t')
        line_2 = lines[temperature_2][index].split('\t')
        energies.append(float(line_1[0]))
        intensities[temperature_1].append(float(line_1[1]))
        intensities[temperature_2].append(float(line_2[1]))
        intensities['diff'].append(float(line_1[1]) - float(line_2[1]))

    parameters['T'] = f'{temperature_1}-{temperature_2}'
    file_name = get_paths(constants.PATH_TO_SPECTRA_DATAFILES, 'spectrum', 'dat',
                          material=material, parameters=parameters)
    remove_if_exists(file_name)
    print(f'Saving file {file_name}...')
    file = open(file_name, 'a', encoding='utf-8')
    for index, energy in enumerate(energies):
        for value in (energy, intensities[temperature_1][index],
                      intensities[temperature_2][index]):
            file.write(value_to_write(value, '\t'))
        file.write(value_to_write(intensities['diff'][index], '\n'))
    print(f'File "{file_name}" is saved')
    return energies, intensities


@get_time_of_execution
def save_susceptibility(material: dict, parameters: dict):
    """Saves temperature dependence of magnetic susceptibilities to file."""
    temperatures = linspace(0.1, 100.0, 300)
    print('\nSaving magnetic susceptibilities')
    common_file_name = get_paths(constants.PATH_TO_SUSCEPTIBILITY_DATAFILES,
                                 'susceptibility', 'dat',
                                 material=material, parameters=parameters)
    cef_object = get_object_with_parameters(material, parameters)
    chi_curie, chi_van_vleck, chi = cef_object.chi_s(temperatures)
    for axis in ('z', 'x', 'total'):
        file_name = common_file_name.replace('.dat', f'_chi_{axis}.dat')
        remove_if_exists(file_name)
        my_file = open(file_name, 'a', encoding='utf-8')
        print(f'Saving file "{file_name}"...')
        if axis in ('z', 'x'):
            my_file.write(f'T(Kelvin)\tchi_curie_{axis}\tchi_van_vleck_{axis}\tchi_{axis}\n')
        else:
            my_file.write('T(Kelvin)\tchi_total\tinverse_chi\n')
        for i, temperature in enumerate(temperatures):
            my_file.write(value_to_write(temperature, '\t'))
            if axis in ('z', 'x'):
                my_file.write(value_to_write(chi_curie[axis][i], '\t'))
                my_file.write(value_to_write(chi_van_vleck[axis][i], '\t'))
                my_file.write(value_to_write(chi[axis][i], '\n'))
            else:
                my_file.write(value_to_write(chi['total'][i], '\t'))
                my_file.write(value_to_write(chi['inverse'][i], '\n'))
        my_file.close()
        print(f'File "{file_name}" is saved')


@get_time_of_execution
def get_ratios(material: dict, parameters: dict):
    """Saves the dependence of transition energies ratio on parameter x to file."""
    ratio_file_name = get_paths(constants.PATH_TO_RATIO_DATAFILES,
                                'ratio', 'dat',
                                material=material, parameters=parameters)
    energy_file_name = get_paths(constants.PATH_TO_ENERGY_DATAFILES,
                                 'energy', 'dat',
                                 material=material, parameters=parameters)
    remove_if_exists(ratio_file_name)
    ratio_file = open(ratio_file_name, 'a', encoding='utf-8')
    energy_file = open(energy_file_name, 'r', encoding='utf-8')
    print('\n', 'Saving ratio datafiles',
          f'Saving file "{ratio_file_name}"...',
          'It will take some time...', sep='\n')
    for line in energy_file:
        line = line.rstrip('\n')
        energies = [float(energy) for energy in line.split('\t')]
        if len(energies) < 6:
            for _ in range(len(energies), 6):
                energies.append(0)
        ratios = [energies[0]]
        for low in range(1, 6):
            for high in range(low + 1, 6):
                if energies[low] == 0:
                    ratios.append(0)
                else:
                    ratios.append(energies[high] / energies[low])

        for level, ratio in enumerate(ratios):
            if level == len(ratios) - 1:
                ratio_file.write(value_to_write(ratio, '\n'))
            else:
                ratio_file.write(value_to_write(ratio, '\t'))

    ratio_file.close()
    energy_file.close()
    print(f'File "{ratio_file_name}" is saved.')


def check_ratios(numbers, points, material: dict, parameters: dict):
    """Checks the array of ratios
    if one of them is approximately equal to the given value."""
    ratios = numbers[1:]
    for index, ratio in enumerate(ratios):
        if abs(parameters['value'] - ratio) < parameters['accuracy']:
            current = CrossPoint(rare_earth=material['rare_earth'],
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
                    points[-1] = CrossPoint(rare_earth=material['rare_earth'],
                                            w=parameters['w'],
                                            x=current_x,
                                            difference=0,
                                            ratio_name=constants.RATIOS_NAMES[index])
                else:
                    points.append(current)
    return points


def find_cross(experimental_value, material: dict, parameters: dict, accuracy=0.005):
    """Returns points of cross experimental and calculated curves."""
    points = []
    for w_parameter in (abs(parameters['w']), -abs(parameters['w'])):
        parameters['w'] = w_parameter
        ratio_file_name = get_paths(constants.PATH_TO_RATIO_DATAFILES,
                                    'ratio', 'dat',
                                    material=material, parameters=parameters)
        ratio_file = open(ratio_file_name, 'r', encoding='utf-8')
        for line in ratio_file:
            line = line.rstrip('\n')
            numbers = [float(number) for number in line.split('\t')]
            if any(abs(experimental_value - value) < accuracy for value in numbers[1:]):
                points = check_ratios(numbers, points, material,
                                      parameters={'value': experimental_value,
                                                  'w': w_parameter,
                                                  'accuracy': accuracy})
    return points


def recalculation(points, experimental_energy, material: dict):
    """Returns points of cross experimental and calculated curves with correct value of W."""
    recalculated = []
    for point in points:
        x_parameter = point.x
        old_w = point.w
        level = int(point.ratio_name[-1]) - 1
        cef = get_object_with_parameters(material, {'w': old_w, 'x': x_parameter})
        old_energy = cef.get_energies()[level]
        new_w = experimental_energy / old_energy
        new_w = -new_w if old_w < 0 else new_w
        recalculated.append(CrossPoint(rare_earth=material['rare_earth'],
                                       w=new_w,
                                       x=x_parameter,
                                       difference=point.difference,
                                       ratio_name=point.ratio_name))
    return recalculated


if __name__ == '__main__':
    MATERIAL = {'crystal': 'YNi2', 'rare_earth': 'Tm'}
    PARAMETERS = {'w': 1}
    CROSSES = find_cross(4.59, MATERIAL, PARAMETERS)
    RECALCULATED_CROSSES = recalculation(CROSSES, 0.45606, MATERIAL)
    print(*CROSSES, sep='\n')
    print(*RECALCULATED_CROSSES, sep='\n')
