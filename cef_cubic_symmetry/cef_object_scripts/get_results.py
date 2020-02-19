"""The module contains functions for printing or saving results."""
from collections import namedtuple
import numpy
from cef_object_scripts import common
from cef_object_scripts.tabular_information import F4
from cef_object_scripts.get_cef_object import CF

CrossPoint = namedtuple('CrossPoint', ['rare_earth', 'w', 'x', 'ratio_name', 'difference'])


def get_object_with_parameters(crystal, rare_earth, w_parameter, x_parameter):
    """Returns CF object for crystals with cubic symmetry"""
    cef_object = CF(name=f'{crystal}: {rare_earth}3+', rare_earth=rare_earth)
    if rare_earth in ['Ce', 'Sm', 'Eu']:
        print('This element is not supported. The default object will be returned.')
    else:
        cef_object.parameters['B40'] = w_parameter * x_parameter / F4
        cef_object.parameters['B44'] = 5 * cef_object.parameters['B40']
        cef_object.parameters['B60'] = (w_parameter * (1 - abs(x_parameter)) /
                                        cef_object.rare_earth.f_6)
        cef_object.parameters['B64'] = -21 * cef_object.parameters['B60']
    return cef_object


def get_one_dot(crystal, rare_earth, w_parameter, x_parameter):
    """Prints information about RE ion with specified parameters"""
    cef_object = get_object_with_parameters(crystal, rare_earth, w_parameter, x_parameter)
    energies = cef_object.get_energies()
    print(f'x:\t\t{x_parameter: 9.3f}')
    for i, energy in enumerate(energies):
        print(f'E[{i + 1}]:\t{energy: 9.3f} meV')


def load_data(crystal, rare_earth, w_parameter, x_parameter):
    """Loads CF object from file"""
    file_name = common.get_paths(common.PATH_TO_SAVED_OBJECTS,
                                 'parameters', 'cfg', crystal, rare_earth, w_parameter, x_parameter)
    return CF(name=f'{crystal}: {rare_earth}3+', rare_earth=rare_earth, par_file=file_name)


@common.get_time_of_execution
def save_energy_dat(crystal, rare_earth, w_parameter, number_of_intervals):
    """Saves the dependence of transition energies on parameter x to file."""
    x_space = numpy.linspace(-1, 1, number_of_intervals + 1)
    file_name = common.get_paths(common.PATH_TO_ENERGY_DATAFILES,
                                 'energy', 'dat', crystal, rare_earth, w_parameter)
    common.remove_if_exists(file_name)
    my_file = open(file_name, 'a', encoding='utf-8')
    print(f'\nSaving energy datafiles\nSaving file "{file_name}"...\nIt will take some time...')
    for x_parameter in x_space:
        cef_object = get_object_with_parameters(crystal, rare_earth, w_parameter, x_parameter)
        energies = cef_object.get_energies()
        my_file.write(common.value_to_write(x_parameter, '\t'))
        for level, energy in enumerate(energies):
            if level == len(energies) - 1:
                my_file.write(common.value_to_write(energy, '\n'))
            else:
                my_file.write(common.value_to_write(energy, '\t'))
    my_file.close()
    print(f'File "{file_name}" is saved.')


@common.get_time_of_execution
def save_parameters(crystal, rare_earth, w_parameter, x_parameter):
    """Saves CEF parameters to file."""
    cef_object = get_object_with_parameters(crystal, rare_earth, w_parameter, x_parameter)
    file_name = common.get_paths(common.PATH_TO_SAVED_OBJECTS,
                                 'parameters', 'cfg', crystal, rare_earth, w_parameter, x_parameter)
    print(f'\nSaving CEF parameters\nSaving file "{file_name}"...')
    cef_object.save_to_par_file(file_name)
    print(f'File "{file_name}" is saved.')


@common.get_time_of_execution
def save_spectra_with_one_temperature(crystal, rare_earth, w_parameter, x_parameter,
                                      temperature, gamma):
    """Saves inelastic neutron scattering spectra at specified temperature to file."""
    energies = numpy.linspace(-5, 30, 10001)
    print('\nSaving neutron inelastic scattering spectra')
    cef_object = get_object_with_parameters(crystal, rare_earth, w_parameter, x_parameter)
    spectrum = cef_object.spectrum(energy=energies, gamma=gamma, temperature=temperature)
    file_name = common.get_paths(common.PATH_TO_SPECTRA_DATAFILES, 'spectrum', 'dat',
                                 crystal, rare_earth, w_parameter, x_parameter, temperature)
    common.check_path(file_name)
    print(f'Saving file {file_name}...')
    numpy.savetxt(file_name, common.create_table(energies, spectrum), delimiter='\t')
    print(f'File "{file_name}" is saved')


# @common.get_time_of_execution
def save_spectra_with_two_temperatures(crystal, rare_earth, w_parameter, x_parameter,
                                       temperature_1, temperature_2):
    """Saves inelastic neutron scattering spectra at two specified temperatures to file."""
    print('\nSaving neutron inelastic scattering spectra')
    energies = []
    intensities = {'diff': []}
    lines = {}
    for temperature in (temperature_1, temperature_2):
        intensities[temperature] = []
        file_name = common.get_paths(common.PATH_TO_SPECTRA_DATAFILES, 'spectrum', 'dat',
                                     crystal, rare_earth, w_parameter, x_parameter, temperature)
        common.check_path(file_name)
        with open(file_name, 'r', encoding='utf-8') as file_name:
            lines[temperature] = list(file_name)

    for index in range(len(lines[temperature_1])):
        line_1 = lines[temperature_1][index].split('\t')
        line_2 = lines[temperature_2][index].split('\t')
        energies.append(float(line_1[0]))
        intensities[temperature_1].append(float(line_1[1]))
        intensities[temperature_2].append(float(line_2[1]))
        intensities['diff'].append(float(line_1[1]) - float(line_2[1]))

    file_name = common.get_paths(common.PATH_TO_SPECTRA_DATAFILES, 'spectrum', 'dat',
                                 crystal, rare_earth, w_parameter, x_parameter,
                                 f'{temperature_1}-{temperature_2}')
    common.remove_if_exists(file_name)
    print(f'Saving file {file_name}...')
    file = open(file_name, 'a', encoding='utf-8')
    for index, energy in enumerate(energies):
        for value in (energy, intensities[temperature_1][index],
                      intensities[temperature_2][index]):
            file.write(common.value_to_write(value, '\t'))
        file.write(common.value_to_write(intensities['diff'][index], '\n'))
    print(f'File "{file_name}" is saved')
    return energies, intensities


@common.get_time_of_execution
def save_susceptibility(crystal, rare_earth, w_parameter, x_parameter):
    """Saves temperature dependence of magnetic susceptibilities to file."""
    temperatures = numpy.linspace(0.1, 100.0, 300)
    print('\nSaving magnetic susceptibilities')
    common_file_name = common.get_paths(common.PATH_TO_SUSCEPTIBILITY_DATAFILES,
                                        'susceptibility', 'dat',
                                        crystal, rare_earth, w_parameter, x_parameter)
    cef_object = get_object_with_parameters(crystal, rare_earth, w_parameter, x_parameter)
    chi_s = cef_object.chi_s(temperatures)
    for axis in ('z', 'x', 'total'):
        file_name = common_file_name.replace('.dat', f'_chi_{axis}.dat')
        common.remove_if_exists(file_name)
        my_file = open(file_name, 'a', encoding='utf-8')
        print(f'Saving file "{file_name}"...')
        if axis in ('z', 'x'):
            my_file.write(f'T(Kelvin)\tchi_curie_{axis}\tchi_van_vleck_{axis}\tchi_{axis}\n')
        else:
            my_file.write('T(Kelvin)\tchi_total\tinverse_chi\n')
        for i, temperature in enumerate(temperatures):
            my_file.write(common.value_to_write(temperature, '\t'))
            if axis in ('z', 'x'):
                my_file.write(common.value_to_write(chi_s[f'chi_curie_{axis}'][i], '\t'))
                my_file.write(common.value_to_write(chi_s[f'chi_van_vleck_{axis}'][i], '\t'))
                my_file.write(common.value_to_write(chi_s[f'chi_{axis}'][i], '\n'))
            else:
                my_file.write(common.value_to_write(chi_s['chi'][i], '\t'))
                my_file.write(common.value_to_write(chi_s['inverse_chi'][i], '\n'))
        my_file.close()
        print(f'File "{file_name}" is saved')


@common.get_time_of_execution
def get_ratios(crystal, rare_earth, w_parameter):
    """Saves the dependence of transition energies ratio on parameter x to file."""
    ratio_file_name = common.get_paths(common.PATH_TO_RATIO_DATAFILES,
                                       'ratio', 'dat', crystal, rare_earth, w_parameter)
    energy_file_name = common.get_paths(common.PATH_TO_ENERGY_DATAFILES,
                                        'energy', 'dat', crystal, rare_earth, w_parameter)
    common.remove_if_exists(ratio_file_name)
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
                ratio_file.write(common.value_to_write(ratio, '\n'))
            else:
                ratio_file.write(common.value_to_write(ratio, '\t'))

    ratio_file.close()
    energy_file.close()
    print(f'File "{ratio_file_name}" is saved.')


def find_cross(experimental_value, crystal, rare_earth, w_parameter):
    """Returns points of cross experimental and calculated curves."""
    accuracy = 0.005
    points = []
    for _w in (abs(w_parameter), -abs(w_parameter)):
        ratio_file_name = common.get_paths(common.PATH_TO_RATIO_DATAFILES,
                                           'ratio', 'dat', crystal, rare_earth, _w)
        ratio_file = open(ratio_file_name, 'r', encoding='utf-8')
        for line in ratio_file:
            line = line.rstrip('\n')
            numbers = [float(number) for number in line.split('\t')]
            x_parameter = numbers[0]
            ratios = numbers[1:]
            if any(abs(experimental_value - value) < accuracy for value in ratios):
                for index, ratio in enumerate(ratios):
                    if abs(experimental_value - ratio) < accuracy:
                        current = CrossPoint(rare_earth=rare_earth,
                                             w=_w,
                                             x=x_parameter,
                                             ratio_name=common.RATIOS_NAMES[index],
                                             difference=experimental_value - ratio)
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
                                points[-1] = CrossPoint(rare_earth=rare_earth,
                                                        w=_w,
                                                        x=current_x,
                                                        difference=0,
                                                        ratio_name=common.RATIOS_NAMES[index])
                            else:
                                points.append(current)

    return points


def recalculation(points, experimental_energy, crystal, rare_earth):
    """Returns points of cross experimental and calculated curves with correct value of W."""
    recalculated = []
    for point in points:
        x_parameter = point.x
        old_w = point.w
        level = int(point.ratio_name[-1]) - 1
        cef = get_object_with_parameters(crystal, rare_earth, old_w, x_parameter)
        old_energy = cef.get_energies()[level]
        new_w = experimental_energy / old_energy
        new_w = -new_w if old_w < 0 else new_w
        recalculated.append(CrossPoint(rare_earth=rare_earth,
                                       w=new_w,
                                       x=x_parameter,
                                       difference=point.difference,
                                       ratio_name=point.ratio_name))
    return recalculated


if __name__ == '__main__':
    CRYSTAL_NAME = 'YNi2'
    RARE_EARTH_AME = 'Tm'
    W = 1
    CROSSES = find_cross(4.59, CRYSTAL_NAME, RARE_EARTH_AME, W)
    RECALCULATED_CROSSES = recalculation(CROSSES, 0.45606, CRYSTAL_NAME, RARE_EARTH_AME)
    print(*CROSSES, sep='\n')
    print(*RECALCULATED_CROSSES, sep='\n')
