import numpy
from cef_object_scripts import common
from cef_object_scripts.tabular_information import f_4
from cef_object_scripts.get_cef_object_new import CF
from collections import namedtuple

cross_point = namedtuple('cross_point', ['rare_earth', 'w', 'x', 'ratio_name', 'difference'])


def get_object_with_parameters(crystal, rare_earth, w, x):
    cef_object = CF(name=f'{crystal}: {rare_earth}3+', rare_earth=rare_earth)
    if rare_earth in ['Ce', 'Sm', 'Eu']:
        print('This element is not supported. The default object will be returned.')
    else:
        cef_object.parameters['B40'] = w * x / f_4
        cef_object.parameters['B44'] = 5 * cef_object.parameters['B40']
        cef_object.parameters['B60'] = w * (1 - abs(x)) / cef_object.rare_earth.f_6
        cef_object.parameters['B64'] = -21 * cef_object.parameters['B60']
    return cef_object


def get_one_dot(crystal, rare_earth, w, x):
    cef_object = get_object_with_parameters(crystal, rare_earth, w, x)
    energies = cef_object.get_energies()
    print(f'x:\t\t{x: 9.3f}')
    for i in range(len(energies)):
        print(f'E[{i + 1}]:\t{energies[i]: 9.3f} meV')


def load_data(crystal, rare_earth, w, x):
    file_name = common.get_paths(common.path_to_saved_objects, 'parameters', 'cfg', crystal, rare_earth, w, x)
    return CF(name=f'{crystal}: {rare_earth}3+', rare_earth=rare_earth, par_file=file_name)


@common.get_time_of_execution
def save_energy_dat(crystal, rare_earth, w, number_of_intervals):
    x_space = numpy.linspace(-1, 1, number_of_intervals + 1)
    file_name = common.get_paths(common.path_to_energy_datafiles, 'energy', 'dat', crystal, rare_earth, w)
    common.remove_if_exists(file_name)
    my_file = open(file_name, 'a', encoding='utf-8')
    print(f'\nSaving energy datafiles\nSaving file "{file_name}"...\nIt will take some time...')
    for x in x_space:
        cef_object = get_object_with_parameters(crystal, rare_earth, w, x)
        energies = cef_object.get_energies()
        my_file.write(common.value_to_write(x, '\t'))
        for level in range(len(energies)):
            if level == len(energies) - 1:
                my_file.write(common.value_to_write(energies[level], '\n'))
            else:
                my_file.write(common.value_to_write(energies[level], '\t'))
    my_file.close()
    print(f'File "{file_name}" is saved.')


@common.get_time_of_execution
def save_parameters(crystal, rare_earth, w, x):
    cef_object = get_object_with_parameters(crystal, rare_earth, w, x)
    file_name = common.get_paths(common.path_to_saved_objects, 'parameters', 'cfg', crystal, rare_earth, w, x)
    print(f'\nSaving CEF parameters\nSaving file "{file_name}"...')
    cef_object.save_to_par_file(file_name)
    print(f'File "{file_name}" is saved.')


@common.get_time_of_execution
def save_spectra_with_one_temperature(crystal, rare_earth, w, x, temperature, gamma):
    energies = numpy.linspace(-5, 30, 10001)
    print('\nSaving neutron inelastic scattering spectra')
    cef_object = get_object_with_parameters(crystal, rare_earth, w, x)
    spectrum = cef_object.spectrum(energy=energies, gamma=gamma, temperature=temperature)
    file_name = common.get_paths(common.path_to_spectra_datafiles, 'spectrum', 'dat',
                                 crystal, rare_earth, w, x, temperature)
    common.check_path(file_name)
    print(f'Saving file {file_name}...')
    numpy.savetxt(file_name, common.create_table(energies, spectrum), delimiter='\t')
    print(f'File "{file_name}" is saved')


# @common.get_time_of_execution
def save_spectra_with_two_temperatures(crystal, rare_earth, w, x, temperature_1, temperature_2):
    print('\nSaving neutron inelastic scattering spectra')
    energies = []
    intensities = {'diff': []}
    lines = {}
    for temperature in (temperature_1, temperature_2):
        intensities[temperature] = []
        file_name = common.get_paths(common.path_to_spectra_datafiles, 'spectrum', 'dat',
                                     crystal, rare_earth, w, x, temperature)
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

    file_name = common.get_paths(common.path_to_spectra_datafiles, 'spectrum', 'dat',
                                 crystal, rare_earth, w, x, f'{temperature_1}-{temperature_2}')
    common.remove_if_exists(file_name)
    print(f'Saving file {file_name}...')
    file = open(file_name, 'a', encoding='utf-8')
    for index in range(len(energies)):
        for value in (energies[index], intensities[temperature_1][index], intensities[temperature_2][index]):
            file.write(common.value_to_write(value, '\t'))
        file.write(common.value_to_write(intensities['diff'][index], '\n'))
    print(f'File "{file_name}" is saved')
    return energies, intensities


@common.get_time_of_execution
def save_susceptibility(crystal, rare_earth, w, x):
    temperatures = numpy.linspace(0.1, 100.0, 300)
    print('\nSaving magnetic susceptibilities')
    common_file_name = common.get_paths(common.path_to_susceptibility_datafiles, 'susceptibility', 'dat',
                                        crystal, rare_earth, w, x)
    cef_object = get_object_with_parameters(crystal, rare_earth, w, x)
    chi_s = cef_object.chi_s(temperatures)
    for axis in ['z', 'x', 'total']:
        file_name = common_file_name.replace('.dat', f'_chi_{axis}.dat')
        common.remove_if_exists(file_name)
        my_file = open(file_name, 'a', encoding='utf-8')
        print(f'Saving file "{file_name}"...')
        if axis == 'z' or axis == 'x':
            my_file.write(f'T(Kelvin)\tchi_curie_{axis}\tchi_van_vleck_{axis}\tchi_{axis}\n')
        else:
            my_file.write('T(Kelvin)\tchi_total\tinverse_chi\n')
        for i in range(len(temperatures)):
            my_file.write(common.value_to_write(temperatures[i], '\t'))
            if axis == 'z' or axis == 'x':
                my_file.write(common.value_to_write(chi_s[f'chi_curie_{axis}'][i], '\t'))
                my_file.write(common.value_to_write(chi_s[f'chi_van_vleck_{axis}'][i], '\t'))
                my_file.write(common.value_to_write(chi_s[f'chi_{axis}'][i], '\n'))
            else:
                my_file.write(common.value_to_write(chi_s['chi'][i], '\t'))
                my_file.write(common.value_to_write(chi_s['inverse_chi'][i], '\n'))
        my_file.close()
        print(f'File "{file_name}" is saved')


@common.get_time_of_execution
def get_ratios(crystal, rare_earth, w):
    ratio_file_name = common.get_paths(common.path_to_ratio_datafiles, 'ratio', 'dat', crystal, rare_earth, w)
    energy_file_name = common.get_paths(common.path_to_energy_datafiles, 'energy', 'dat', crystal, rare_earth, w)
    common.remove_if_exists(ratio_file_name)
    ratio_file = open(ratio_file_name, 'a', encoding='utf-8')
    energy_file = open(energy_file_name, 'r', encoding='utf-8')
    print(f'\nSaving ratio datafiles\nSaving file "{ratio_file_name}"...\nIt will take some time...')
    for line in energy_file:
        line = line.rstrip('\n')
        energies = [float(energy) for energy in line.split('\t')]
        if len(energies) < 6:
            for energy in range(len(energies), 6):
                energies.append(0)
        ratios = [energies[0]]
        for low in range(1, 6):
            for high in range(low + 1, 6):
                if energies[low] == 0:
                    ratios.append(0)
                else:
                    ratios.append(energies[high] / energies[low])

        for level in range(len(ratios)):
            if level == len(ratios) - 1:
                ratio_file.write(common.value_to_write(ratios[level], '\n'))
            else:
                ratio_file.write(common.value_to_write(ratios[level], '\t'))

    ratio_file.close()
    energy_file.close()
    print(f'File "{ratio_file_name}" is saved.')


def find_cross(experimental_value, crystal, rare_earth, w):
    accuracy = 0.005
    points = []
    for _w in (abs(w), -abs(w)):
        ratio_file_name = common.get_paths(common.path_to_ratio_datafiles, 'ratio', 'dat', crystal, rare_earth, _w)
        ratio_file = open(ratio_file_name, 'r', encoding='utf-8')
        for line in ratio_file:
            line = line.rstrip('\n')
            numbers = [float(number) for number in line.split('\t')]
            x = numbers[0]
            ratios = numbers[1:]
            if any(abs(experimental_value - value) < accuracy for value in ratios):
                for index, ratio in enumerate(ratios):
                    if abs(experimental_value - ratio) < accuracy:
                        current = cross_point(rare_earth=rare_earth, w=_w, x=x, ratio_name=common.ratios_names[index],
                                              difference=experimental_value - ratio)
                        if not points:
                            points.append(current)
                        else:
                            previous = points[len(points) - 1]
                            if (current.rare_earth == previous.rare_earth and current.w == previous.w and
                                    current.ratio_name == previous.ratio_name and abs(current.x - previous.x) < 1e-3):
                                current_x = (current.x * previous.difference - previous.x * current.difference) / \
                                            (previous.difference - current.difference)
                                points[-1] = cross_point(rare_earth=rare_earth, w=_w, x=current_x, difference=0,
                                                         ratio_name=common.ratios_names[index])
                            else:
                                points.append(current)

    return points


def recalculation(points, experimental_energy, crystal, rare_earth):
    recalculated = []
    for point in points:
        x = point.x
        old_w = point.w
        level = int(point.ratio_name[-1]) - 1
        cef = get_object_with_parameters(crystal, rare_earth, old_w, x)
        old_energy = cef.get_energies()[level]
        new_w = experimental_energy / old_energy
        new_w = -new_w if old_w < 0 else new_w
        recalculated.append(cross_point(rare_earth=rare_earth, w=new_w, x=x, difference=point.difference,
                                        ratio_name=point.ratio_name))
    return recalculated


if __name__ == '__main__':
    crystal_name = 'YNi2'
    rare_earth_name = 'Tm'
    current_w = 1
    crosses = find_cross(4.59, crystal_name, rare_earth_name, current_w)
    recalculated_crosses = recalculation(crosses, 0.45606, crystal_name, rare_earth_name)
    print(*crosses, sep='\n')
    print(*recalculated_crosses, sep='\n')
