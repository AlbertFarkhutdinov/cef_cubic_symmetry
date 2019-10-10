import numpy
import os
from datetime import *
from cef_object_scripts import core
from cef_object_scripts.get_cef_object import CF


def get_object_with_parameters(w, x, rare_earth):
    f4 = 60
    f6 = {'Pr': 1260, 'Nd': 2520, 'Pm': 1260, 'Gd': 1260, 'Tb': 7560, 'Dy': 13860, 'Ho': 13860, 'Er': 13860, 'Tm': 7560,
          'Yb': 1260}
    cef_object = CF(name=f'YNi2: {rare_earth}3+', rare_earth=rare_earth)
    cef_object.B40 = w * x / f4
    cef_object.B44 = 5 * cef_object.B40
    cef_object.B60 = w * (1 - abs(x)) / f6[rare_earth]
    cef_object.B64 = -21 * cef_object.B60
    return cef_object


def get_one_dot(w, x, rare_earth):
    cef_object = get_object_with_parameters(w, x, rare_earth)
    energies = cef_object.get_energies()
    print(f'x:\t\t{x: 9.3f}')
    for i in range(len(energies)):
        print(f'E[{i + 1}]:\t{energies[i]: 9.3f} meV')


def load_data(crystal, rare_earth, w, x):
    file_name = core.get_paths(core.path_to_saved_objects, 'parameters', 'cfg', crystal, rare_earth, w, x)
    return CF(name=f'{crystal}: {rare_earth}3+', rare_earth=rare_earth, par_file=file_name)


def save_energy_dat(crystal, rare_earth, w, number_of_intervals):
    x_space = numpy.linspace(-1, 1, number_of_intervals + 1)
    start_time = datetime.now()
    print('\nSaving energy datafiles')
    file_name = core.get_paths(core.path_to_energy_datafiles, 'energy', 'dat', crystal, rare_earth, w)
    core.check_path(file_name)
    if os.path.exists(file_name):
        os.remove(file_name)
    my_file = open(file_name, 'a', encoding='utf-8')
    print(f'Saving file "{file_name}"...')
    print('It will take some time...')
    for x in x_space:
        cef_object = get_object_with_parameters(w, x, rare_earth)
        energies = cef_object.get_energies()
        my_file.write(core.value_to_write(x, '\t'))
        for level in range(len(energies)):
            if level == len(energies) - 1:
                my_file.write(core.value_to_write(energies[level], '\n'))
            else:
                my_file.write(core.value_to_write(energies[level], '\t'))
    my_file.close()
    print(f'File "{file_name}" is saved')
    finish_time = datetime.now()
    print(f'Saving time: {finish_time - start_time}')


def save_parameters(crystal, rare_earth, w, x):
    start_time = datetime.now()
    print('\nSaving CEF parameters')
    cef_object = get_object_with_parameters(w, x, rare_earth)
    file_name = core.get_paths(core.path_to_saved_objects, 'parameters', 'cfg', crystal, rare_earth, w, x)
    print(f'Saving file "{file_name}"...')
    cef_object.save_to_par_file(file_name)
    print(f'File "{file_name}" is saved')
    finish_time = datetime.now()
    print(f'Saving time: {finish_time - start_time}')


def save_spectra(crystal, rare_earth, w, x):
    energies = numpy.linspace(-5, 30, 10001)
    temperatures = [5, 25, 50]
    start_time = datetime.now()
    print('\nSaving neutron inelastic scattering spectra')
    cef_object = get_object_with_parameters(w, x, rare_earth)
    for temperature in temperatures:
        spectrum = cef_object.spectrum(energy=energies, gamma=0.5, temperature=temperature)
        file_name = core.get_paths(core.path_to_spectra_datafiles, 'spectrum', 'dat',
                                   crystal, rare_earth, w, x, temperature)
        core.check_path(file_name)
        core.save_numpy(file_name, core.create_table(energies, spectrum))
        print(f'File "{file_name}" is saved')
    finish_time = datetime.now()
    print(f'Saving time: {finish_time - start_time}')


def save_susceptibility(crystal, rare_earth, w, x):
    temperatures = numpy.linspace(0.1, 100.0, 300)
    start_time = datetime.now()
    print('\nSaving magnetic susceptibilities')
    common_file_name = core.get_paths(core.path_to_susceptibility_datafiles, 'susceptibility', 'dat',
                                      crystal, rare_earth, w, x)
    cef_object = get_object_with_parameters(w, x, rare_earth)
    chi_s = cef_object.chi_s(temperatures)
    for axis in ['z', 'x', 'total']:
        file_name = common_file_name.replace('.dat', f'_chi_{axis}.dat')
        core.check_path(file_name)
        if os.path.exists(file_name):
            os.remove(file_name)
        my_file = open(file_name, 'a', encoding='utf-8')
        print(f'Saving file "{file_name}"...')
        if axis == 'z' or axis == 'x':
            my_file.write(f'T(Kelvin)\tchi_curie_{axis}\tchi_van_vleck_{axis}\tchi_{axis}\n')
        else:
            my_file.write('T(Kelvin)\tchi_total\tinverse_chi\n')
        for i in range(len(temperatures)):
            my_file.write(core.value_to_write(temperatures[i], '\t'))
            if axis == 'z' or axis == 'x':
                my_file.write(core.value_to_write(chi_s[f'chi_curie_{axis}'][i], '\t'))
                my_file.write(core.value_to_write(chi_s[f'chi_van_vleck_{axis}'][i], '\t'))
                my_file.write(core.value_to_write(chi_s[f'chi_{axis}'][i], '\n'))
            else:
                my_file.write(core.value_to_write(chi_s['chi'][i], '\t'))
                my_file.write(core.value_to_write(chi_s['inverse_chi'][i], '\n'))
        my_file.close()
        print(f'File "{file_name}" is saved')
    finish_time = datetime.now()

    print(f'Saving time: {finish_time - start_time}')


if __name__ == '__main__':
    crystal_name = 'YNi2'
    for rare_earth_name in ['Tb', 'Er', 'Ho', 'Tm']:
        save_energy_dat(crystal_name, rare_earth_name, 1, 5000)
        save_energy_dat(crystal_name, rare_earth_name, -1, 5000)
    current_w = 0.0931
    current_x = 0.3048
    save_parameters(crystal_name, rare_earth_name, current_w, current_x)
    save_spectra(crystal_name, rare_earth_name, current_w, current_x)
    save_susceptibility(crystal_name, rare_earth_name, current_w, current_x)
