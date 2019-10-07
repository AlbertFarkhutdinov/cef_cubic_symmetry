import numpy
<<<<<<< HEAD
import os
=======
>>>>>>> origin/master
from datetime import *
from cef_object_scripts import core
from cef_object_scripts import work_with_cef_object as cef
from cef_object_scripts.get_cef_object import CF


<<<<<<< HEAD
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


=======
>>>>>>> origin/master
def get_one_dot():
    rare_earth = core.check_input('rare')
    w = core.check_input('w')
    x = core.check_input('x')
<<<<<<< HEAD
    cef_object = get_object_with_parameters(w, x, rare_earth)
    energies = cef_object.get_energies()
    print(f'x:\t\t{x: 9.3f}')
    for i in range(len(energies)):
        print(f'E[{i + 1}]:\t{energies[i]: 9.3f} meV')
=======
    cef.ratios_to_console(x, cef.get_ratios(w, x, rare_earth, 1), cef.get_ratios(-w, x, rare_earth, 1))
>>>>>>> origin/master


def load_data(crystal, rare_earth, w, x):
    file_name = core.get_paths(core.path_to_saved_objects, 'parameters', 'cfg', crystal, rare_earth, w, x)
    return CF(name=f'{crystal}: {rare_earth}3+', rare_earth=rare_earth, par_file=file_name)


def save_energy_dat(crystal, rare_earth, w, number_of_intervals):
    start_time = datetime.now()
    print('\nSaving energy datafiles')
    file_name = core.get_paths(core.path_to_energy_datafiles, 'energy', 'dat', crystal, rare_earth, w)
    core.check_path(file_name)
<<<<<<< HEAD
    if os.path.exists(file_name):
        os.remove(file_name)
=======
>>>>>>> origin/master
    my_file = open(file_name, 'a', encoding='utf-8')
    print(f'Saving file "{file_name}"...')
    for dot in range(number_of_intervals + 1):
        x = -1 + 2 * dot / number_of_intervals
<<<<<<< HEAD
        cef_object = get_object_with_parameters(w, x, rare_earth)
        energies = cef_object.get_energies()
=======
        energies = cef.get_energies(w, x, rare_earth)
>>>>>>> origin/master
        my_file.write(core.value_to_write(x, '\t'))
        for level in range(len(energies)):
            if level == len(energies) - 1:
                my_file.write(core.value_to_write(energies[level], '\n'))
            else:
                my_file.write(core.value_to_write(energies[level], '\t'))
    my_file.close()
    print(f'File "{file_name}" is saved')
    finish_time = datetime.now()
<<<<<<< HEAD
    print(f'Saving time: {finish_time - start_time}')


def save_parameters(crystal, rare_earth, w, x):
    start_time = datetime.now()
    print('\nSaving CEF parameters')
    cef_object = get_object_with_parameters(current_w, current_x, rare_earth_name)
=======
    print(f'Execution time: {finish_time - start_time}')


def save_parameters(crystal, rare_earth, w, x):
    print('\nSaving CEF parameters')
    cef_object = cef.get_object_with_parameters(current_w, current_x, rare_earth_name)
>>>>>>> origin/master
    file_name = core.get_paths(core.path_to_saved_objects, 'parameters', 'cfg', crystal, rare_earth, w, x)
    print(f'Saving file "{file_name}"...')
    cef_object.save_to_par_file(file_name)
    print(f'File "{file_name}" is saved')
<<<<<<< HEAD
    finish_time = datetime.now()
    print(f'Saving time: {finish_time - start_time}')


def save_spectra(crystal, rare_earth, w, x, energies, temperatures):
    start_time = datetime.now()
    print('\nSaving neutron inelastic scattering spectra')
    cef_object = get_object_with_parameters(current_w, current_x, rare_earth_name)
=======


def save_spectra(crystal, rare_earth, w, x, energies, temperatures):
    print('\nSaving neutron inelastic scattering spectra')
    cef_object = cef.get_object_with_parameters(current_w, current_x, rare_earth_name)
>>>>>>> origin/master
    for temperature in temperatures:
        spectrum = cef_object.spectrum(energy=energies, gamma=0.5, temperature=temperature)
        file_name = core.get_paths(core.path_to_spectra_datafiles, 'spectrum', 'dat',
                                   crystal, rare_earth, w, x, temperature)
        core.check_path(file_name)
<<<<<<< HEAD
        core.save_numpy(file_name, core.create_table(energies, spectrum))
    finish_time = datetime.now()
    print(f'Saving time: {finish_time - start_time}')


def save_susceptibility(crystal, rare_earth, w, x, temperatures):
    start_time = datetime.now()
    print('\nSaving magnetic susceptibilities')
    cef_object = get_object_with_parameters(current_w, current_x, rare_earth_name)
=======
        core.save_numpy(file_name, energies, spectrum)


def save_susceptibility(crystal, rare_earth, w, x, temperatures):
    print('\nSaving magnetic susceptibilities')
    cef_object = cef.get_object_with_parameters(current_w, current_x, rare_earth_name)
>>>>>>> origin/master
    file_name = core.get_paths(core.path_to_susceptibility_datafiles, 'susceptibility', 'dat',
                               crystal, rare_earth, w, x)
    core.check_path(file_name)
    for axis in ['z', 'x']:
<<<<<<< HEAD
        table = core.create_table(temperatures,
                                  cef_object.chi_s(temperatures)[f'chi_curie_{axis}'],
                                  cef_object.chi_s(temperatures)[f'chi_van_vleck_{axis}'],
                                  cef_object.chi_s(temperatures)[f'chi_{axis}'])
        core.save_numpy(file_name.replace('.dat', f'_chi_{axis}.dat'), table)
    table = core.create_table(temperatures,
                              cef_object.chi_s(temperatures)['chi'],
                              cef_object.chi_s(temperatures)['inverse_chi'])
    core.save_numpy(file_name.replace('.dat', '_chi_total.dat'), table)
    finish_time = datetime.now()
    print(f'Saving time: {finish_time - start_time}')
=======
        start_time = datetime.now()
        core.save_numpy(file_name.replace('.dat', f'_chi_{axis}.dat'),
                        temperatures,
                        cef_object.chi_s(temperatures)[f'chi_curie_{axis}'],
                        cef_object.chi_s(temperatures)[f'chi_van_vleck_{axis}'],
                        cef_object.chi_s(temperatures)[f'chi_{axis}'])
        finish_time = datetime.now()
        print(f'Execution time: {finish_time - start_time}')

    start_time = datetime.now()
    core.save_numpy(file_name.replace('.dat', '_chi_total.dat'),
                    temperatures,
                    cef_object.chi_s(temperatures)['chi'],
                    cef_object.chi_s(temperatures)['inverse_chi'])
    finish_time = datetime.now()
    print(f'Execution time: {finish_time - start_time}')
>>>>>>> origin/master


if __name__ == '__main__':
    crystal_name = 'YNi2'
<<<<<<< HEAD
    for rare_earth_name in ['Tb', 'Er', 'Ho', 'Tm']:
        save_energy_dat(crystal_name, rare_earth_name, 1, 5000)
        save_energy_dat(crystal_name, rare_earth_name, -1, 5000)

=======
    rare_earth_name = 'Tb'
    # save_energy_dat(crystal_name, rare_earth_name, 1, 5000)
>>>>>>> origin/master
    current_w = 0.0931
    current_x = 0.3048
    energy_space = numpy.linspace(-5, 30, 10001)
    temperature_array = [5, 25, 50]
    temperature_space = numpy.linspace(0.01, 100.0, 300)
    save_parameters(crystal_name, rare_earth_name, current_w, current_x)
    save_spectra(crystal_name, rare_earth_name, current_w, current_x, energy_space, temperature_array)
    save_susceptibility(crystal_name, rare_earth_name, current_w, current_x, temperature_space)
