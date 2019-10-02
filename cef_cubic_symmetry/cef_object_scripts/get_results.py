import numpy
from cef_object_scripts import core
from cef_object_scripts.get_cef_object import CF
from cef_object_scripts import work_with_cef_object as cef 


def get_energy_dat(crystal, rare_earth, w, number_of_intervals):
    for dot in range(number_of_intervals + 1):
        x = -1 + 2 * dot / number_of_intervals
        energies = cef.get_energies(w, x, rare_earth)
        cef.energy_to_console(x, '\t')
        cef.energy_to_file(crystal, rare_earth, w, x, '\t')
        for level in range(len(energies)):
            if level == len(energies) - 1:
                cef.energy_to_console(energies[level], '\n')
                cef.energy_to_file(crystal, rare_earth, w, energies[level], '\n')
            else:
                cef.energy_to_console(energies[level], '\t')
                cef.energy_to_file(crystal, rare_earth, w, energies[level], '\t')

    print('Done!')


def get_one_dot():
    rare_earth = core.check_input('rare')
    w = core.check_input('w')
    x = core.check_input('x')
    cef.ratios_to_console(x, cef.get_ratios(w, x, rare_earth, 1), cef.get_ratios(-w, x, rare_earth, 1))


def load_data(crystal, rare_earth, w, x):
    par_file_path = core.get_paths(crystal, rare_earth, w, x)['config_file_path']
    return CF(name=f'{crystal}: {rare_earth}3+', rare_earth=rare_earth, par_file=par_file_path)


def save_data(crystal, rare_earth, w, x, temperatures_list):
    cef_object = cef.get_object_with_parameters(w, x, rare_earth)
    cef_object.save_to_par_file(core.get_paths(crystal, rare_earth, w, x)['config_file_path'])
    energies = numpy.linspace(-5, 30, 10001)
    for temperature in temperatures_list:
        spectrum_file_path = core.get_paths(crystal, rare_earth, w, x, temperature)['spectrum_file_path']
        core.check_path(spectrum_file_path)
        x_array = cef_object.nx_spectrum(energies, gamma=0.5, temperature=temperature).data.energy_transfer.nxdata
        y_array = cef_object.nx_spectrum(energies, gamma=0.5, temperature=temperature).data.intensity.nxdata
        data = numpy.transpose([x_array, y_array])
        numpy.savetxt(spectrum_file_path, data, delimiter='\t')
    # root = NXroot()
    # for temperature in temperatures_list:
    #     root.NXspectrum[temperature] = cef_object.nx_spectrum(energies, gamma=0.5, temperature=temperature)
    #
    # temperatures = numpy.linspace(0.01, 100.0, 300)
    # root.chi = cef_object.nx_chi(temperatures)
    # nxtree.cef_object = root


if __name__ == '__main__':
    crystal_name = 'YNi2'
    rare_earth_name = 'Tb'
    current_w = 0.0931
    current_x = 0.3048
    temperatures_array = [5, 25, 50]
    save_data(crystal_name, rare_earth_name, current_w, current_x, temperatures_array)
