import numpy
import os
from datetime import *
BASE_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
path_to_energy_datafiles = os.path.join(BASE_DIR, 'energy_datafiles')
path_to_ratio_datafiles = os.path.join(BASE_DIR, 'ratio_datafiles')
path_to_saved_objects = os.path.join(BASE_DIR, 'saved_objects')
path_to_spectra_datafiles = os.path.join(BASE_DIR, 'spectra')
path_to_susceptibility_datafiles = os.path.join(BASE_DIR, 'susceptibilities')
path_to_graphs = os.path.join(BASE_DIR, 'graphs')
ratios_names = ('E2/E1', 'E3/E1', 'E4/E1', 'E5/E1', 'E3/E2', 'E4/E2', 'E5/E2', 'E4/E3', 'E5/E3', 'E5/E4')


def get_sign(value):
    return '-' if value < 0 else '+'


def value_to_write(value, text_separator):
    return f'{value:10.5f}{text_separator}'


def get_paths(directory, data_name, format_name, crystal, rare_earth, w=None, x=None, temperature=None):
    os.chdir(BASE_DIR)
    short_name = f'{crystal}_{rare_earth}'
    full_name = short_name
    if w:
        full_name += f'_w{get_sign(w)}' + f'{abs(w): .3f}'.lstrip(' ')
    if x:
        full_name += f'_x{get_sign(x)}' + f'{abs(x): .3f}'.lstrip(' ')
    if temperature:
        full_name += f'_T{temperature}'
    return os.path.join(directory, short_name, f'{data_name}_{full_name}.{format_name}')


def check_path(path_to_check):
    condition = False
    levels = []
    while not condition:
        if os.path.exists(path_to_check):
            condition = True
        else:
            levels.insert(0, path_to_check)
            path_to_check = os.path.dirname(path_to_check)
    for level in range(len(levels) - 1):
        os.mkdir(levels[level])


def remove_if_exists(path_to_check):
    check_path(path_to_check)
    if os.path.exists(path_to_check):
        os.remove(path_to_check)


def check_input(choice):
    result = 0
    condition = False
    while not condition:
        if choice == 'rare':
            result = input('Input the name of RE ion (Pr, Nd, Pm, Sm, Gd, Tb, Dy, Ho, Er, Tm, Yb): ').capitalize()
            condition = (result in ['Pr', 'Nd', 'Pm', 'Sm', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb'])
        try:
            if choice == 'w':
                result = float(input('Input |W| > 0: '))
                condition = (result > 0)
            if choice == 'x':
                result = float(input('Input -1 <= x <= 1: '))
                condition = (-1 <= result <= 1)
            if choice == 'intervals':
                result = int(input('Input the number of intervals > 0: '))
                condition = (result > 0)
        except ValueError:
            condition = False
    return result


def thermodynamics(temperature, energies=None):
    thermal_dict = {'temperature': temperature / 11.6045}
    if energies.any() and thermal_dict['temperature'] > 0:
        thermal_dict['bolzmann'] = numpy.exp(-energies / thermal_dict['temperature'])
    return thermal_dict


def get_temperature(old_value, new_value):
    return new_value if (old_value is None) else old_value


def gauss(argument, center, sigma):
    return numpy.exp(-(argument - center) ** 2 / (2 * sigma ** 2)) / (sigma * numpy.sqrt(2 * numpy.pi))


def lorentz(argument, center, gamma):
    return (gamma / numpy.pi) / ((argument - center) ** 2 + gamma ** 2)


def pseudo_voigt(argument, center, sigma, gamma):
    width_gauss = 2 * sigma * numpy.sqrt(2 * numpy.log(2))
    width_lorentz = 2 * gamma
    full_width_at_half_maximum = (width_gauss ** 5 + 2.69269 * width_gauss ** 4 * width_lorentz +
                                  2.42843 * width_gauss ** 3 * width_lorentz ** 2 +
                                  4.47163 * width_gauss ** 2 + width_lorentz ** 3 +
                                  0.07842 * width_lorentz ** 4 * width_gauss + width_lorentz ** 5) ** 0.2
    ratio = width_lorentz / full_width_at_half_maximum
    eta = 1.36603 * ratio - 0.47719 * ratio ** 2 + 0.11116 * ratio ** 3
    return (1 - eta) * gauss(argument, center, sigma) + eta * lorentz(argument, center, gamma)


def create_table(*arrays):
    data = []
    for array in arrays:
        data.append(array)
    return numpy.transpose(data)


def user_input():
    crystal = input('Input the name of crystal (e.g. "YNi2"): ')
    rare_earth = check_input('rare')
    w = check_input('w')
    x = check_input('x')
    return crystal, rare_earth, w, x


def get_empty_matrix(size, dimension=2):
    if dimension == 2:
        return numpy.zeros((size, size), dtype='float64')
    if dimension == 1:
        return numpy.zeros(size, dtype='float64')


def get_time_of_execution(function):
    def wrapper(*args, **kwargs):
        start_time = datetime.now()
        function(*args, **kwargs)
        finish_time = datetime.now()
        print(f'Saving time: {finish_time - start_time}')
    return wrapper


def saving_file(function):
    def wrapper(file, object_to_save, *args, **kwargs):
        print(f'\n***Saving {object_to_save}***:\nSaving file "{file}"...\nIt may take some time...')
        function(file, object_to_save, *args, **kwargs)
        print(f'File "{file}" is saved.')
    return wrapper


if __name__ == '__main__':
    print(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
    print(os.getcwd())
    print(path_to_energy_datafiles)
