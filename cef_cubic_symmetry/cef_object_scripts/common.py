"""The module contains some common constants and functions that used in this project."""
import os
from datetime import datetime
import numpy
from cef_object_scripts.tabular_information import ACCEPTABLE_RARE_EARTHS

BASE_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
PATH_TO_ENERGY_DATAFILES = os.path.join(BASE_DIR, 'energy_datafiles')
PATH_TO_RATIO_DATAFILES = os.path.join(BASE_DIR, 'ratio_datafiles')
PATH_TO_SAVED_OBJECTS = os.path.join(BASE_DIR, 'saved_objects')
PATH_TO_SPECTRA_DATAFILES = os.path.join(BASE_DIR, 'spectra')
PATH_TO_SUSCEPTIBILITY_DATAFILES = os.path.join(BASE_DIR, 'susceptibilities')
PATH_TO_GRAPHS = os.path.join(BASE_DIR, 'graphs')
RATIOS_NAMES = ('E2/E1', 'E3/E1', 'E4/E1', 'E5/E1',
                'E3/E2', 'E4/E2', 'E5/E2',
                'E4/E3', 'E5/E3',
                'E5/E4')


def get_sign(value):
    """Returns minus, if argument is negative,
    else it returns plus."""
    return '-' if value < 0 else '+'


def get_value(value):
    """Returns float number as a string with sign plus or minus."""
    return f'{get_sign(value)}' + f'{abs(value): .3f}'.lstrip(' ')


def value_to_write(value, text_separator):
    """Returns float number as a string
    with tabulation symbol or newline one in the end."""
    return f'{value:10.5f}{text_separator}'


def get_paths(directory, data_name, format_name, material: dict = None, parameters: dict = None):
    """Returns path of the file that will be saved."""
    os.chdir(BASE_DIR)
    short_name = ''
    if material:
        short_name = f'{material["crystal"]}_{material["rare_earth"]}'
    full_name = short_name
    if parameters:
        for key, value in parameters.items():
            if key != 'T':
                full_name += f'_{key}{get_value(value)}'
            else:
                full_name += f'_{key}{value}'
    return os.path.join(directory, short_name, f'{data_name}_{full_name}.{format_name}')


def check_path(path_to_check):
    """Makes parent directories for argument, if they do not exist."""
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
    """Makes remove file, if it exists."""
    check_path(path_to_check)
    if os.path.exists(path_to_check):
        os.remove(path_to_check)


def check_input(choice):
    """Checks a value inputted by user, returns it,
    if it satisfies the condition, else requests input again."""
    result = 0
    condition = False
    while not condition:
        if choice == 'rare':
            request = f'Input the name of RE ion ({", ". join(ACCEPTABLE_RARE_EARTHS)}): '
            result = input(request).capitalize()
            condition = (result in ACCEPTABLE_RARE_EARTHS)
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
    """Returns dictionary including value of temperature in meV and Bolzmann factor."""
    thermal_dict = {'temperature': temperature / 11.6045}
    if energies is not None and thermal_dict['temperature'] > 0:
        zeros = numpy.zeros(len(energies))
        thermal_dict['bolzmann'] = numpy.exp(zeros - energies / thermal_dict['temperature'])
    return thermal_dict


def get_temperature(old_value, new_value):
    """Returns new_value if old_value is None, else it returns old_value."""
    return new_value if (old_value is None) else old_value


def gauss(argument, center, sigma):
    """Returns value of Gauss function."""
    return (numpy.exp(-(argument - center) ** 2 / (2 * sigma ** 2)) /
            (sigma * numpy.sqrt(2 * numpy.pi)))


def lorentz(argument, center, gamma):
    """Returns value of Lorentz function."""
    return (gamma / numpy.pi) / ((argument - center) ** 2 + gamma ** 2)


def pseudo_voigt(argument, center, sigma, gamma):
    """Returns value of pseudo-Voigt function."""
    width_gauss = 2 * sigma * numpy.sqrt(2 * numpy.log(2))
    width_lorentz = 2 * gamma
    full_width_at_half_maximum = (width_gauss ** 5 +
                                  2.69269 * width_gauss ** 4 * width_lorentz +
                                  2.42843 * width_gauss ** 3 * width_lorentz ** 2 +
                                  4.47163 * width_gauss ** 2 +
                                  width_lorentz ** 3 +
                                  0.07842 * width_lorentz ** 4 * width_gauss +
                                  width_lorentz ** 5) ** 0.2
    ratio = width_lorentz / full_width_at_half_maximum
    eta = 1.36603 * ratio - 0.47719 * ratio ** 2 + 0.11116 * ratio ** 3
    return (1 - eta) * gauss(argument, center, sigma) + eta * lorentz(argument, center, gamma)


def create_table(*arrays):
    """Return transposed array of arrays."""
    data = []
    for array in arrays:
        data.append(array)
    return numpy.transpose(data)


def user_input():
    """Returns data inputted by user."""
    crystal = input('Input the name of crystal (e.g. "YNi2"): ')
    rare_earth = check_input('rare')
    w_parameter = check_input('w')
    x_parameter = check_input('x')
    return crystal, rare_earth, w_parameter, x_parameter


def get_empty_matrix(size, dimension=2):
    """Returns 1D or 2D array filled by zeros."""
    empty_matrix = None
    if dimension == 2:
        empty_matrix = numpy.zeros((size, size), dtype='float64')
    if dimension == 1:
        empty_matrix = numpy.zeros(size, dtype='float64')
    return empty_matrix


def get_time_of_execution(function):
    """Prints time of function's execution."""
    def wrapper(*args, **kwargs):
        start_time = datetime.now()
        function(*args, **kwargs)
        finish_time = datetime.now()
        print(f'Saving time: {finish_time - start_time}')
    return wrapper


def saving_file(function):
    """Prints information about saved file."""
    def wrapper(file, object_to_save, *args, **kwargs):
        print('\n',
              f'***Saving {object_to_save}***:',
              f'Saving file "{file}"...',
              'It may take some time...', sep='\n')
        function(file, object_to_save, *args, **kwargs)
        print(f'File "{file}" is saved.')
    return wrapper


if __name__ == '__main__':
    print(get_paths('directory', 'data_name', 'format_name',
                    material={'crystal': 'YNi2', 'rare_earth': 'Tb'},
                    parameters={'w': 0.3243, 'x': -0.4687, 'T': 5}))
