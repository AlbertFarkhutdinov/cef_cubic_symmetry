import numpy
import os
BASE_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
path_to_energy_datafiles = os.path.join(BASE_DIR, 'energy_datafiles')
path_to_saved_objects = os.path.join(BASE_DIR, 'saved_objects')
path_to_spectra_datafiles = os.path.join(BASE_DIR, 'spectra')


def get_sign(value):
    return '-' if value < 0 else '+'


def value_to_write(value, text_separator):
    return f'{value:9.5f}{text_separator}'


def get_paths(crystal, rare_earth, w=None, x=None, temperature=None):
    config_file = 'parameters_'
    energy_file = 'energy_'
    spectrum_file = 'spectrum_'
    paths = {'short_name': f'{crystal}_{rare_earth}'}
    paths['full_name'] = paths['short_name'] 
    if w:
        paths['full_name'] += f'_w{get_sign(w)}{abs(w)}'
        energy_file += f'_w{get_sign(w)}{abs(w)}'
    if x:
        paths['full_name'] += f'_x{get_sign(x)}{abs(x)}'
    if temperature:
        paths['full_name'] += f'_T{temperature}'
    paths['config_file_path'] = os.path.join(path_to_saved_objects, paths['short_name'], 
                                             config_file + paths['full_name'] + '.cfg')
    paths['spectrum_file_path'] = os.path.join(path_to_spectra_datafiles, paths['short_name'], 
                                               spectrum_file + paths['full_name'] + '.dat')
    paths['energy_file_path'] = os.path.join(path_to_energy_datafiles, 
                                             energy_file + paths['full_name'] + '.dat')
    return paths


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


def kelvin_to_mev(temperature):
    return temperature / 11.6045


def bolzman(energies, temperature):
    return numpy.exp(-energies / kelvin_to_mev(temperature))


def get_temperature(checked, new):
    if checked is None:
        checked = new
    return checked


def gauss(argument, center, sigma):
    return numpy.exp(-(argument - center) ** 2 / (2 * sigma ** 2)) / (sigma * numpy.sqrt(2 * numpy.pi))


def lorentz(argument, center, gamma):
    return (gamma / numpy.pi) / ((argument - center) ** 2 + gamma ** 2)


def pseudo_voigt(argument, center, sigma, gamma):
    gamma_g = 2 * numpy.sqrt(2 * numpy.log(2)) * sigma
    gamma_l = 2 * gamma
    full_width_at_half_maximum = (gamma_g ** 5 + 2.69269 * gamma_g ** 4 * gamma_l +
                                  2.42843 * gamma_g ** 3 * gamma_l ** 2 +
                                  4.47163 * gamma_g ** 2 + gamma_l ** 3 +
                                  0.07842 * gamma_g ** 4 * gamma_l + gamma_l ** 5) ** 0.2
    ratio = gamma_l / full_width_at_half_maximum
    fraction = 1.36603 * ratio - 0.47719 * ratio ** 2 + 0.11116 * ratio ** 3
    return (1 - fraction) * gauss(argument, center, sigma) + fraction * lorentz(argument, center, gamma)


if __name__ == '__main__':
    print('Done!')
