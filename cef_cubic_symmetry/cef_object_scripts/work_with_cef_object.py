from cef_object_scripts import core
from cef_object_scripts.get_cef_object import CF


def get_object_with_parameters(w, x, rare_earth):
    f4 = 60
    f6 = {'Pr': 1260, 'Nd': 2520, 'Pm': 1260, 'Gd': 1260, 'Tb': 7560, 'Dy': 13860, 'Ho': 13860, 'Er': 13860, 'Tm': 7560,
          'Yb': 1260}  # add check for Ce, Sm, Eu
    cef_object = CF(name=f'YNi2: {rare_earth}3+', rare_earth=rare_earth)
    cef_object.B40 = w * x / f4
    cef_object.B44 = 5 * cef_object.B40
    cef_object.B60 = w * (1 - abs(x)) / f6[rare_earth]
    cef_object.B64 = -21 * cef_object.B60
    return cef_object


def get_energies(w, x, rare_earth):
    cef_object = get_object_with_parameters(w, x, rare_earth)
    return [peak[0] for peak in cef_object.get_peaks()]


def get_intensities(w, x, rare_earth):
    cef_object = get_object_with_parameters(w, x, rare_earth)
    return [peak[1] for peak in cef_object.get_peaks()]


def get_ratios(w, x, rare_earth, number_of_ratio):
    energies = get_energies(w, x, rare_earth)
    level = 0
    while energies[level] < 2:  # ? minimum energy due to elastic scattering peak
        level += 1

    if len(energies) == level + 1:
        return 0
    else:
        ratios = [energies[level + i] / energies[level] for i in range(len(energies) - level)]
        return ratios[number_of_ratio]


def energy_to_console(value, text_separator):
    print(core.value_to_write(value, text_separator), end='')


def ratios_to_console(value1, value2, value3):
    print(core.value_to_write(value1, '\t'), core.value_to_write(value2, '\t'), core.value_to_write(value3, '\n'))


if __name__ == '__main__':
    print('Done!')
