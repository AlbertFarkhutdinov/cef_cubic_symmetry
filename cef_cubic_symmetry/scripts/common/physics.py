"""The module contains physics functions that used in this project."""
from numpy import zeros, exp, sqrt, pi, log


def thermodynamics(temperature, energies=None):
    """Returns dictionary including value of temperature in meV and Bolzmann factor."""
    thermal_dict = {'temperature': temperature / 11.6045}
    if energies is not None and thermal_dict['temperature'] > 0:
        zero_array = zeros(len(energies))
        thermal_dict['bolzmann'] = exp(zero_array - energies / thermal_dict['temperature'])
    return thermal_dict


def gauss(argument, center, sigma):
    """Returns value of Gauss function."""
    return (exp(-(argument - center) ** 2 / (2 * sigma ** 2)) /
            (sigma * sqrt(2 * pi)))


def lorentz(argument, center, gamma):
    """Returns value of Lorentz function."""
    return (gamma / pi) / ((argument - center) ** 2 + gamma ** 2)


def pseudo_voigt(argument, center, sigma, gamma):
    """Returns value of pseudo-Voigt function."""
    width_gauss = 2 * sigma * sqrt(2 * log(2))
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


def lowering_operator(initial_number, squared_j, degree):
    """Сalculates the result of the lowering operator’s action
    on the wave function with quantum number initial_number"""
    result = 1
    for step in range(degree):
        result *= (squared_j -
                   (initial_number - step) *
                   (initial_number - step - 1))
    return sqrt(result)


def steven_operators(key, squared_j, mqn_1, mqn_2=None):
    """"""
    result = {
        'o20': lambda: 3 * mqn_1[2] - squared_j,
        'o40': lambda: (35 * mqn_1[4] -
                        30 * squared_j * mqn_1[2] +
                        25 * mqn_1[2] -
                        6 * squared_j +
                        3 * squared_j ** 2),
        'o60': lambda: (231 * mqn_1[1] ** 6 -
                        315 * squared_j * mqn_1[4] +
                        735 * mqn_1[4] +
                        105 * squared_j ** 2 * mqn_1[2] -
                        525 * squared_j * mqn_1[2] +
                        294 * mqn_1[2] -
                        5 * squared_j ** 3 +
                        40 * squared_j ** 2 -
                        60 * squared_j),
    }
    if mqn_2:
        result['o22'] = lambda: 0.5 * lowering_operator(mqn_2[1], squared_j, 2)
        result['o43'] = lambda: (0.25 * lowering_operator(mqn_2[1], squared_j, 3) *
                                 (mqn_1[1] + mqn_2[1]))
        result['o63'] = lambda: (0.25 * (11 * (mqn_1[1] ** 3 + mqn_2[3]) -
                                         3 * (mqn_1[1] + mqn_2[1]) * squared_j -
                                         59 * (mqn_1[1] + mqn_2[1])) *
                                 lowering_operator(mqn_2[1], squared_j, 3))
        result['o44'] = lambda: 0.5 * lowering_operator(mqn_2[1], squared_j, 4)
        result['o66'] = lambda: 0.5 * lowering_operator(mqn_2[1], squared_j, 6)
        result['o42'] = lambda: ((3.5 * (mqn_1[2] + mqn_2[2]) - squared_j - 5) *
                                 0.5 * lowering_operator(mqn_2[1], squared_j, 2))
        result['o62'] = lambda: (16.5 * (mqn_1[4] + mqn_2[4]) -
                                 9 * (mqn_1[2] + mqn_2[2]) * squared_j -
                                 61.5 * (mqn_1[2] + mqn_2[2]) +
                                 squared_j ** 2 +
                                 10 * squared_j +
                                 102) * (0.5 * lowering_operator(mqn_2[1], squared_j, 2))
        result['o64'] = lambda: ((5.5 * (mqn_1[2] + mqn_2[2]) - squared_j - 38) *
                                 0.5 * lowering_operator(mqn_2[1], squared_j, 4))
    return result.get(key, lambda: None)()


if __name__ == '__main__':
    for value in ('20', '40', '60', '22', '42', '62', '43', '63', '44', '64', '66'):
        print(steven_operators(f'o{value}', 56, [1, 1, 1, 1, 1], [1, 1, 1, 1, 1]))
