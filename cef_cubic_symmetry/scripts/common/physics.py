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
