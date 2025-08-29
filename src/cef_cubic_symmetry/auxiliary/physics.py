"""The module contains physics functions that used in this project."""

from numpy import exp, log, pi, sqrt, zeros


def gaussian_normalized(
    argument: float,
    center: float,
    sigma: float,
) -> float:
    """
    Return value of normalized Gauss function.

    Integral with infinite limits equals 1.
    Full width at half-maximum (FWHM) is 2 * sqrt(2*ln(2)) * sigma.
    Value at the maximum is 1 / (sigma * sqrt(2 * pi)).

    """
    return (
        exp(-(argument - center) ** 2 / (2 * sigma ** 2))
        / (sigma * sqrt(2 * pi))
    )


def lorentzian_normalized(
    argument: float,
    center: float,
    gamma: float,
) -> float:
    """
    Return value of normalized Lorentz function.

    Integral with infinite limits equals 1.
    Full width at half-maximum (FWHM) is 2 * gamma.
    Value at the maximum is 1 / (pi * gamma).

    """
    return (gamma / pi) / ((argument - center) ** 2 + gamma ** 2)


def gaussian(
    arg: float,
    center: float,
    width: float,
    amplitude: float,
) -> float:
    """
    Return value of Gauss function.

    Integral with infinite limits equals amplitude.
    Value at the maximum is amplitude / (sigma * sqrt(2 * pi)).

    """
    return amplitude * gaussian_normalized(arg, center, width)


def lorentzian(
    arg: float,
    center: float,
    width: float,
    amplitude: float,
) -> float:
    """
    Return value of Lorentz function.

    Integral with infinite limits equals amplitude.
    Value at the maximum is amplitude / (pi * gamma).

    """
    return amplitude * lorentzian_normalized(arg, center, width)


def pseudo_voigt_normalized(
    argument: float,
    center: float,
    sigma: float,
    gamma: float,
) -> float:
    """
    Return value of normalized pseudo-Voigt function.

    Integral with infinite limits equals 1.
    Full width at half-maximum (FWHM) is fwhm_total.

    """
    fwhm_gauss = 2 * sigma * sqrt(2 * log(2))
    fwhm_lorentz = 2 * gamma
    fwhm_total = (
        (
            fwhm_gauss ** 5
            + 2.69269 * fwhm_gauss ** 4 * fwhm_lorentz
            + 2.42843 * fwhm_gauss ** 3 * fwhm_lorentz ** 2
            + 4.47163 * fwhm_gauss ** 2 * fwhm_lorentz ** 3
            + 0.07842 * fwhm_gauss * fwhm_lorentz ** 4
            + fwhm_lorentz ** 5
        ) ** 0.2
    )
    ratio = fwhm_lorentz / fwhm_total
    eta = 1.36603 * ratio - 0.47719 * ratio ** 2 + 0.11116 * ratio ** 3
    return (
        (1 - eta) * gaussian_normalized(argument, center, sigma)
        + eta * lorentzian_normalized(argument, center, gamma)
    )


def multi_peak(
    function,
    arg: float,
    *parameters,
) -> float:
    """Return value of several peaks sum."""
    background = parameters[0]
    peaks_parameters = parameters[1:]
    peaks = [
        function(arg, *peaks_parameters[index: index + 3])
        for index in range(0, len(peaks_parameters), 3)
    ]
    return background + sum(peaks)


def multi_lorentzian(
    arg: float,
    *parameters,
) -> float:
    """Return value of multi_peak function for lorentzian."""
    return multi_peak(lorentzian, arg, *parameters)


def multi_gaussian(
    arg: float,
    *parameters,
) -> float:
    """Return value of multi_peak function for gaussian."""
    return multi_peak(gaussian, arg, *parameters)


def thermodynamics(
    temperature: float,
    energies=None,
) -> dict[str, float]:
    """Return temperature value in meV and Boltzmann factor."""
    thermal_dict = {'temperature': temperature / 11.6045}
    if energies is not None and thermal_dict['temperature'] > 0:
        zero_array = zeros(len(energies))
        thermal_dict['boltzmann'] = exp(
            zero_array - energies / thermal_dict['temperature'],
        )
    return thermal_dict


def lowering_operator(
    initial_number: float,
    squared_j: float,
    degree: int,
) -> float:
    """Return the result of the lowering operator's."""
    result = 1
    for step in range(degree):
        result *= (
            squared_j
            - (initial_number - step)
            * (initial_number - step - 1)
        )
    return sqrt(result)


def steven_operators(
    key: str,
    squared_j: float,
    mqn1,
    mqn2=None,
) -> float:
    """Return the result of the Stevens operators'."""
    result = {
        'o20': lambda: 3 * mqn1[2] - squared_j,
        'o40': lambda: (
            35 * mqn1[4]
            - 30 * squared_j * mqn1[2]
            + 25 * mqn1[2]
            - 6 * squared_j
            + 3 * squared_j ** 2
        ),
        'o60': lambda: (
            231 * mqn1[1] ** 6
            - 315 * squared_j * mqn1[4]
            + 735 * mqn1[4]
            + 105 * squared_j ** 2 * mqn1[2]
            - 525 * squared_j * mqn1[2]
            + 294 * mqn1[2]
            - 5 * squared_j ** 3
            + 40 * squared_j ** 2
            - 60 * squared_j
        ),
    }
    if mqn2:
        result['o22'] = lambda: 0.5 * lowering_operator(mqn2[1], squared_j, 2)
        result['o43'] = lambda: (
            0.25 * lowering_operator(mqn2[1], squared_j, 3)
            * (mqn1[1] + mqn2[1])
        )
        result['o63'] = lambda: (
            0.25
            * (
                11 * (mqn1[1] ** 3 + mqn2[3])
                - 3 * (mqn1[1] + mqn2[1]) * squared_j
                - 59 * (mqn1[1] + mqn2[1])
            )
            * lowering_operator(mqn2[1], squared_j, 3)
        )
        result['o44'] = lambda: 0.5 * lowering_operator(mqn2[1], squared_j, 4)
        result['o66'] = lambda: 0.5 * lowering_operator(mqn2[1], squared_j, 6)
        result['o42'] = lambda: (
            (3.5 * (mqn1[2] + mqn2[2]) - squared_j - 5)
            * 0.5 * lowering_operator(mqn2[1], squared_j, 2)
        )
        result['o62'] = lambda: (
            (
                16.5 * (mqn1[4] + mqn2[4])
                - 9 * (mqn1[2] + mqn2[2]) * squared_j
                - 61.5 * (mqn1[2] + mqn2[2])
                + squared_j ** 2
                + 10 * squared_j + 102
            ) * (0.5 * lowering_operator(mqn2[1], squared_j, 2))
        )
        result['o64'] = lambda: (
            (5.5 * (mqn1[2] + mqn2[2]) - squared_j - 38)
            * 0.5 * lowering_operator(mqn2[1], squared_j, 4)
        )
    return result.get(key, lambda: None)()
