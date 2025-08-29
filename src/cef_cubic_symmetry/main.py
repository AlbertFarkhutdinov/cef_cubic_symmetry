"""The module contains the procedure of CEF parameters calculation."""

import pandas as pd

from cef_cubic_symmetry.auxiliary import utils as ut
from cef_cubic_symmetry.core.sample import Sample
from cef_cubic_symmetry.scripts import plot_objects as gg
from cef_cubic_symmetry.scripts.cubic_cef_object import Cubic
from cef_cubic_symmetry.scripts.experiment_object import Experiment


@ut.get_time_of_execution
def get_fixed_results(
    rare_earth: str,
    properties: dict,
    crystal: str,
    choice: int = 0,
    *,
    only_plots: bool = True,
) -> None:
    """
    Calculate and save the fixed results.

    This function saves the dependence of transition energies, their ratio
    on parameter x to file and its graphs for specified RE ions.

    """
    material = Sample(
        crystal=crystal,
        rare_earth=rare_earth,
    )
    if not only_plots:
        for w_parameter in (1, -1):
            cubic_object = Cubic(
                material,
                llw_parameters={'w': w_parameter},
            )
            cubic_object.save_peak_dat(
                number_of_intervals=5000,
                choice=choice,
            )
            cubic_object.get_ratios(choice=choice)
    y_max = (
        properties['max_energy']
        if choice == 0
        else properties['max_intensity']
    )
    y_major = (
        properties['energy_locator']
        if choice == 0
        else properties['intensity_locator']
    )
    gg.get_llw_plot(
        material=material,
        y_max=y_max,
        y_major=y_major,
        y_minor=y_major // 5,
        choice=choice,
    )


@ut.get_time_of_execution
def main(rare_earth: str, properties: dict) -> None:
    """Procedure of CEF parameters defining with spectra saving."""
    experiment = Experiment(
        material=Sample(
            crystal='YNi2',
            rare_earth=rare_earth,
        ),
        experimental_energies=properties['experimental_energies'],
        temperatures=properties['temperatures'],
    )
    experiment.get_llw_ratios_plot(**properties['ratios'])
    try:
        for _key, _value in properties['experiment'].items():
            experiment.get_spectrum_experiment(
                spectrometer=_key,
                **_value,
            )
        for _key, _value in properties['experiment_diff'].items():
            experiment.get_spectrum_differences(
                spectrometer=_key,
                **_value,
            )
    except (FileNotFoundError, IndexError):
        pass
    recalculated_crosses = experiment.get_cross_points()
    experiment.get_spectrum_theory(
        recalculated_crosses,
        **properties['theory'],
    )
    experiment.get_intensity_on_temperature(
        crosses=recalculated_crosses,
        **properties['intensities'],
    )


def get_scheme() -> None:
    """Print level scheme for specified parameters."""
    results = {
        'Pr': {
            'w': -0.105,
            'x': -0.46,
        },
        'Nd': {
            'w': 0.147,
            'x': -0.748,
        },
    }
    for _key, _value in results.items():
        print(
            Cubic(
                Sample(
                    rare_earth=_key,
                    crystal='YNi2',
                ),
                llw_parameters=_value,
            ),
        )


def run(is_recalculated: bool = False) -> None:
    if is_recalculated:
        fixed_properties = ut.get_json_object('json/fixed.json')
        for key, value in fixed_properties.items():
            get_fixed_results(
                rare_earth=key,
                properties=value,
                crystal='YNi2',
                only_plots=True,
                choice=0,
            )

    rare_earths_properties = ut.get_json_object('json/properties.json')

    for key, value in rare_earths_properties.items():
        main(rare_earth=key, properties=value)

    if is_recalculated:
        get_scheme()


if __name__ == '__main__':
    run(is_recalculated=True)
