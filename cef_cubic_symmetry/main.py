"""
The module contains the procedure of CEF parameters defining
with spectra saving.

"""


from scripts.common.constants import Material
from scripts.common.utils import get_time_of_execution, get_json_object
from scripts.cubic_cef_object import Cubic
from scripts.experiment_object import Experiment
from scripts import plot_objects as gg


@get_time_of_execution
def get_fixed_results(
        rare_earth: str,
        properties: dict,
        crystal: str,
        only_plots=True,
        choice=0,
):
    """Saves the dependence of transition energies, their ratio
    on parameter x to file and its graphs for specified RE ions"""
    material = Material(
        crystal=crystal,
        rare_earth=rare_earth,
    )
    if not only_plots:
        for w_parameter in (1, -1):
            cubic_object = Cubic(
                material,
                llw_parameters={'w': w_parameter}
            )
            cubic_object.save_peak_dat(
                number_of_intervals=5000,
                choice=choice
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


@get_time_of_execution
def main(rare_earth: str, properties: dict):
    """Procedure of CEF parameters defining with spectra saving"""
    experiment = Experiment(
        material=Material(
            crystal='YNi2',
            rare_earth=rare_earth
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
        **properties['theory']
    )
    # experiment.get_intensity_on_temperature(
    #     crosses=recalculated_crosses,
    #     **properties['intensities']
    # )


def get_scheme():
    """Prints level scheme for specified parameters"""
    results = {
        'Pr': {
            'w': -0.105,
            'x': -0.460,
        },
        'Nd': {
            'w': 0.147,
            'x': -0.748,
        }
    }
    for _key, _value in results.items():
        print(
            Cubic(
                Material(
                    rare_earth=_key,
                    crystal='YNi2'
                ),
                llw_parameters=_value
            )
        )


if __name__ == '__main__':
    # FIXED_PROPS = get_json_object('fixed.json')
    # for key, value in FIXED_PROPS.items():
    #     get_fixed_results(
    #         rare_earth=key,
    #         properties=value,
    #         crystal='YNi2',
    #         only_plots=True,
    #         choice=0,
    #     )

    PROPS = get_json_object('properties.json')

    for key, value in PROPS.items():
        main(rare_earth=key, properties=value)

    # get_scheme()
