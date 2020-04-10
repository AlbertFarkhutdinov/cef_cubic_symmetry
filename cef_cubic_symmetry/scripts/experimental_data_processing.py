"""The module contains the procedure of CEF parameters defining with spectra saving."""
import os
from copy import deepcopy
from scripts.common.constants import DATA_PATHS, Material, Data, Scale
from scripts.common.utils import get_time_of_execution
from scripts.common.fitting_utils import get_data_from_file
from scripts.cubic_cef_object import Cubic
from scripts import plot_objects as gg


class Experiment:
    """Class contains experimental parameters"""
    def __init__(self, material: Material, experimental_energies: tuple, temperatures: tuple):
        """Initialization of class Experiment"""
        self.material = material
        self.experimental_energies = experimental_energies
        self.experimental_ratio = experimental_energies[1] / experimental_energies[0]
        self.temperatures = temperatures
        self.cubic_object = Cubic(material=self.material, llw_parameters={'w': 1})

    def __repr__(self):
        """String representation of class"""
        return ''.join([
            f"{self.__class__.__name__}(",
            f"material={self.material!r}, ",
            f"experimental_energies={self.experimental_energies!r}, ",
            f"temperatures={self.temperatures!r})",
        ])

    def get_llw_ratios_plot(self, min_value: float, max_value: float,
                            y_minor: float, y_major: float):
        """Method saves the plot for LLW diagram of energies ratio"""
        gg.get_llw_ratios_plot(
            material=self.material,
            experimental_value=self.experimental_ratio,
            y_limits={
                'max_value': max_value,
                'min_value': min_value,
            },
            y_ticks={
                'y_major': y_major,
                'y_minor': y_minor,
            }
        )

    def get_spectrum_experiment(self, limits: dict, locators: dict,
                                spectrometer='PSI', initial_energy=3):
        """Method saves the plot for experimental spectrum"""
        data = []
        for _temperature in self.temperatures:
            data.append(
                get_data_from_file(
                    os.path.join(
                        DATA_PATHS['experiment'],
                        '_'.join(
                            [spectrometer,
                             self.material.rare_earth,
                             self.material.crystal,
                             f'{initial_energy}meV',
                             f'{_temperature}K.txt'],
                        )
                    )
                )
            )
        gg.get_spectrum_experiment(
            material=self.material,
            temperatures=self.temperatures,
            data=tuple(data),
            scale=Scale(limits=limits, locators=locators),
        )

    def get_cross_points(self):
        """Method returns cross points for experimental and theoretic curves"""
        self.cubic_object.llw_parameters = {'w': 1}
        crosses = self.cubic_object.find_cross(
            experimental_value=self.experimental_ratio,
            experimental_energy=self.experimental_energies[0],
        )
        print('Cross points:')
        for point in crosses:
            print(f'w = {point.w : 6.3f};\tx = {point.x : .3f};')
        return crosses

    def get_spectrum_theory(self, recalculated_crosses, limits: dict, locators: dict, gamma=0.16):
        """Method saves the plot for theoretical spectrum"""
        for point in recalculated_crosses:
            self.cubic_object.llw_parameters = {
                'w': point.w,
                'x': point.x,
            }
            for temperature in self.temperatures:
                self.cubic_object.save_spectra_with_one_temperature(
                    gamma=gamma,
                    temperature=temperature
                )
            spectra = self.cubic_object.save_spectra_with_two_temperatures(
                gamma=gamma,
                temperature_1=self.temperatures[0],
                temperature_2=self.temperatures[1]
            )
            intensities = deepcopy(spectra)
            del intensities['energies']
            gg.get_spectrum_theory(
                material=self.material,
                parameters={
                    'w': point.w,
                    'x': point.x
                },
                temperatures=self.temperatures,
                data=Data(
                    x=spectra['energies'],
                    y_set=intensities,
                    errors=None,
                    legend={
                        self.temperatures[0]: f'{self.temperatures[0]} K',
                        self.temperatures[1]: f'{self.temperatures[1]} K',
                        'diff': f'{self.temperatures[0]} K - {self.temperatures[1]} K'
                    }
                ),
                scale=Scale(
                    limits=limits,
                    locators=locators,
                ),
            )


@get_time_of_execution
def get_fixed_results(crystal: str):
    """Saves the dependence of transition energies, their ratio
    on parameter x to file and its graphs for specified RE ions"""
    rare_earths = ('Tb', 'Tm', 'Er', 'Ho')
    max_energies = (350, 350, 1000, 1500)
    locators = (50, 50, 200, 500)
    for index, value in enumerate(rare_earths):
        material = Material(crystal=crystal, rare_earth=value)
        for w_parameter in (1, -1):
            cubic_object = Cubic(material, llw_parameters={'w': w_parameter})
            cubic_object.save_energy_dat(number_of_intervals=5000)
            cubic_object.get_ratios()
        gg.get_llw_energies_plot(material=material,
                                 max_energy=max_energies[index],
                                 y_major=locators[index],
                                 y_minor=locators[index] // 5)


if __name__ == '__main__':
    # get_fixed_results('YNi2')

    EXPERIMENT = Experiment(
        material=Material(crystal='YNi2', rare_earth='Er'),
        experimental_energies=(0.6, 1.3),
        temperatures=(1.4, 15)
    )

    EXPERIMENT.get_llw_ratios_plot(
        max_value=3,
        min_value=1,
        y_major=0.5,
        y_minor=0.1,
    )
    EXPERIMENT.get_spectrum_experiment(
        limits={
            'x_min': 0,
            'x_max': 2.2,
            'y_min': 0,
            'y_max': 20,
        },
        locators={
            'x_major': 1,
            'x_minor': 0.2,
            'y_major': 5,
            'y_minor': 1,
        }
    )
    RECALCULATED_CROSSES = EXPERIMENT.get_cross_points()
    EXPERIMENT.get_spectrum_theory(
        RECALCULATED_CROSSES,
        limits={
            'x_min': 0,
            'x_max': 2.2,
            'y_min': 0,
            'y_max': 4000,
        },
        locators={
            'x_major': 1,
            'x_minor': 0.2,
            'y_major': 500,
            'y_minor': 100,
        }
    )
