"""The module contains Experiment class."""

from copy import deepcopy
from pathlib import Path

from cef_cubic_symmetry.auxiliary import constants as con
from cef_cubic_symmetry.auxiliary.paths import DATA_PATHS
from cef_cubic_symmetry.auxiliary.utils import get_repr
from cef_cubic_symmetry.core.sample import Sample
from cef_cubic_symmetry.fitting.fitting_procedures import get_data_from_file
from cef_cubic_symmetry.scripts import plot_objects as gg
from cef_cubic_symmetry.scripts.cubic_cef_object import Cubic


class Experiment:
    """Class contains experimental parameters."""

    def __init__(
        self,
        material: Sample,
        experimental_energies: tuple,
        temperatures: tuple,
    ) -> None:
        """Initialize class Experiment."""
        self.material = material
        self.experimental_energies = experimental_energies
        self.experimental_ratio = (
            experimental_energies[1] / experimental_energies[0]
        )
        self.temperatures = temperatures
        self.cubic_object = Cubic(
            material=self.material,
            llw_parameters={'w': 1},
        )

    def __repr__(self) -> str:
        """Return string representation of the Experiment object."""
        return get_repr(
            self,
            'material',
            'experimental_energies',
            'temperatures',
        )

    def get_llw_ratios_plot(
        self,
        limits: dict,
        ticks: dict,
    ) -> None:
        """Save the plot for LLW diagram of energies ratio."""
        gg.get_llw_ratios_plot(
            material=self.material,
            experimental_value=self.experimental_ratio,
            limits=limits,
            ticks=ticks,
        )

    def get_spectrum_experiment(
        self,
        limits: dict,
        locators: dict,
        spectrometer: str,
        initial_energy: float,
    ) -> tuple:
        """Save the plot for experimental spectrum."""
        data, temperatures = self._get_spectrum_experiment(
            spectrometer=spectrometer,
            initial_energy=initial_energy,
        )
        gg.get_spectrum_experiment(
            material=self.material,
            temperatures=tuple(temperatures),
            data=tuple(data),
            parameters={
                'setup': spectrometer,
            },
            scale=con.Scale(limits=limits, locators=locators),
        )
        return data, temperatures

    def get_spectrum_differences(
        self,
        limits: dict,
        locators: dict,
        spectrometer: str,
        initial_energy: float,
    ) -> None:
        """Save the plot for experimental spectrum."""
        diff_data, differences = self._get_spectrum_differences(
            spectrometer=spectrometer,
            initial_energy=initial_energy,
        )
        gg.get_spectrum_experiment(
            material=self.material,
            temperatures=tuple(differences),
            data=tuple(diff_data),
            parameters={
                'setup': spectrometer,
                'T': differences[0],
            },
            scale=con.Scale(limits=limits, locators=locators),
        )

    def get_cross_points(self) -> list:
        """Return cross points for experimental and theoretic curves."""
        self.cubic_object.llw_parameters = {'w': 1}
        crosses = self.cubic_object.find_cross(
            experimental_value=self.experimental_ratio,
            experimental_energy=self.experimental_energies[0],
        )
        print('Cross points:')
        for point in crosses:
            print(
                f'w = {point.w: 6.3f};\tx = {point.x: .3f}; '
                f'Ratio: {point.ratio_name}',
            )
        return crosses

    def get_intensity_on_temperature(self, crosses, y_max: float) -> None:
        """
        Save the plot for dependence of transition intensities on temperature.

        Parameters
        ----------
        crosses
            Cross points.
        y_max : float
            Max y value.

        """
        gg.get_intensity_on_temperature(
            material=self.material,
            crosses=crosses,
            y_max=y_max,
        )

    def get_spectrum_theory(
        self,
        recalculated_crosses,
        limits: dict,
        locators: dict,
        gamma=0.16,
    ) -> None:
        """Save the plot for theoretical spectrum."""
        for point in recalculated_crosses:
            self.cubic_object.llw_parameters = {
                'w': point.w,
                'x': point.x,
            }
            for temperature in self.temperatures:
                self.cubic_object.save_spectra_with_one_temperature(
                    gamma=gamma,
                    temperature=temperature,
                )
            spectra = self.cubic_object.save_spectra_with_many_temperatures(
                gamma=gamma,
                temperatures=self.temperatures,
            )
            intensities = deepcopy(spectra)
            del intensities['energies']
            gg.get_spectrum_theory(
                material=self.material,
                parameters={
                    'w': point.w,
                    'x': point.x,
                },
                data=con.Data(
                    x=spectra['energies'],
                    y_set=intensities,
                    errors=None,
                    legend={
                        temperature: f'{temperature} K'
                        for temperature in self.temperatures
                    },
                ),
                scale=con.Scale(
                    limits=limits,
                    locators=locators,
                ),
            )

    def _get_spectrum_experiment(
        self,
        spectrometer: str,
        initial_energy: float,
    ) -> tuple:
        """Return data for experimental spectra."""
        data = []
        temperatures = []
        for _temperature in self.temperatures:
            try:
                data.append(
                    get_data_from_file(
                        Path(DATA_PATHS['experiment']).joinpath(
                            f'{self.material.crystal.name}_'
                            f'{self.material.rare_earth.identifier}',
                            '_'.join(
                                [spectrometer,
                                 self.material.rare_earth.identifier,
                                 self.material.crystal.name,
                                 f'{initial_energy}meV',
                                 f'{_temperature}K.dat'],
                            ),
                        ),
                    ),
                )
                temperatures.append(_temperature)
            except FileNotFoundError:
                pass
        return data, temperatures

    def _get_spectrum_differences(
        self,
        spectrometer: str,
        initial_energy: float,
    ) -> tuple:
        """Return data for experimental spectra differences."""
        data, temperatures = self._get_spectrum_experiment(
            spectrometer=spectrometer,
            initial_energy=initial_energy,
        )
        diff_data = []
        differences = []
        for index, _temperature in enumerate(temperatures[1:]):
            result = {
                'x': [],
                'y': [],
                'errors': [],
            }
            for _index, _ in enumerate(data[0]['x']):
                result['x'].append(data[0]['x'][_index])
                result['y'].append(
                    data[0]['y'][_index] - data[index + 1]['y'][_index],
                )
                result['errors'].append(
                    data[0]['errors'][_index]
                    + data[index + 1]['errors'][_index],
                )
            differences.append(f'{temperatures[0]} K - {_temperature}')
            diff_data.append(result)
        return diff_data, differences
