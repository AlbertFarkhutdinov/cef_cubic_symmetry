"""The module contains functions for plotting graphs."""


from collections import OrderedDict

import matplotlib.pyplot as plt
from cycler import cycler

from common import constants as con
from common import utils as ut
from common.path_utils import get_paths
from scripts.cubic_cef_object import Cubic


def _set_plot_parameters():
    """Setting of rcParams"""
    custom_parameters = ut.get_json_object('plot_parameters.json')
    custom_parameters[
        'axes.prop_cycle'
    ] = (cycler(color=custom_parameters['axes.prop_cycle']['color']) +
         cycler(linestyle=custom_parameters['axes.prop_cycle']['linestyle']))
    custom_parameters['figure.figsize'] = [i / 2.54 for i in (10, 10)]
    for _key, _value in custom_parameters.items():
        plt.rcParams[_key] = _value
    tick_parameters = {
        'direction': 'in',
        'major.pad': 3,
        'major.size': 6,
        'major.width': 2,
        'minor.size': 3,
        'minor.width': 2,
    }
    for tick in ('xtick', 'ytick'):
        for _key, _value in tick_parameters.items():
            plt.rcParams[f'{tick}.{_key}'] = _value


class CustomPlot:
    """Description of Plot object"""

    def __init__(self, data, dpi=300):
        """Initialization of Plot object"""
        self.data = data
        self.dpi = dpi
        self.limits = {
            _key: None
            for _key
            in ('x_min', 'x_max', 'y_min', 'y_max')
        }
        self.fig = None
        self._ax = None

    def __enter__(self):
        """Method for entrance to context manager"""
        _set_plot_parameters()
        self.fig, self._ax = plt.subplots(dpi=self.dpi)
        return self

    def set_labels(self,
                   xlabel='x',
                   ylabel='y',
                   title=None):
        """Sets labels of axis and plot"""
        args = locals()
        del args['self']
        if self._ax:
            for _key, _value in args.items():
                self._ax.__getattribute__(f'set_{_key}')(_value)

    def set_limits(self,
                   x_min=None,
                   x_max=None,
                   y_min=None,
                   y_max=None):
        """Sets limits of x and y intervals"""
        y_set = self.data.y_set.values()
        _limits = {
            'x_min': (x_min, min(self.data.x)),
            'x_max': (x_max, max(self.data.x)),
            'y_min': (y_min, min((min(value) for value in y_set))),
            'y_max': (y_max, max((max(value) for value in y_set))),
        }
        for _key, _value in _limits.items():
            self.limits[_key] = ut.get_default(*_value)
        if self._ax:
            for axis in ('x', 'y'):
                axis_limits = {
                    f'{axis}{lim}': self.limits[f'{axis}_{lim}']
                    for lim in ('min', 'max')
                }
                self._ax.__getattribute__(f'set_{axis}lim')(**axis_limits)

    def set_locators(self,
                     x_major=None,
                     x_minor=None,
                     y_major=None,
                     y_minor=None):
        """Sets major and minor ticks for plot"""
        majors = (
            ut.get_default(
                x_major,
                (self.limits['x_max'] - self.limits['x_min']) // 5,
            ),
            ut.get_default(
                y_major,
                (self.limits['y_max'] - self.limits['y_min']) // 5,
            )
        )
        minors = (
            ut.get_default(x_minor, majors[0] / 5),
            ut.get_default(y_minor, majors[1] / 5)
        )
        if self._ax:
            for i, axis in enumerate((self._ax.xaxis, self._ax.yaxis)):
                axis.set_major_locator(plt.MultipleLocator(majors[i]))
                axis.set_minor_locator(plt.MultipleLocator(minors[i]))

    def make_plot(self,
                  mode='plot',
                  text: con.Text = None):
        """Draws the plot at specified mode"""
        if self.fig and self._ax:
            functions = {
                key: self._ax.__getattribute__(key)
                for key in ('plot', 'scatter', 'errorbar')
            }
            for key, y_data in self.data.y_set.items():
                args = (self.data.x, y_data)
                kwargs = {
                    'label': self.data.legend[key]
                }
                if mode == 'errorbar':
                    kwargs['yerr'] = self.data.errors[key]
                    kwargs['fmt'] = 'o'
                    kwargs['elinewidth'] = 1
                functions[mode](*args, **kwargs)
            self._ax.legend()
            if text:
                plt.text(x=text.x, y=text.y, s=text.string)

    def save_or_show(self, filename=None, form=None):
        """Saves or shows the plot"""
        if self.fig and self._ax:
            if form:
                self.fig.savefig(f'{filename}.{form}')
                print(f'Graph {filename}.{form} is saved.')
            else:
                plt.show()

    def save_in_two_forms(self,
                          filename: str,
                          form_1='png',
                          form_2='eps'):
        """Saves or shows the plot"""
        self.save_or_show(filename=filename, form=form_1)
        self.save_or_show(filename=filename, form=form_2)

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Method for exit from context manager"""
        if self.fig and self._ax:
            plt.close('all')

    def __repr__(self):
        """Method returns string representation of the Plot object."""
        return ut.get_repr(self, 'data', 'dpi')


class CubicPlot(CustomPlot):
    """Description of CubicPlot object"""

    def __init__(self,
                 data,
                 material:
                 con.Material,
                 dpi=300):
        """Initialization of Plot object"""
        super().__init__(data=data, dpi=dpi)
        self.material = material

    def __enter__(self):
        """Method for entrance to context manager"""
        super().__enter__()
        return self

    def get_graph_file_name(self,
                            data_name: str,
                            parameters: dict = None):
        """Returns path for plot saving"""
        return get_paths(
            data_name=data_name,
            material=self.material,
            parameters=parameters,
            is_graph=True,
        )

    def __repr__(self):
        """Method returns string representation of the Plot object."""
        return ut.get_repr(self, 'data', 'material', 'dpi')


@ut.get_time_of_execution
def get_llw_plot(material: con.Material,
                 y_max,
                 y_major,
                 y_minor,
                 choice=0):
    """Returns graphs for dependence of transition energies
    or intensities on CEF parameters"""
    data_name = 'energies' if choice == 0 else 'intensities'
    for w_parameter in (1, -1):
        data = {
            'x': [],
            'y_set': OrderedDict(),
            'legend': OrderedDict()
        }
        for level in range(1, 7):
            data['y_set'][ut.get_label(level, choice)] = []
            data['legend'][ut.get_label(level, choice)] = ut.get_label(
                level, choice
            )

        parameters = {'w': w_parameter}
        peak_file_name = get_paths(
            data_name=data_name,
            material=material,
            parameters=parameters,
        )
        with ut.UTF8File(peak_file_name) as file:
            for line in file:
                array = line.rstrip('\n').split('\t')
                array = [float(value.strip()) for value in array]
                data['x'].append(array[0])
                for level in range(1, 7):
                    try:
                        data['y_set'][
                            ut.get_label(level, choice)
                        ].append(
                            con.INFINITY
                            if ((choice != 0) and (level == 1))
                            else array[level]
                        )
                    except IndexError:
                        data['y_set'][
                            ut.get_label(level, choice)
                        ].append(con.INFINITY)
        ut.data_popping(data, lambda arg: len(arg) == 0)
        data = con.Data(
            x=data['x'],
            y_set=data['y_set'],
            legend=data['legend'],
            errors=None,
        )
        with CubicPlot(data=data, material=material) as plot:
            plot.set_labels(
                xlabel=r'$x$',
                ylabel=(
                    con.ENERGY_TRANSFER
                    if choice == 0
                    else con.TRANSITION_INTENSITY
                ),
                title=fr'{material.rare_earth}, $W={parameters["w"]}$ мэВ',
            )
            if y_max:
                plot.set_limits(
                    x_min=-1,
                    x_max=1,
                    y_min=-10 if choice == 0 else 0,
                    y_max=y_max,
                )
            if y_major and y_minor:
                plot.set_locators(
                    x_major=0.5,
                    x_minor=0.1,
                    y_major=y_major,
                    y_minor=y_minor,
                )
            plot.make_plot(mode='scatter')
            plot.save_in_two_forms(
                filename=plot.get_graph_file_name(
                    data_name=data_name,
                    parameters=parameters,
                )
            )


@ut.get_time_of_execution
def get_llw_ratios_plot(material: con.Material,
                        experimental_value,
                        limits: dict,
                        ticks: dict,
                        choice=0):
    """Returns graphs for dependence of transition energies or intensities
     ratios on CEF parameters"""
    data_name = 'ratios_energies' if choice == 0 else 'ratios_intensities'
    for w_parameter in (1, -1):
        data = {
            'x': [],
            'y_set': OrderedDict({'Experiment': []}),
            'legend': OrderedDict({'Experiment': 'Experiment'})
        }
        for name in ut.get_ratios_names(choice):
            data['y_set'][name] = []
            data['legend'][name] = name
        parameters = {'w': w_parameter}
        ratio_file_name = get_paths(
            data_name=data_name,
            material=material,
            parameters=parameters,
        )
        with ut.UTF8File(ratio_file_name) as file:
            for line in file:
                array = line.rstrip('\n').split('\t')
                array = [float(value.strip()) for value in array]
                data['x'].append(array[0])
                data['y_set']['Experiment'].append(experimental_value)
                for level, name in enumerate(ut.get_ratios_names(choice)):
                    if array[level + 1] == 0:
                        data['y_set'][name].append(con.INFINITY)
                    else:
                        data['y_set'][name].append(array[level + 1])
        ut.data_popping(
            data,
            lambda arg: (
                len(arg) == 0 or
                min(arg) > experimental_value or
                max(arg) < experimental_value
            )
        )
        parameters['exp'] = experimental_value
        data = con.Data(
            x=data['x'],
            y_set=data['y_set'],
            legend=data['legend'],
            errors=None,
        )
        with CubicPlot(data=data, material=material) as plot:
            plot.set_labels(
                xlabel=r'$x$',
                ylabel=(
                    con.ENERGY_TRANSFER_RATIO
                    if choice == 0
                    else con.TRANSITION_INTENSITY_RATIO
                ),
                title=fr'{material.rare_earth}, $W={parameters["w"]}$ мэВ',
            )
            plot.set_limits(**limits)
            plot.set_locators(**ticks)
            plot.make_plot(mode='scatter')
            plot.save_in_two_forms(
                filename=plot.get_graph_file_name(
                    data_name=data_name,
                    parameters=parameters,
                )
            )


def get_spectrum_theory(material: con.Material,
                        parameters: dict,
                        data: con.Data = None,
                        scale: con.Scale = None):
    """Returns inelastic neutron scattering spectrum
    that calculated with specified parameters"""
    with CubicPlot(data=data, material=material) as plot:
        plot.set_labels(**con.SPECTRUM_LABELS)
        plot.set_limits(**scale.limits)
        plot.set_locators(**scale.locators)
        plot.make_plot(mode='plot')
        plot.save_in_two_forms(
            filename=plot.get_graph_file_name(
                data_name='spectra',
                parameters=parameters,
            )
        )


def get_spectrum_experiment(material: con.Material,
                            temperatures: tuple,
                            data: tuple,
                            parameters: dict,
                            scale=None):
    """Returns inelastic neutron scattering spectrum from experiment"""
    data_kwargs = {
        'x': data[0]['x'],
        'y_set': OrderedDict(),
        'errors': OrderedDict(),
        'legend': OrderedDict(),
    }
    for i, temperature in enumerate(temperatures):
        data_kwargs['y_set'][temperature] = data[i]['y']
        data_kwargs['errors'][temperature] = data[i]['errors']
        data_kwargs['legend'][temperature] = f'{temperature} K'
    with CubicPlot(data=con.Data(**data_kwargs),
                   material=material) as plot:
        plot.set_labels(**con.SPECTRUM_LABELS)
        plot.set_limits(**scale.limits)
        plot.set_locators(**scale.locators)
        plot.make_plot(mode='errorbar')
        plot.save_in_two_forms(
            filename=plot.get_graph_file_name(
                data_name='spectra',
                parameters=parameters,
            )
        )


def get_intensity_on_temperature(
        material: con.Material,
        crosses: con.CrossPoint,
        y_max: float,
):
    """
    Returns graphs for dependence of transition intensities on temperature.

    """
    data_kwargs = {
        'x': [],
        'y_set': OrderedDict(),
        'errors': OrderedDict(),
        'legend': OrderedDict(),
    }
    data_name = 'intensities_on_temperature'
    for cross_number, cross in enumerate(crosses):
        llw = {'w': cross.w, 'x': cross.x}
        label = f'$W = {cross.w:.3f}, x = {cross.x:.3f}$'
        data_kwargs['y_set'][label] = []
        data_kwargs['errors'][label] = None
        data_kwargs['legend'][label] = label
        cubic_object = Cubic(
            material=material,
            llw_parameters=llw,
        )
        cubic_object.save_intensities()
        file_name = get_paths(
            data_name=data_name,
            material=material,
            parameters=llw,
        )
        with ut.UTF8File(file_name) as file:
            for line in file:
                row = line.rstrip('\n').split('\t')
                row = [float(value) for value in row]
                if cross_number == 0:
                    data_kwargs['x'].append(row[0])
                data_kwargs['y_set'][label].append(row[1] / row[2])
    with CubicPlot(data=con.Data(**data_kwargs), material=material) as plot:
        plot.set_labels(
            xlabel='Temperature, K',
            ylabel=con.TRANSITION_INTENSITY_RATIO,
        )
        plot.set_limits(
            y_min=0,
            y_max=y_max,
        )
        plot.set_locators()
        plot.make_plot(mode='scatter')
        plot.save_in_two_forms(
            filename=plot.get_graph_file_name(
                data_name=data_name,
            )
        )


if __name__ == '__main__':
    X_ARRAY = list(range(11))
    Y_ARRAY = {
        '4': [x ** 4 for x in X_ARRAY],
        '3': [x ** 3 for x in X_ARRAY],
        '2': [x ** 2 for x in X_ARRAY],
        '1': [x ** 1 for x in X_ARRAY],
    }
    LEGEND = {
        '4': '$x^4$',
        '3': '$x^3$',
        '2': '$x^2$',
        '1': '$x^1$',
    }
    DATA = con.Data(
        x=X_ARRAY,
        y_set=Y_ARRAY,
        legend=LEGEND,
        errors=None,
    )
    with CustomPlot(data=DATA, dpi=100) as custom_plot:
        custom_plot.set_limits(
            x_min=2,
            x_max=8,
            y_min=0,
            y_max=5000,
        )
        custom_plot.set_labels(
            xlabel='x_test',
            ylabel='y_test',
            title='test',
        )
        custom_plot.set_locators()
        custom_plot.make_plot()
        custom_plot.save_or_show()
