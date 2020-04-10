"""The module contains functions for plotting graphs."""
from collections import OrderedDict
from types import SimpleNamespace
import matplotlib.pyplot as plt
from scripts.common import constants as con
from scripts.common.path_utils import get_paths
from scripts.common.utils import get_time_of_execution, data_popping, OpenedFile, get_default

plt.rcParams['axes.labelpad'] = 0
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['axes.prop_cycle'] = con.DEFAULT_CYCLER
plt.rcParams["errorbar.capsize"] = 1
plt.rcParams['figure.dpi'] = 500
plt.rcParams['figure.figsize'] = [i / 2.54 for i in con.FIG_SIZE]
plt.rcParams['figure.subplot.bottom'] = 0.11
plt.rcParams['figure.subplot.hspace'] = 0
plt.rcParams['figure.subplot.left'] = 0.18
plt.rcParams['figure.subplot.right'] = 0.97
plt.rcParams['figure.subplot.top'] = 0.90
plt.rcParams['figure.subplot.wspace'] = 0
plt.rcParams['font.family'] = con.FONT_FAMILY
plt.rcParams['font.size'] = 14
plt.rcParams['legend.frameon'] = False
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['lines.markersize'] = 1.5
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['mathtext.it'] = con.FONT_FAMILY

for tick in ('xtick', 'ytick'):
    plt.rcParams[f'{tick}.direction'] = 'in'
    plt.rcParams[f'{tick}.major.pad'] = 3
    plt.rcParams[f'{tick}.major.size'] = 6
    plt.rcParams[f'{tick}.major.width'] = 2
    plt.rcParams[f'{tick}.minor.size'] = 3
    plt.rcParams[f'{tick}.minor.width'] = 2

for tick in ('xtick.minor.bottom', 'xtick.minor.top', 'xtick.top',
             'ytick.minor.left', 'ytick.minor.right', 'ytick.right'):
    plt.rcParams[tick] = True


class CustomPlot:
    """Description of Plot object"""

    def __init__(self, data, dpi=300):
        """Initialization of Plot object"""
        self.data = data
        self.dpi = dpi
        self.limits = SimpleNamespace(
            x_min=None,
            x_max=None,
            y_min=None,
            y_max=None
        )
        self.fig = None
        self._ax = None

    def __enter__(self):
        """Method for entrance to context manager"""
        self.fig, self._ax = plt.subplots(dpi=self.dpi)
        return self

    def set_labels(self, x_label='x', y_label='y', title=None):
        """Sets labels of axis and plot"""
        if self._ax:
            self._ax.set_title(title)
            self._ax.set_xlabel(x_label)
            self._ax.set_ylabel(y_label)

    def set_limits(self, x_min=None, x_max=None, y_min=None, y_max=None):
        """Sets limits of x and y intervals"""
        self.limits.x_min = get_default(x_min, min(self.data.x))
        self.limits.x_max = get_default(x_max, max(self.data.x))
        self.limits.y_min = get_default(
            y_min,
            min((min(value) for value in self.data.y_set.values()))
        )
        self.limits.y_max = get_default(
            y_max,
            max((max(value) for value in self.data.y_set.values()))
        )
        if self._ax:
            self._ax.set_xlim(xmin=self.limits.x_min,
                              xmax=self.limits.x_max)
            self._ax.set_ylim(ymin=self.limits.y_min,
                              ymax=self.limits.y_max)

    def set_locators(self, x_major=None, x_minor=None, y_major=None, y_minor=None):
        """Sets major and minor ticks for plot"""
        majors = (
            get_default(x_major, (self.limits.x_max - self.limits.x_min) // 5),
            get_default(y_major, (self.limits.y_max - self.limits.y_min) // 5)
        )
        minors = (
            get_default(x_minor, majors[0] / 5),
            get_default(y_minor, majors[1] / 5)
        )
        if self._ax:
            for i, axis in enumerate((self._ax.xaxis, self._ax.yaxis)):
                axis.set_major_locator(plt.MultipleLocator(majors[i]))
                axis.set_minor_locator(plt.MultipleLocator(minors[i]))

    def make_plot(self, filename=None, mode='plot', form='show', text: con.Text = None):
        """Draws the plot at specified mode"""
        if self.fig and self._ax:
            functions = {
                'plot': self._ax.plot,
                'scatter': self._ax.scatter,
                'errorbar': self._ax.errorbar,
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
            if form == 'show':
                plt.show()
            else:
                self.fig.savefig(f'{filename}.{form}')
                print(f'Graph {filename}.{form} is saved.')

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Method for exit from context manager"""
        if self.fig and self._ax:
            plt.close('all')

    def __repr__(self):
        """String representation of Plot object"""
        return f"{self.__class__.__name__}(data={self.data!r}, dpi={self.dpi!r})"


def get_spectrum_theory(material: con.Material, parameters: dict, temperatures: tuple,
                        data: con.Data = None, scale: con.Scale = None):
    """Returns inelastic neutron scattering spectrum
    that calculated with specified parameters"""
    parameters['T'] = f'{temperatures[0]}-{temperatures[1]}'
    file_name = get_paths(
        data_name='spectra',
        material=material,
        parameters=parameters,
        is_graph=True,
    )
    text = con.Text(x=1, y=2300,
                    string=fr'$W={parameters["w"]: .3f}$' + '\n' + fr'$x = {parameters["x"]: .3f}$')
    with CustomPlot(data=data) as plot:
        plot.set_labels(x_label=con.ENERGY_TRANSFER,
                        y_label=con.SCATTERING)
        plot.set_limits(**scale.limits)
        plot.set_locators(**scale.locators)
        plot.make_plot(filename=file_name,
                       text=text,
                       mode='plot',
                       form=con.FIG_FORMAT)


def get_spectrum_experiment(material: con.Material, temperatures: tuple,
                            data: tuple, scale=None):
    """Returns inelastic neutron scattering spectrum from experiment"""
    parameters = {'T': f'{temperatures[0]}-{temperatures[1]}'}
    file_name = get_paths(
        data_name='spectra',
        material=material,
        parameters=parameters,
        is_graph=True
    )
    data_kwargs = {
        'x': data[0][0],
        'y_set': OrderedDict(),
        'errors': OrderedDict(),
        'legend': OrderedDict(),
    }
    for i in range(2):
        data_kwargs['y_set'][temperatures[i]] = data[i][1]
        data_kwargs['errors'][temperatures[i]] = data[i][2]
        data_kwargs['legend'][temperatures[i]] = f'{temperatures[i]} K'
    data_kwargs['y_set']['diff'] = data[0][1] - data[1][1]
    data_kwargs['errors']['diff'] = data[0][2] + data[1][2]
    data_kwargs['legend']['diff'] = f'{temperatures[0]} K - {temperatures[1]} K'
    with CustomPlot(data=con.Data(**data_kwargs)) as plot:
        plot.set_labels(x_label=con.ENERGY_TRANSFER,
                        y_label=con.SCATTERING)
        plot.set_limits(**scale.limits)
        plot.set_locators(**scale.locators)
        plot.make_plot(filename=file_name,
                       mode='errorbar',
                       form=con.FIG_FORMAT)


@get_time_of_execution
def get_llw_energies_plot(material: con.Material, max_energy, y_major, y_minor):
    """Returns graphs for dependence of transfer energies on CEF parameters"""
    levels_number = 7
    for w_parameter in (1, -1):
        data = {
            'x': [],
            'y_set': OrderedDict(),
            'legend': OrderedDict()
        }
        for level in range(1, levels_number):
            data['y_set'][fr'$E_{level}$'] = []
            data['legend'][fr'$E_{level}$'] = fr'$E_{level}$'

        parameters = {'w': w_parameter}
        energy_file_name = get_paths(
            data_name='energies',
            material=material,
            parameters=parameters,
        )
        with OpenedFile(energy_file_name) as file:
            for line in file:
                array = line.rstrip('\n').split('\t')
                array = [float(value.strip()) for value in array]
                data['x'].append(array[0])
                for level in range(1, levels_number):
                    try:
                        data['y_set'][fr'$E_{level}$'].append(array[level])
                    except IndexError:
                        data['y_set'][fr'$E_{level}$'].append(con.INFINITY)
        data_popping(data, lambda arg: len(arg) == 0)
        graph_file_name = get_paths(
            data_name='energies',
            material=material,
            parameters=parameters,
            is_graph=True
        )
        data = con.Data(x=data['x'],
                        y_set=data['y_set'],
                        legend=data['legend'],
                        errors=None)
        with CustomPlot(data=data) as plot:
            plot.set_labels(x_label=r'$x$',
                            y_label=con.ENERGY_TRANSFER,
                            title=fr'{material.rare_earth}, $W={parameters["w"]}$ meV')
            plot.set_limits(x_min=-1,
                            x_max=1,
                            y_min=-10,
                            y_max=max_energy)
            plot.set_locators(x_major=0.5,
                              x_minor=0.1,
                              y_major=y_major,
                              y_minor=y_minor)
            plot.make_plot(filename=graph_file_name,
                           mode='scatter',
                           form=con.FIG_FORMAT)


@get_time_of_execution
def get_llw_ratios_plot(material: con.Material, experimental_value, y_limits, y_ticks):
    """Returns graphs for dependence of transfer energies on CEF parameters"""
    for w_parameter in (1, -1):
        data = {'x': [],
                'y_set': OrderedDict({'Experiment': []}),
                'legend': OrderedDict({'Experiment': 'Experiment'})
                }
        for name in con.RATIOS_NAMES:
            data['y_set'][name] = []
            data['legend'][name] = name
        parameters = {'w': w_parameter}
        ratio_file_name = get_paths(
            data_name='ratios',
            material=material,
            parameters=parameters,
        )
        with OpenedFile(ratio_file_name) as file:
            for line in file:
                array = line.rstrip('\n').split('\t')
                array = [float(value.strip()) for value in array]
                data['x'].append(array[0])
                data['y_set']['Experiment'].append(experimental_value)
                for level, name in enumerate(con.RATIOS_NAMES):
                    if array[level + 1] == 0:
                        data['y_set'][name].append(con.INFINITY)
                    else:
                        data['y_set'][name].append(array[level + 1])

        data_popping(data, lambda arg: (len(arg) == 0 or
                                        min(arg) > experimental_value or
                                        max(arg) < experimental_value))
        parameters['exp'] = experimental_value
        graph_file_name = get_paths(
            data_name='ratios',
            material=material,
            parameters=parameters,
            is_graph=True,
        )
        data = con.Data(x=data['x'],
                        y_set=data['y_set'],
                        legend=data['legend'],
                        errors=None)
        with CustomPlot(data=data) as plot:
            plot.set_labels(x_label=r'$x$',
                            y_label='Energy Transfer Ratio',
                            title=fr'$W={parameters["w"]}$')
            plot.set_limits(x_min=-1,
                            x_max=1,
                            y_min=y_limits['min_value'],
                            y_max=y_limits['max_value'])
            plot.set_locators(x_major=0.5,
                              x_minor=0.1,
                              y_major=y_ticks['y_major'],
                              y_minor=y_ticks['y_minor'])
            plot.make_plot(filename=graph_file_name,
                           mode='scatter',
                           form=con.FIG_FORMAT)


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
    DATA = con.Data(x=X_ARRAY, y_set=Y_ARRAY, legend=LEGEND, errors=None)
    with CustomPlot(data=DATA, dpi=100) as custom_plot:
        custom_plot.set_limits(x_min=2, x_max=8, y_min=0, y_max=5000)
        custom_plot.set_labels(x_label='x_test', y_label='y_test', title='test')
        custom_plot.set_locators()
        custom_plot.make_plot()
