"""The module contains functions for plotting graphs."""
from collections import OrderedDict
from types import SimpleNamespace
import matplotlib.pyplot as plt
from scripts.common import constants as con
from scripts.common.path_utils import get_paths
from scripts.common.utils import get_time_of_execution, data_popping, OpenedFile
from scripts.common.named_tuples import Data, Text, Scale, Material

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
# plt.rc('text', usetex=True)

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
    """Description of plot object"""
    def __init__(self, data, dpi=300):
        """Initialization of Plot object"""
        self.data = data
        self.dpi = dpi
        self.limits = SimpleNamespace(x_min=None,
                                      x_max=None,
                                      y_min=None,
                                      y_max=None)
        self.fig, self._ax = plt.subplots(dpi=self.dpi)

    def set_labels(self, x_label='x', y_label='y', title=None):
        """Sets labels of axis and plot"""
        self._ax.set_title(title)
        self._ax.set_xlabel(x_label)
        self._ax.set_ylabel(y_label)

    def set_limits(self, x_min=None, x_max=None, y_min=None, y_max=None):
        """Sets limits of x and y intervals"""
        self.limits.x_min = min(self.data.x) if x_min is None else x_min
        self.limits.x_max = max(self.data.x) if x_max is None else x_max
        self.limits.y_min = (min([min(value) for value in self.data.y_set.values()])
                             if y_min is None else y_min)
        self.limits.y_max = (max([max(value) for value in self.data.y_set.values()])
                             if y_max is None else y_max)
        self._ax.set_xlim(xmin=self.limits.x_min, xmax=self.limits.x_max)
        self._ax.set_ylim(ymin=self.limits.y_min, ymax=self.limits.y_max)

    def set_locators(self, x_major=None, x_minor=None, y_major=None, y_minor=None):
        """Sets major and minor ticks for plot"""
        x_major = (self.limits.x_max - self.limits.x_min) // 5 if x_major is None else x_major
        y_major = (self.limits.y_max - self.limits.y_min) // 5 if y_major is None else y_major
        x_minor = x_major / 5 if x_minor is None else x_minor
        y_minor = y_major / 5 if y_minor is None else y_minor
        self._ax.xaxis.set_major_locator(plt.MultipleLocator(x_major))
        self._ax.yaxis.set_major_locator(plt.MultipleLocator(y_major))
        self._ax.xaxis.set_minor_locator(plt.MultipleLocator(x_minor))
        self._ax.yaxis.set_minor_locator(plt.MultipleLocator(y_minor))

    def make_plot(self, filename=None, mode='line', form='show', text: Text = None):
        """Draws the plot at specified mode"""
        for key, y_data in self.data.y_set.items():
            if mode == 'line':
                self._ax.plot(self.data.x, y_data, label=self.data.legend[key])
            elif mode == 'scatter':
                self._ax.scatter(self.data.x, y_data, label=self.data.legend[key])
            elif mode == 'errorbar':
                self._ax.errorbar(self.data.x, y_data, yerr=self.data.errors[key], fmt='o',
                                  elinewidth=1, label=self.data.legend[key])
        self._ax.legend()
        if text:
            plt.text(x=text.x, y=text.y, s=text.string)
        if form == 'show':
            plt.show()
        else:
            self.fig.savefig(f'{filename}.{form}')
        plt.close('all')

    def __repr__(self):
        """String representation of Plot object"""
        return f"{self.__class__.__name__}(data={self.data!r}, dpi={self.dpi!r})"


def get_spectrum_theory(material: Material, parameters: dict, temperatures: tuple,
                        data: Data = None, scale: Scale = None):
    """Returns inelastic neutron scattering spectrum that calculated with specified parameters"""
    parameters['T'] = f'{temperatures[0]}-{temperatures[1]}'
    file_name = get_paths(con.PATH_TO_SPECTRA_GRAPHS, 'spectrum', format_name='',
                          material=material, parameters=parameters)
    text = Text(x=1, y=2300,
                string=fr'$W={parameters["w"]: .3f}$' + '\n' + fr'$x = {parameters["x"]: .3f}$')
    plot = CustomPlot(data=data)
    plot.set_labels(x_label=con.ENERGY_TRANSFER, y_label=con.SCATTERING)
    plot.set_limits(**scale.limits)
    plot.set_locators(**scale.locators)
    plot.make_plot(filename=file_name, text=text, mode='line', form=con.FIG_FORMAT)
    print(f'Graph {file_name} saved')


def get_spectrum_experiment(material: Material, temperatures: tuple,
                            data: tuple, scale=None):
    """Returns inelastic neutron scattering spectrum from experiment"""
    parameters = {'T': f'{temperatures[0]}-{temperatures[1]}'}
    file_name = get_paths(con.PATH_TO_SPECTRA_GRAPHS, 'spectrum', format_name='',
                          material=material, parameters=parameters)
    data = Data(x=data[0][0],
                y_set={temperatures[0]: data[0][1],
                       temperatures[1]: data[1][1],
                       'diff': data[0][1] - data[1][1],
                       },
                errors={temperatures[0]: data[0][2],
                        temperatures[1]: data[1][2],
                        'diff': data[0][2] + data[1][2],
                        },
                legend={temperatures[0]: f'{temperatures[0]} K',
                        temperatures[1]: f'{temperatures[1]} K',
                        'diff': f'{temperatures[0]} K - {temperatures[1]} K'}
                )
    plot = CustomPlot(data=data)
    plot.set_labels(x_label=con.ENERGY_TRANSFER, y_label=con.SCATTERING)
    plot.set_limits(**scale.limits)
    plot.set_locators(**scale.locators)
    plot.make_plot(filename=file_name, mode='errorbar', form=con.FIG_FORMAT)


@get_time_of_execution
def get_llw_energies_plot(material: Material, max_energy, y_major, y_minor):
    """Returns graphs for dependence of transfer energies on CEF parameters"""
    levels_number = 7
    for w_parameter in (1, -1):
        data = {'x': [], 'y_set': OrderedDict(), 'legend': OrderedDict()}
        for level in range(1, levels_number):
            data['y_set'][fr'$E_{level}$'] = []
            data['legend'][fr'$E_{level}$'] = fr'$E_{level}$'

        parameters = {'w': w_parameter}
        energy_file_name = get_paths(con.PATH_TO_ENERGY_DATAFILES, 'energy',
                                     material=material, parameters=parameters)
        with OpenedFile(energy_file_name) as file:
            lines = file.readlines()
            for line in lines:
                line = line.rstrip('\n').split('\t')
                line = [float(value.lstrip(' ')) for value in line]
                data['x'].append(line[0])
                for level in range(1, levels_number):
                    try:
                        data['y_set'][fr'$E_{level}$'].append(line[level])
                    except IndexError:
                        data['y_set'][fr'$E_{level}$'].append(con.INFINITY)

        data_popping(data, lambda arg: len(arg) == 0)
        graph_file_name = get_paths(con.PATH_TO_ENERGY_GRAPHS, 'energies', format_name='',
                                    material=material, parameters=parameters)
        data = Data(x=data['x'], y_set=data['y_set'],
                    legend=data['legend'], errors=None)
        plot = CustomPlot(data=data)
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
        plot.make_plot(filename=graph_file_name, mode='errorbar', form=con.FIG_FORMAT)


def get_all_llw_energies():
    """Returns graphs for dependence of transfer energies
    on CEF parameters for 4 RE ions"""
    rare_earths = ('Tb', 'Tm', 'Er', 'Ho')
    max_energies = (350, 350, 1000, 1500)
    locators = (50, 50, 200, 500)
    for index, value in enumerate(rare_earths):
        get_llw_energies_plot(material=Material(crystal='YNi2', rare_earth=value),
                              max_energy=max_energies[index],
                              y_major=locators[index],
                              y_minor=locators[index] // 5)


@get_time_of_execution
def get_llw_ratios_plot(material: Material, experimental_value, y_limits, y_ticks):
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
        ratio_file_name = get_paths(con.PATH_TO_RATIO_DATAFILES, 'ratio',
                                    material=material, parameters=parameters)
        with OpenedFile(ratio_file_name) as file:
            lines = file.readlines()
            for line in lines:
                line = line.rstrip('\n').split('\t')
                line = [float(value.lstrip(' ')) for value in line]
                data['x'].append(line[0])
                data['y_set']['Experiment'].append(experimental_value)
                for level, name in enumerate(con.RATIOS_NAMES):
                    if line[level + 1] == 0:
                        data['y_set'][name].append(con.INFINITY)
                    else:
                        data['y_set'][name].append(line[level + 1])

        data_popping(data, lambda arg: (len(arg) == 0 or
                                        min(arg) > experimental_value or
                                        max(arg) < experimental_value))
        parameters['exp'] = experimental_value
        graph_file_name = get_paths(con.PATH_TO_RATIO_GRAPHS, 'ratios', format_name='',
                                    material=material, parameters=parameters)
        data = Data(x=data['x'], y_set=data['y_set'],
                    legend=data['legend'], errors=None)
        plot = CustomPlot(data=data)
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
        plot.make_plot(filename=graph_file_name, mode='scatter', form=con.FIG_FORMAT)


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
    DATA = Data(x=X_ARRAY, y_set=Y_ARRAY, legend=LEGEND, errors=None)
    TEST_PLOT = CustomPlot(data=DATA, dpi=100)
    TEST_PLOT.set_limits(x_min=2, x_max=8, y_min=0, y_max=5000)
    TEST_PLOT.set_labels(x_label='x_test', y_label='y_test', title='test')
    TEST_PLOT.set_locators()
    TEST_PLOT.make_plot()

    # get_all_llw_energies()
