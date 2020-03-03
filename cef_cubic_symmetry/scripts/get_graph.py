"""The module contains functions for plotting graphs."""
from collections import namedtuple, OrderedDict
import matplotlib.pyplot as plt
from scripts.common import constants as con
from scripts.common.path_utils import check_path, get_paths
from scripts.common.utils import get_time_of_execution, data_popping

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


Data = namedtuple('Data', ['x', 'y_set', 'errors', 'legend'])
Labels = namedtuple('Labels', ['x', 'y'])
Limits = namedtuple('Limits', ['x_min', 'x_max', 'y_min', 'y_max'])
Locators = namedtuple('Locators', ['x_major', 'x_minor', 'y_major', 'y_minor'])
Text = namedtuple('Text', ['x', 'y', 'string'])
GraphInfo = namedtuple('GraphInfo', ['filename', 'title', 'mode', 'form'])
Scale = namedtuple('Scale', ['limits', 'locators'])


def get_plot(info=GraphInfo(filename='graph', title=None, mode='line', form='show'),
             labels=Labels(x=r'$x$', y=r'$y$'), scale=None, data=None,
             text=Text(x=0, y=0, string='')):
    """Returns graph with specified parameters"""
    if scale is None and data is not None:
        limits = Limits(x_min=min(data.x),
                        x_max=max(data.x),
                        y_min=min([min(value) for value in data.y_set.values()]),
                        y_max=max([max(value) for value in data.y_set.values()])
                        )
        locators = Locators(x_major=(limits.x_max - limits.x_min) // 5,
                            x_minor=(limits.x_max - limits.x_min) // 5 / 5,
                            y_major=(limits.y_max - limits.y_min) // 5,
                            y_minor=(limits.y_max - limits.y_min) // 5 / 5)
        scale = Scale(limits=limits, locators=locators)
    dpi = 100 if info.form == 'show' else 500
    fig, _ax = plt.subplots(dpi=dpi)
    _ax.set_title(info.title)
    _ax.set_xlabel(labels.x)
    _ax.set_ylabel(labels.y)
    _ax.xaxis.set_major_locator(plt.MultipleLocator(scale.locators.x_major))
    _ax.xaxis.set_minor_locator(plt.MultipleLocator(scale.locators.x_minor))
    _ax.yaxis.set_major_locator(plt.MultipleLocator(scale.locators.y_major))
    _ax.yaxis.set_minor_locator(plt.MultipleLocator(scale.locators.y_minor))
    _ax.set_xlim(xmin=scale.limits.x_min, xmax=scale.limits.x_max)
    _ax.set_ylim(ymin=scale.limits.y_min, ymax=scale.limits.y_max)
    for key, y_data in data.y_set.items():
        if info.mode == 'line':
            _ax.plot(data.x, y_data, label=data.legend[key])
        elif info.mode == 'scatter':
            _ax.scatter(data.x, y_data, label=data.legend[key])
        elif info.mode == 'errorbar':
            _ax.errorbar(data.x, y_data, yerr=data.errors[key], fmt='o',
                         elinewidth=1, label=data.legend[key])
    _ax.legend()
    plt.text(x=text.x, y=text.y, s=text.string)
    if info.form == 'show':
        plt.show()
    else:
        fig.savefig(f'{info.filename}.{info.form}')
    plt.close('all')


def get_spectrum_theory(material: dict, parameters: dict, temperatures: tuple,
                        data=None, scale=None):
    """Returns inelastic neutron scattering spectrum that calculated with specified parameters"""
    parameters['T'] = f'{temperatures[0]}-{temperatures[1]}'
    file_name = get_paths(con.PATH_TO_SPECTRA_GRAPHS, 'spectrum', format_name='',
                          material=material, parameters=parameters)
    check_path(file_name)
    text = Text(x=1, y=2300,
                string=fr'$W={parameters["w"]: .3f}$' + '\n' + fr'$x = {parameters["x"]: .3f}$')
    get_plot(data=data,
             scale=scale, text=text,
             labels=Labels(x=con.ENERGY_TRANSFER, y=con.SCATTERING),
             info=GraphInfo(filename=file_name, 
                            title=None, 
                            form=con.FIG_FORMAT, 
                            mode='line')
             )


def get_spectrum_experiment(material: dict, temperatures: tuple,
                            data: tuple, scale=None):
    """Returns inelastic neutron scattering spectrum from experiment"""
    parameters = {'T': f'{temperatures[0]}-{temperatures[1]}'}
    file_name = get_paths(con.PATH_TO_SPECTRA_GRAPHS, 'spectrum', format_name='',
                          material=material, parameters=parameters)
    check_path(file_name)
    get_plot(data=Data(x=data[0][0],
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
                       ),
             scale=scale,
             labels=Labels(x=con.ENERGY_TRANSFER, y=con.SCATTERING),
             info=GraphInfo(filename=file_name, 
                            title=None, 
                            form=con.FIG_FORMAT, 
                            mode='errorbar')
             )


@get_time_of_execution
def get_llw_energies_plot(material: dict, max_energy, y_major, y_minor):
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
        with open(energy_file_name, 'r', encoding=con.ENCODING) as file:
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
        check_path(graph_file_name)
        data = Data(x=data['x'], y_set=data['y_set'],
                    legend=data['legend'], errors=None)
        get_plot(data=data,
                 scale=Scale(limits=Limits(x_min=-1, x_max=1, y_min=-10, y_max=max_energy),
                             locators=Locators(x_major=0.5, x_minor=0.1,
                                               y_major=y_major, y_minor=y_minor)),
                 labels=Labels(x=r'$x$', y=con.ENERGY_TRANSFER),
                 info=GraphInfo(filename=graph_file_name,
                                title=fr'{material["rare_earth"]}, $W={parameters["w"]}$ meV',
                                form=con.FIG_FORMAT, 
                                mode='scatter')
                 )


def get_all_llw_energies():
    """Returns graphs for dependence of transfer energies
    on CEF parameters for 4 RE ions"""
    rare_earths = ('Tb', 'Tm', 'Er', 'Ho')
    max_energies = (350, 350, 1000, 1500)
    locators = (50, 50, 200, 500)
    for index, value in enumerate(rare_earths):
        get_llw_energies_plot(material={'crystal': 'YNi2', 'rare_earth': value},
                              max_energy=max_energies[index],
                              y_major=locators[index],
                              y_minor=locators[index] // 5)


@get_time_of_execution
def get_llw_ratios_plot(material: dict, experimental_value, min_value, max_value, y_major, y_minor):
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
        with open(ratio_file_name, 'r', encoding=con.ENCODING) as file:
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
        check_path(graph_file_name)
        data = Data(x=data['x'], y_set=data['y_set'],
                    legend=data['legend'], errors=None)
        get_plot(data=data,
                 scale=Scale(limits=Limits(x_min=-1, x_max=1, y_min=min_value, y_max=max_value),
                             locators=Locators(x_major=0.5, x_minor=0.1,
                                               y_major=y_major, y_minor=y_minor)),
                 labels=Labels(x=r'$x$', y='Energy Transfer Ratio'),
                 info=GraphInfo(filename=graph_file_name,
                                title=fr'$W={parameters["w"]}$',
                                form=con.FIG_FORMAT, 
                                mode='scatter')
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

    get_all_llw_energies()
