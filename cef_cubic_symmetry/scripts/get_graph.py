"""The module contains functions for plotting graphs."""
from collections import namedtuple
import matplotlib.pyplot as plt
from scripts.common import constants
from scripts.common.path_utils import check_path, get_paths

plt.rcParams['axes.labelpad'] = 0
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['axes.prop_cycle'] = constants.DEFAULT_CYCLER
plt.rcParams['figure.dpi'] = 500
plt.rcParams['figure.figsize'] = [i / 2.54 for i in constants.FIG_SIZE]
plt.rcParams['figure.subplot.bottom'] = 0.11
plt.rcParams['figure.subplot.hspace'] = 0
plt.rcParams['figure.subplot.left'] = 0.18
plt.rcParams['figure.subplot.right'] = 0.97
plt.rcParams['figure.subplot.top'] = 0.90
plt.rcParams['figure.subplot.wspace'] = 0
plt.rcParams['font.family'] = constants.FONT_FAMILY
plt.rcParams['font.size'] = 14
plt.rcParams['legend.frameon'] = False
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['lines.markersize'] = 1.5
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['mathtext.it'] = constants.FONT_FAMILY

for tick in ('xtick', 'ytick'):
    plt.rcParams[f'{tick}.direction'] = 'in'
    plt.rcParams[f'{tick}.major.pad'] = 3
    plt.rcParams[f'{tick}.major.size'] = 6
    plt.rcParams[f'{tick}.major.width'] = 2
    plt.rcParams[f'{tick}.minor.size'] = 3
    plt.rcParams[f'{tick}.minor.width'] = 2

plt.rcParams['xtick.minor.bottom'] = True
plt.rcParams['xtick.minor.top'] = True
plt.rcParams['xtick.top'] = True
plt.rcParams['ytick.minor.left'] = True
plt.rcParams['ytick.minor.right'] = True
plt.rcParams['ytick.right'] = True


Data = namedtuple('Data', ['x', 'y_set', 'legend'])
Labels = namedtuple('Labels', ['x', 'y'])
Limits = namedtuple('Limits', ['x_min', 'x_max', 'y_min', 'y_max'])
Locators = namedtuple('Locators', ['x_major', 'x_minor', 'y_major', 'y_minor'])
Text = namedtuple('Text', ['x', 'y', 'string'])
GraphInfo = namedtuple('GraphInfo', ['filename', 'title', 'mode', 'form'])
Scale = namedtuple('Scale', ['limits', 'locators'])


def get_plot(info=GraphInfo(filename='graph', title=None, mode='line', form='show'),
             labels=Labels(x=r'$x$', y=r'$y$'), scale=None, data=None,
             text=Text(x=0, y=0, string='Test')):
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
    _ax.set_ylim(ymin=-float(scale.locators.y_minor), ymax=scale.limits.y_max)
    for key, y_data in data.y_set.items():
        if info.mode == 'line':
            _ax.plot(data.x, y_data, label=data.legend[key])
        if info.mode == 'scatter':
            _ax.scatter(data.x, y_data, label=data.legend[key])
    _ax.legend()
    plt.text(x=text.x, y=text.y, s=text.string)
    if info.form == 'show':
        plt.show()
    elif info.form == 'eps':
        fig.savefig(f'{info.filename}.eps')
    elif info.form == 'png':
        fig.savefig(f'{info.filename}.png')
    else:
        print(f'Format {info.form} is not supported.')
    plt.close('all')


def get_energy_transfer_plot(material: dict, parameters: dict, temperatures: tuple,
                             data=None, scale=None):
    """Returns graph for inelastic neutron scattering spectrum with specified parameters"""
    form = 'png'
    parameters['T'] = f'{temperatures[0]}-{temperatures[1]}'
    file_name = get_paths(constants.PATH_TO_GRAPHS, 'graph', form,
                          material=material, parameters=parameters)
    check_path(file_name)
    text = Text(x=2.2, y=2300,
                string=fr'$W={parameters["w"]: .3f}$' + '\n' + fr'$x = {parameters["x"]: .3f}$')
    get_plot(data=data, scale=scale, text=text,
             labels=Labels(x='Energy Transfer, meV', y=r'$S(\omega)$, arb.u.'),
             info=GraphInfo(filename=file_name.rstrip(f'.{form}'),
                            title=None, form=form, mode='line')
             )


def get_llw_energies_plot(material: dict, max_energy, y_major, y_minor):
    """Returns graphs for dependence of transfer energies on CEF parameters"""
    form = 'png'
    levels_number = 7
    for w_parameter in (1, -1):
        data = {'x': [], 'y_set': {}, 'legend': {}}
        for level in range(1, levels_number):
            data['y_set'][fr'$E_{level}$'] = []
            data['legend'][fr'$E_{level}$'] = fr'$E_{level}$'

        parameters = {'w': w_parameter}
        energy_file_name = get_paths(constants.PATH_TO_ENERGY_DATAFILES, 'energy', 'dat',
                                     material=material, parameters=parameters)
        with open(energy_file_name, 'r', encoding='utf-8') as file:
            lines = file.readlines()
            for line in lines:
                line = line.rstrip('\n').split('\t')
                line = [float(value.lstrip(' ')) for value in line]
                data['x'].append(line[0])
                for level in range(1, levels_number):
                    key = fr'$E_{level}$'
                    try:
                        data['y_set'][key].append(line[level])
                    except IndexError:
                        data['y_set'][key].append(float('inf'))

        for key, array in data['y_set'].items():
            if all([value == float('inf') for value in array]):
                data['legend'][key] = None
        graph_file_name = get_paths(constants.PATH_TO_GRAPHS, 'energies', form,
                                    material=material, parameters=parameters)
        check_path(graph_file_name)
        text = Text(x=2.2, y=2300, string=fr'$W={parameters["w"]: .3f}$')
        _data = Data(x=data['x'], y_set=data['y_set'],
                     legend=data['legend'])
        get_plot(data=_data, text=text,
                 scale=Scale(limits=Limits(x_min=-1, x_max=1, y_min=-10, y_max=max_energy),
                             locators=Locators(x_major=0.5, x_minor=0.1,
                                               y_major=y_major, y_minor=y_minor)),
                 labels=Labels(x=r'$x$', y='Energy Transfer, meV'),
                 info=GraphInfo(filename=graph_file_name.rstrip(f'.{form}'),
                                title=None, form=form, mode='scatter')
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
    get_llw_energies_plot(material={'crystal': 'YNi2', 'rare_earth': 'Tb'},
                          max_energy=350, y_major=50, y_minor=10)
    get_llw_energies_plot(material={'crystal': 'YNi2', 'rare_earth': 'Tm'},
                          max_energy=350, y_major=50, y_minor=10)
    get_llw_energies_plot(material={'crystal': 'YNi2', 'rare_earth': 'Er'},
                          max_energy=1000, y_major=200, y_minor=40)
    get_llw_energies_plot(material={'crystal': 'YNi2', 'rare_earth': 'Ho'},
                          max_energy=1500, y_major=500, y_minor=100)