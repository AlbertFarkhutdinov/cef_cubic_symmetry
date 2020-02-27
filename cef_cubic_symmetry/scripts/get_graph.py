"""The module contains functions for plotting graphs."""
from collections import namedtuple
import matplotlib.pyplot as plt
from .common.constants import PATH_TO_GRAPHS, DEFAULT_CYCLER, FIG_SIZE, FONT_FAMILY
from .common.path_utils import check_path, get_paths

plt.rcParams['axes.labelpad'] = 0
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['axes.prop_cycle'] = DEFAULT_CYCLER
plt.rcParams['figure.dpi'] = 500
plt.rcParams['figure.figsize'] = [i / 2.54 for i in FIG_SIZE]
plt.rcParams['figure.subplot.bottom'] = 0.11
plt.rcParams['figure.subplot.hspace'] = 0
plt.rcParams['figure.subplot.left'] = 0.19
plt.rcParams['figure.subplot.right'] = 0.98
plt.rcParams['figure.subplot.top'] = 0.90
plt.rcParams['figure.subplot.wspace'] = 0
plt.rcParams['font.family'] = FONT_FAMILY
plt.rcParams['font.size'] = 14
plt.rcParams['legend.frameon'] = False
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['mathtext.it'] = FONT_FAMILY

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
GraphInfo = namedtuple('GraphInfo', ['filename', 'title', 'mode'])
Scale = namedtuple('Scale', ['limits', 'locators'])


def get_plot(info=GraphInfo(filename='graph', title=None, mode='show'),
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
    dpi = 100 if info.mode == 'show' else 500
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
        _ax.plot(data.x, y_data, label=data.legend[key])
    _ax.legend()
    plt.text(x=text.x, y=text.y, s=text.string)
    if info.mode == 'show':
        plt.show()
    elif info.mode == 'eps':
        fig.savefig(f'{info.filename}.eps')
    elif info.mode == 'png':
        fig.savefig(f'{info.filename}.png')
    else:
        print(f'Mode {info.mode} is not supported.')
    plt.close('all')


def get_energy_transfer_plot(material: dict, parameters: dict, temperatures: tuple,
                             data=None, scale=None):
    """Returns graph for inelastic neutron scattering spectrum with specified parameters"""
    mode = 'png'
    parameters['T'] = f'{temperatures[0]}-{temperatures[1]}'
    file_name = get_paths(PATH_TO_GRAPHS, 'graph', mode,
                          material=material, parameters=parameters)
    check_path(file_name)
    text = Text(x=2.2, y=2300,
                string=fr'$W={parameters["w"]: .3f}$' + '\n' + fr'$x = {parameters["x"]: .3f}$')
    get_plot(data=data, scale=scale, text=text,
             labels=Labels(x='Energy Transfer, meV', y=r'$S(\omega)$, arb.u.'),
             info=GraphInfo(filename=file_name.rstrip(f'.{mode}'),
                            title=None, mode=mode)
             )


if __name__ == '__main__':
    MAXIMUM = 11
    X_ARRAY = list(range(MAXIMUM))
    Y_ARRAY = {
        '4': [i ** 4 for i in range(MAXIMUM)],
        '3': [i ** 3 for i in range(MAXIMUM)],
        '2': [i ** 2 for i in range(MAXIMUM)],
        '1': [i ** 1 for i in range(MAXIMUM)],
    }
    COLORS = {
        '4': 'black',
        '3': 'red',
        '2': 'green',
        '1': 'blue',
    }
    LEGEND = {
        '4': '$x^4$',
        '3': '$x^3$',
        '2': '$x^2$',
        '1': '$x^1$',
    }
    get_plot(data=Data(x=X_ARRAY, y_set=Y_ARRAY, legend=LEGEND))
