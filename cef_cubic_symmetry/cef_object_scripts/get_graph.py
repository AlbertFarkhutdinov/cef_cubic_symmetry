"""The module contains functions for plotting graphs."""
from collections import namedtuple
import matplotlib.pyplot as plt
from cef_object_scripts import common

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 14
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['xtick.top'] = True
plt.rcParams['ytick.right'] = True
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['xtick.minor.width'] = 2
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['ytick.minor.width'] = 2
plt.rcParams['xtick.major.size'] = 6
plt.rcParams['xtick.minor.size'] = 3
plt.rcParams['ytick.major.size'] = 6
plt.rcParams['ytick.minor.size'] = 3
plt.rcParams['figure.subplot.left'] = 0.19
plt.rcParams['figure.subplot.right'] = 0.98
plt.rcParams['figure.subplot.top'] = 0.90
plt.rcParams['figure.subplot.bottom'] = 0.11
plt.rcParams['figure.subplot.wspace'] = 0
plt.rcParams['figure.subplot.hspace'] = 0
FIG_SIZE = (10, 10)
plt.rcParams['figure.figsize'] = [i / 2.54 for i in FIG_SIZE]

print(plt.rcParams)


Data = namedtuple('Data', ['x', 'y_set'])
Labels = namedtuple('Labels', ['x', 'y'])
Limits = namedtuple('Limits', ['x_min', 'x_max', 'y_min', 'y_max'])
Locators = namedtuple('Locators', ['x_major', 'x_minor', 'y_major', 'y_minor'])
Text = namedtuple('Text', ['x', 'y', 'string'])
Legend = namedtuple('Legend', ['label_set', 'color_set'])
GraphInfo = namedtuple('GraphInfo', ['filename', 'title', 'mode'])


def get_plot(info=GraphInfo(filename='graph', title=None, mode='show'),
             labels=Labels(x='x_label', y='y_label'), label_pad=0,
             limits=None, locators=None, data=None,
             legend=Legend(label_set=None, color_set=None),
             text=Text(x=0, y=0, string='Test')):
    """Returns graph with specified parameters"""
    if limits is None and data is not None:
        limits = Limits(x_min=min(data.x),
                        x_max=max(data.x),
                        y_min=min([min(value) for value in data.y_set.values()]),
                        y_max=max([max(value) for value in data.y_set.values()])
                        )
    dpi = 100 if info.mode == 'show' else 500
    fig, _ax = plt.subplots(dpi=dpi)
    _ax.set_title(info.title)
    _ax.set_xlabel(labels.x, labelpad=label_pad)
    _ax.set_ylabel(labels.y, labelpad=label_pad)
    if locators is None and limits is not None:
        locators = Locators(x_major=(limits.x_max - limits.x_min) // 5,
                            x_minor=(limits.x_max - limits.x_min) // 5 / 5,
                            y_major=(limits.y_max - limits.y_min) // 5,
                            y_minor=(limits.y_max - limits.y_min) // 5 / 5)
    _ax.xaxis.set_major_locator(plt.MultipleLocator(locators.x_major))
    _ax.xaxis.set_minor_locator(plt.MultipleLocator(locators.x_minor))
    _ax.yaxis.set_major_locator(plt.MultipleLocator(locators.y_major))
    _ax.yaxis.set_minor_locator(plt.MultipleLocator(locators.y_minor))
    _ax.set_xlim(xmin=limits.x_min, xmax=limits.x_max)
    _ax.set_ylim(ymin=-float(locators.y_minor), ymax=limits.y_max)
    for key, y_data in data.y_set.items():
        _ax.plot(data.x, y_data, color=legend.color_set[key],
                 label=legend.label_set[key], linewidth=2)
    _ax.legend(frameon=False)
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


def get_energy_transfer_plot(material: dict, parameters: dict, temperature_1, temperature_2, mode,
                             data=None, limits=None, locators=None,
                             legend=Legend(label_set=None, color_set=None),
                             text=Text(x=0, y=0, string='')):
    """Returns graph for inelastic neutron scattering spectrum with specified parameters"""
    parameters['T'] = f'{temperature_1}-{temperature_2}'
    file_name = common.get_paths(common.PATH_TO_GRAPHS, 'graph', mode,
                                 material=material, parameters=parameters)
    common.check_path(file_name)
    text.string = fr'$W={parameters["w"]: .3f}$' + '\n' + fr'$x = {parameters["x"]: .3f}$'
    get_plot(data=data, legend=legend,
             labels=Labels(x='Energy Transfer, meV', y=r'$S(\omega)$, arb.u.'),
             info=GraphInfo(filename=file_name.rstrip(f'.{mode}'),
                            title=None, mode=mode),
             limits=limits, text=text, locators=locators)


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
    get_plot(data=Data(x=X_ARRAY, y_set=Y_ARRAY),
             legend=Legend(color_set=COLORS, label_set=LEGEND))
