"""The module contains some common constants and named_tuples that used in this project."""
from os.path import join, dirname, realpath
from collections import namedtuple
from cycler import cycler

BASE_DIR = dirname(dirname(dirname(realpath(__file__))))

DATAFILES_DIR = join(BASE_DIR, 'datafiles')
GRAPHS_DIR = join(BASE_DIR, 'graphs')

DATA_NAMES = (
    'energies',
    'experiment',
    'ratios',
    'parameters',
    'spectra',
    'susceptibilities',
)

DATA_PATHS = {key: join(DATAFILES_DIR, key) for key in DATA_NAMES}
GRAPHS_PATHS = {key: join(GRAPHS_DIR, key) for key in DATA_NAMES}

RATIOS_NAMES = ('$E_2/E_1$', '$E_3/E_1$', '$E_4/E_1$', '$E_5/E_1$', '$E_6/E_1$',
                '$E_3/E_2$', '$E_4/E_2$', '$E_5/E_2$', '$E_6/E_2$',
                '$E_4/E_3$', '$E_5/E_3$', '$E_6/E_3$',
                '$E_5/E_4$', '$E_6/E_4$',
                '$E_6/E_5$')
FIG_SIZE = (10, 10)
FONT_FAMILY = 'Times New Roman'
DEFAULT_CYCLER = (cycler(color=['black', 'red', 'green', 'blue',
                                'magenta', 'purple', 'Olive', 'DarkCyan',
                                'DarkOrange', 'pink']) +
                  cycler(linestyle=['-', '--', ':', '-.',
                                    '-', '--', ':', '-.',
                                    '-', '--', ]))

FIG_FORMAT = 'png'
ENCODING = 'utf-8'
ENERGY_TRANSFER = 'Energy Transfer, meV'
SCATTERING = r'$S(\omega)$, arb.u.'
INFINITY = float('inf')
PM = chr(177)

Element = namedtuple('Element', ['name',
                                 'total_momentum_ground',
                                 'lande_factor',
                                 'matrix_size',
                                 'f_6',
                                 'radial_integrals',
                                 'stevens_factors'])
Material = namedtuple('Material', ['rare_earth', 'crystal'])
CrossPoint = namedtuple('CrossPoint', ['rare_earth', 'w', 'x', 'ratio_name', 'difference'])
Data = namedtuple('Data', ['x', 'y_set', 'errors', 'legend'])
Text = namedtuple('Text', ['x', 'y', 'string'])
Scale = namedtuple('Scale', ['limits', 'locators'])
