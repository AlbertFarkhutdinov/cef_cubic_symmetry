"""The module contains some common constants that used in this project."""
from os.path import join, dirname, realpath
from cycler import cycler

BASE_DIR = dirname(dirname(dirname(realpath(__file__))))

DATAFILES_DIR = join(BASE_DIR, 'datafiles')
GRAPHS_DIR = join(BASE_DIR, 'graphs')

PATH_TO_ENERGY_DATAFILES = join(DATAFILES_DIR, 'energies')
PATH_TO_EXPERIMENTAL_DATAFILES = join(DATAFILES_DIR, 'experiment')
PATH_TO_RATIO_DATAFILES = join(DATAFILES_DIR, 'ratios')
PATH_TO_SAVED_OBJECTS = join(DATAFILES_DIR, 'parameters')
PATH_TO_SPECTRA_DATAFILES = join(DATAFILES_DIR, 'spectra')
PATH_TO_SUSCEPTIBILITY_DATAFILES = join(DATAFILES_DIR, 'susceptibilities')

PATH_TO_ENERGY_GRAPHS = join(GRAPHS_DIR, 'energies')
PATH_TO_RATIO_GRAPHS = join(GRAPHS_DIR, 'ratios')
PATH_TO_SPECTRA_GRAPHS = join(GRAPHS_DIR, 'spectra')

RATIOS_NAMES = ('$E_2/E_1$', '$E_3/E_1$', '$E_4/E_1$', '$E_5/E_1$', '$E_6/E_1$',
                '$E_3/E_2$', '$E_4/E_2$', '$E_5/E_2$', '$E_6/E_2$',
                '$E_4/E_3$', '$E_5/E_3$', '$E_6/E_3$',
                '$E_5/E_4$', '$E_6/E_4$',
                '$E_6/E_5$')
FIG_SIZE = (10, 10)
TWO_FIG_SIZE = (20, 10)
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
