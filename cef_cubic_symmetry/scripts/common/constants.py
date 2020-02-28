"""The module contains some common constants that used in this project."""
from os.path import join, dirname, realpath
from cycler import cycler

BASE_DIR = dirname(dirname(dirname(realpath(__file__))))
DATAFILES_DIR = join(BASE_DIR, 'datafiles')
PATH_TO_ENERGY_DATAFILES = join(DATAFILES_DIR, 'energies')
PATH_TO_EXPERIMENTAL_DATAFILES = join(DATAFILES_DIR, 'experiment')
PATH_TO_RATIO_DATAFILES = join(DATAFILES_DIR, 'ratios')
PATH_TO_SAVED_OBJECTS = join(DATAFILES_DIR, 'parameters')
PATH_TO_SPECTRA_DATAFILES = join(DATAFILES_DIR, 'spectra')
PATH_TO_SUSCEPTIBILITY_DATAFILES = join(DATAFILES_DIR, 'susceptibilities')
PATH_TO_GRAPHS = join(DATAFILES_DIR, 'graphs')
RATIOS_NAMES = ('E2/E1', 'E3/E1', 'E4/E1', 'E5/E1',
                'E3/E2', 'E4/E2', 'E5/E2',
                'E4/E3', 'E5/E3',
                'E5/E4')
FIG_SIZE = (10, 10)
TWO_FIG_SIZE = (20, 10)
FONT_FAMILY = 'Times New Roman'
DEFAULT_CYCLER = (cycler(color=['black', 'red', 'green', 'blue', 'magenta', 'DarkOrange']) +
                  cycler(linestyle=['-', '--', ':', '-.', '-', '--']))
