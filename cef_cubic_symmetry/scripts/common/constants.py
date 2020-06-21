"""The module contains some common constants and named_tuples that used in this project."""
from os.path import join, dirname, realpath
from collections import namedtuple

BASE_DIR = dirname(dirname(dirname(realpath(__file__))))

DATAFILES_DIR = join(BASE_DIR, 'datafiles')
GRAPHS_DIR = join(BASE_DIR, 'graphs')
JSON_DIR = join(BASE_DIR, 'json')

DATA_NAMES = (
    'energies',
    'intensities',
    'experiment',
    'ratios_energies',
    'ratios_intensities',
    'parameters',
    'spectra',
    'susceptibilities',
    'intensities_on_temperature',
)

DATA_PATHS = {key: join(DATAFILES_DIR, key) for key in DATA_NAMES}
GRAPHS_PATHS = {key: join(GRAPHS_DIR, key) for key in DATA_NAMES}

X_PARAMETER = r'$x$'
# ENERGY_TRANSFER = 'Energy Transfer, meV'
ENERGY_TRANSFER = 'Передача энергии, мэВ'
# TRANSITION_INTENSITY = 'Transition Intensity, arb.u.'
TRANSITION_INTENSITY = 'Интенсивность перехода, отн.ед.'
# ENERGY_TRANSFER_RATIO = 'Energy Transfer Ratio'
ENERGY_TRANSFER_RATIO = 'Отношение передач энергии'
# TRANSITION_INTENSITY_RATIO = 'Transition Intensity Ratio'
TRANSITION_INTENSITY_RATIO = 'Отношение интенсивностей перехода'
# SCATTERING = r'$S(\omega)$, arb.u.'
SCATTERING = r'$S(\omega)$, отн.ед.'

SPECTRUM_LABELS = {
    'xlabel': ENERGY_TRANSFER,
    'ylabel': SCATTERING,
}

INFINITY = float('inf')
PM = chr(177)

Element = namedtuple(
    typename='Element',
    field_names=[
        'number_of_f_electrons',
        'name',
        'total_momentum_ground',
        'lande_factor',
        'matrix_size',
        'f_6',
        'radial_integrals',
        'stevens_factors',
    ]
)
Material = namedtuple(
    typename='Material',
    field_names=[
        'rare_earth',
        'crystal',
    ]
)
CrossPoint = namedtuple(
    typename='CrossPoint',
    field_names=[
        'rare_earth',
        'w',
        'x',
        'ratio_name',
        'difference',
    ]
)
Data = namedtuple(
    typename='Data',
    field_names=[
        'x',
        'y_set',
        'errors',
        'legend',
    ]
)
Text = namedtuple(
    typename='Text',
    field_names=[
        'x',
        'y',
        'string',
    ]
)
Scale = namedtuple(
    typename='Scale',
    field_names=[
        'limits',
        'locators',
    ]
)
