"""
The module contains some common constants and named_tuples that used
in this project.

"""


from collections import namedtuple
from pathlib import Path

BASE_DIR = Path('.').resolve().parent

DATA_DIR = BASE_DIR / 'data'
PLOTS_DIR = BASE_DIR / 'plots'

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

DATA_PATHS = {key: DATA_DIR / key for key in DATA_NAMES}
PLOT_PATHS = {key: PLOTS_DIR / key for key in DATA_NAMES}

X_PARAMETER = r'$x$'
ENERGY_TRANSFER = 'Energy Transfer, meV'
TRANSITION_INTENSITY = 'Transition Intensity, arb.u.'
ENERGY_TRANSFER_RATIO = 'Energy Transfer Ratio'
TRANSITION_INTENSITY_RATIO = 'Transition Intensity Ratio'
SCATTERING = r'$S(\omega)$, arb.u.'

SPECTRUM_LABELS = {
    'xlabel': ENERGY_TRANSFER,
    'ylabel': SCATTERING,
}

RESOLUTION = 1e-2
THRESHOLD = 1e-4
INFINITY = float('inf')
PM = chr(177)

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
