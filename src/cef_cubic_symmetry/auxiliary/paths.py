"""The module contains paths used in this project."""


from pathlib import Path

ROOT_PATH = Path(__file__).parents[3]

DATA_PATH = ROOT_PATH.joinpath('data')
PLOTS_PATH = ROOT_PATH.joinpath('plots')

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

DATA_PATHS = {key: DATA_PATH / key for key in DATA_NAMES}
PLOT_PATHS = {key: PLOTS_PATH / key for key in DATA_NAMES}


if __name__ == '__main__':
    print(ROOT_PATH)
    print(DATA_PATH)
    print(PLOTS_PATH)
