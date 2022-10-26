"""
The module contains functions working with paths to directories used
in this project.

"""


import os
from pathlib import Path


from common.constants import BASE_DIR, DATA_PATHS, PLOT_PATHS
from common.utils import get_value_with_sign
from core import Sample


class PathProcessor:

    def __init__(self, path: Path) -> None:
        self.path = path

    def create_parent_dirs(self) -> None:
        """Create parent directories for the path, if they do not exist."""
        paths = [self.path]
        while not paths[-1].exists():
            paths.append(paths[-1].parent)
        for path in paths[-2:0:-1]:
            path.mkdir()

    def remove_if_exists(self):
        """Create parent dirs for the file and remove it, if it exists."""
        self.create_parent_dirs()
        if self.path.exists():
            self.path.unlink()


def get_paths(
        data_name: str,
        format_name='.dat',
        is_graph=False,
        sample: Sample = None,
        parameters: dict = None,
):
    """Returns path of the file that will be saved."""
    os.chdir(BASE_DIR)
    short_name = ''
    if sample:
        short_name = f'{sample.crystal.name}_{sample.rare_earth.info.symbol}'

    full_name = short_name
    if parameters:
        for key, value in parameters.items():
            if key in ('w', 'x') and value is not None:
                full_name += f'_{key}{get_value_with_sign(value)}'
            elif key == 'T':
                full_name += f'_{key}{value}'
            elif key == 'setup':
                full_name += f'_{key}_{value}'
            else:
                full_name += f'_{key}{value:.3f}'
    if is_graph:
        result_path = (
            PLOT_PATHS[data_name] / short_name / f'{data_name}_{full_name}'
        )
    else:
        result_path = (
            DATA_PATHS[data_name] / short_name /
            f'{data_name}_{full_name}{format_name}'
        )
    PathProcessor(result_path).create_parent_dirs()
    return result_path
