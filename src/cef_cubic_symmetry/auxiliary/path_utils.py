"""The module for work with paths to directories used in this project."""


import os
from pathlib import Path

from cef_cubic_symmetry.auxiliary import paths
from cef_cubic_symmetry.auxiliary.utils import get_value_with_sign
from cef_cubic_symmetry.core import Sample


class PathProcessor:

    def __init__(self, path: Path) -> None:
        """Initialize self. See help(type(self)) for accurate signature."""
        self.path = path

    def create_parent_dirs(self) -> None:
        """Create parent directories for the path, if they do not exist."""
        paths = [self.path]
        while not paths[-1].exists():
            paths.append(paths[-1].parent)
        for path in paths[-2:0:-1]:
            path.mkdir()

    def remove_if_exists(self) -> None:
        """Create parent dirs for the file and remove it, if it exists."""
        self.create_parent_dirs()
        if self.path.exists():
            self.path.unlink()


def get_paths(
    data_name: str,
    format_name: str = '.dat',
    sample: Sample = None,
    parameters: dict | None = None,
    *,
    is_graph: bool = False,
) -> Path:
    """Return path of the file that will be saved."""
    os.chdir(paths.ROOT_PATH)
    short_name = ''
    if sample:
        short_name = f'{sample.crystal.name}_{sample.rare_earth.info.symbol}'

    full_name_parts = [short_name]
    if parameters:
        for key, value in parameters.items():
            if key in {'w', 'x'} and value is not None:
                full_name_parts.append(f'_{key}{get_value_with_sign(value)}')
            elif key == 'T':
                full_name_parts.append(f'_{key}{value}')
            elif key == 'setup':
                full_name_parts.append(f'_{key}_{value}')
            else:
                full_name_parts.append(f'_{key}{value:.3f}')
    full_name = ''.join(full_name_parts)
    if is_graph:
        result_path = paths.PLOT_PATHS[data_name].joinpath(
            short_name,
            f'{data_name}_{full_name}',
        )
    else:
        result_path = paths.DATA_PATHS[data_name].joinpath(
            short_name,
            f'{data_name}_{full_name}{format_name}',
        )
    PathProcessor(result_path).create_parent_dirs()
    return result_path
