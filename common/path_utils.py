"""
The module contains functions working with paths to directories used
in this project.

"""


import os


from common.constants import BASE_DIR, Material, DATA_PATHS, GRAPHS_PATHS
from common.utils import get_value_with_sign


class PathProcessor:

    def __init__(self, path: str) -> None:
        self.path = path

    def create_parent_dirs(self) -> None:
        """Create parent directories for the path, if they do not exist."""
        paths = [self.path]
        while not os.path.exists(paths[-1]):
            paths.append(os.path.dirname(paths[-1]))
        for path in paths[-2:0:-1]:
            os.mkdir(path)

    def remove_if_exists(self):
        """Create parent dirs for the file and remove it, if it exists."""
        self.create_parent_dirs()
        if os.path.exists(self.path):
            os.remove(self.path)


def get_paths(
        data_name: str,
        format_name='.dat',
        is_graph=False,
        material: Material = None,
        parameters: dict = None,
):
    """Returns path of the file that will be saved."""
    os.chdir(BASE_DIR)
    short_name = ''
    if material:
        if isinstance(material.rare_earth, str):
            short_name = f'{material.crystal}_{material.rare_earth}'
        else:
            short_name = f'{material.crystal}_{material.rare_earth.name}'
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
        result_path = os.path.join(
            GRAPHS_PATHS[data_name],
            short_name,
            f'{data_name}_{full_name}',
        )
    else:
        result_path = os.path.join(
            DATA_PATHS[data_name],
            short_name,
            f'{data_name}_{full_name}{format_name}',
        )
    PathProcessor(result_path).create_parent_dirs()
    return result_path
