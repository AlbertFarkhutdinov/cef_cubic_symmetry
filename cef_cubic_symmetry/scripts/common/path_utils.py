"""The module contains functions working with paths to directories used in this project."""
from os import chdir, mkdir, remove
from os.path import join, exists, dirname
from scripts.common.constants import BASE_DIR, Material, DATA_PATHS, GRAPHS_PATHS
from scripts.common.utils import get_value_with_sign


def check_parent_dirs(path_to_check: str):
    """Makes parent directories for argument, if they do not exist."""
    condition = False
    paths = []
    while not condition:
        if exists(path_to_check):
            condition = True
        else:
            paths.append(path_to_check)
            path_to_check = dirname(path_to_check)
    for path in paths[:0:-1]:
        mkdir(path)


def get_paths(data_name: str, format_name='.dat', is_graph=False,
              material: Material = None, parameters: dict = None):
    """Returns path of the file that will be saved."""
    chdir(BASE_DIR)
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
        result_path = join(
            GRAPHS_PATHS[data_name],
            short_name,
            f'{data_name}_{full_name}',
        )
    else:
        result_path = join(
            DATA_PATHS[data_name],
            short_name,
            f'{data_name}_{full_name}{format_name}',
        )
    check_parent_dirs(result_path)
    return result_path


def remove_if_exists(path_to_check: str):
    """Makes remove file, if it exists."""
    check_parent_dirs(path_to_check)
    if exists(path_to_check):
        remove(path_to_check)
