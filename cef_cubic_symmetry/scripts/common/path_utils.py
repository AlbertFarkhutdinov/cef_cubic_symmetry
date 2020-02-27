"""The module contains functions working with paths to directories used in this project."""
from os import chdir, mkdir, remove
from os.path import join, exists, dirname
from .constants import BASE_DIR
from .utils import get_value_with_sign


def get_paths(directory, data_name, format_name, material: dict = None, parameters: dict = None):
    """Returns path of the file that will be saved."""
    chdir(BASE_DIR)
    short_name = ''
    if material:
        short_name = f'{material["crystal"]}_{material["rare_earth"]}'
    full_name = short_name
    if parameters:
        for key, value in parameters.items():
            if key != 'T':
                full_name += f'_{key}{get_value_with_sign(value)}'
            else:
                full_name += f'_{key}{value}'
    return join(directory, short_name, f'{data_name}_{full_name}.{format_name}')


def check_path(path_to_check):
    """Makes parent directories for argument, if they do not exist."""
    condition = False
    levels = []
    while not condition:
        if exists(path_to_check):
            condition = True
        else:
            levels.insert(0, path_to_check)
            path_to_check = dirname(path_to_check)
    for level in range(len(levels) - 1):
        mkdir(levels[level])


def remove_if_exists(path_to_check):
    """Makes remove file, if it exists."""
    check_path(path_to_check)
    if exists(path_to_check):
        remove(path_to_check)
