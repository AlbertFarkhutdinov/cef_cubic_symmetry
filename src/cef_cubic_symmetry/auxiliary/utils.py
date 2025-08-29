"""The module contains some common functions that used in this project."""

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from numpy import zeros

from cef_cubic_symmetry.auxiliary.constants import DATA_DIR, INFINITY


def get_sign(value: float) -> str:
    """Return minus, if argument is negative, otherwise return plus."""
    return '-' if value < 0 else '+'


def get_value_with_sign(value: float) -> str | None:
    """Return float number as a string with sign plus or minus."""
    if value:
        return f'{get_sign(value)}{abs(value):.3f}'
    return None


def get_default(value, default):
    """Return default if value is None, else it returns value."""
    return default if (value is None) else value


def write_row(file, row) -> None:
    """Write the row of the float numbers to the file."""
    result = ''.join(
        f'{value:11.5f}\t'
        for value in row
    )
    file.write(f'{result.strip()}\n')


def check_input(choice: str) -> float:
    """
    Check input.

    The method inspects a value inputted by user, returns it,
    if it satisfies the condition, else requests input again.

    """
    result = 0
    condition = False
    while not condition:
        if choice == 'rare':
            request = 'Input the name of RE ion'
            result = input(request).capitalize()
            condition = True
        try:
            if choice == 'w':
                result = float(input('Input |W| > 0: '))
                condition = (result > 0)
            if choice == 'x':
                result = float(input('Input -1 <= x <= 1: '))
                condition = (-1 <= result <= 1)
            if choice == 'intervals':
                result = int(input('Input the number of intervals > 0: '))
                condition = (result > 0)
        except ValueError:
            condition = False
    return result


def get_empty_matrix(size: int, dimension=2) -> np.ndarray:
    """Return 1D or 2D array filled by zeros."""
    sizes = size if dimension == 1 else (size, size)
    return zeros(sizes, dtype='float64')


def data_popping(data: dict, condition) -> None:
    """Pop items from data, that satisfy condition."""
    popped_number = 0
    for key, array in data['y_set'].copy().items():
        finite_array = [value for value in array if value != INFINITY]
        if condition(finite_array):
            data['y_set'].move_to_end(key)
            data['legend'].move_to_end(key)
            popped_number += 1
    for _ in range(popped_number):
        data['y_set'].popitem()
        data['legend'].popitem()


def get_time_of_execution(function) -> callable:
    """Print time of function's execution."""

    def wrapper(*args, **kwargs) -> None:
        start_time = datetime.now(tz=timezone.utc)
        function(*args, **kwargs)
        finish_time = datetime.now(tz=timezone.utc) - start_time
        print(f'Saving time: {finish_time}\n')

    return wrapper


def get_label(number: int, choice=0) -> str:
    """Return label for legend."""
    index = int(choice != 0)
    return (fr'$E_{number}$', fr'$I_{number}$')[index]


def get_ratios_names(choice=0) -> list:
    """Return list of ratios names."""
    letter = 'E' if choice == 0 else 'I'
    result = []
    for low in range(1, 7):
        result.extend(
            f'${letter}_{high}/{letter}_{low}$'
            for high in range(low + 1, 7)
        )
    return result


def get_repr(obj, *args) -> str:
    """Return string representation of the object."""
    repr_parts = [f'{obj.__class__.__name__}(']
    repr_parts.extend(
        f'{arg}={getattr(obj, arg)!r}, '
        for arg in args
    )
    result = ''.join(repr_parts)
    return f'{result.rstrip(", ")})'


class UTF8File:
    """Context manager for file opening."""

    def __init__(self, name: str, mode='r') -> None:
        """Initialize class."""
        self.name = name
        self.file = None
        self.mode = mode

    def __enter__(self):
        """Execute entrance to context manager and return self."""
        if self.mode != 'r':
            print(f'Saving file "{self.name}"...')
        self.file = Path(self.name).open(mode=self.mode, encoding='utf-8')
        return self.file

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """Execute exit from context manager."""
        if self.file:
            self.file.close()
            if self.mode != 'r':
                print(f'File "{self.name}" is saved.')


def get_json_object(file_name: str):
    """Return object from JSON file."""
    with UTF8File(str(DATA_DIR / file_name)) as file:
        return json.load(file)
