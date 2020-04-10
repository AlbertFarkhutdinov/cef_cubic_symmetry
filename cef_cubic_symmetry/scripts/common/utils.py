"""The module contains some common functions that used in this project."""
from datetime import datetime
from numpy import zeros
from scripts.common.tabular_information import ACCEPTABLE_RARE_EARTHS
from scripts.common.constants import INFINITY


def get_sign(value: float):
    """Returns minus, if argument is negative,
    else it returns plus."""
    return '-' if value < 0 else '+'


def get_value_with_sign(value: float):
    """Returns float number as a string with sign plus or minus."""
    if value:
        return f'{get_sign(value)}{abs(value):.3f}'
    return None


def get_default(value, default):
    """Returns default if value is None, else it returns value."""
    return default if (value is None) else value


def write_row(file, row):
    """Writes to file the row of the float numbers
    separated with tabulation symbol."""
    result = ''
    for value in row:
        result += f'{value:10.5f}\t'
    file.write(f'{result.strip()}\n')


def check_input(choice: str):
    """Checks a value inputted by user, returns it,
    if it satisfies the condition, else requests input again."""
    result = 0
    condition = False
    while not condition:
        if choice == 'rare':
            request = f'Input the name of RE ion ({", ". join(ACCEPTABLE_RARE_EARTHS)}): '
            result = input(request).capitalize()
            condition = (result in ACCEPTABLE_RARE_EARTHS)
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


def user_input():
    """Returns data inputted by user."""
    crystal = input('Input the name of crystal (e.g. "YNi2"): ')
    rare_earth = check_input('rare')
    w_parameter = check_input('w')
    x_parameter = check_input('x')
    return crystal, rare_earth, w_parameter, x_parameter


def get_empty_matrix(size: int, dimension=2):
    """Returns 1D or 2D array filled by zeros."""
    empty_matrix = None
    if dimension == 2:
        empty_matrix = zeros((size, size), dtype='float64')
    if dimension == 1:
        empty_matrix = zeros(size, dtype='float64')
    return empty_matrix


def data_popping(data: dict, condition):
    """Pops items from data, that satisfy condition"""
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


def get_time_of_execution(function):
    """Prints time of function's execution."""
    def wrapper(*args, **kwargs):
        start_time = datetime.now()
        function(*args, **kwargs)
        finish_time = datetime.now()
        print(f'Saving time: {finish_time - start_time}\n')
    return wrapper


class OpenedFile:
    """Context manager for file opening"""
    def __init__(self, name: str, mode='r'):
        """Initialization of class"""
        self.name = name
        self.file = None
        self.mode = mode

    def __enter__(self):
        """Method for entrance to context manager"""
        if self.mode != 'r':
            print(f'Saving file "{self.name}"...')
        self.file = open(self.name, mode=self.mode, encoding='utf-8')
        return self.file

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Method for exit from context manager"""
        if self.file:
            self.file.close()
            if self.mode != 'r':
                print(f'File "{self.name}" is saved.')


if __name__ == '__main__':
    pass
