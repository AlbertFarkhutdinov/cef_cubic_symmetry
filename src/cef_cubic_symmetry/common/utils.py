"""
The module contains some common functions that used in this project.

"""


from datetime import datetime
from json import load

from numpy import zeros

from common.constants import DATA_DIR, INFINITY


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
        result += f'{value:11.5f}\t'
    file.write(f'{result.strip()}\n')


def check_input(choice: str):
    """Checks a value inputted by user, returns it,
    if it satisfies the condition, else requests input again."""
    result = 0
    condition = False
    while not condition:
        if choice == 'rare':
            request = (
                'Input the name of RE ion '
                # f'({", ". join(ACCEPTABLE_RARE_EARTHS)}): '
            )
            result = input(request).capitalize()
            # condition = (result in ACCEPTABLE_RARE_EARTHS)
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


def get_empty_matrix(size: int, dimension=2):
    """Returns 1D or 2D array filled by zeros."""
    sizes = size if dimension == 1 else (size, size)
    return zeros(sizes, dtype='float64')


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
        print(f'Saving time: {datetime.now() - start_time}\n')
    return wrapper


def get_label(number: int, choice=0):
    """Return label for legend"""
    index = 1 if choice != 0 else choice
    return (fr'$E_{number}$', fr'$I_{number}$')[index]


def get_ratios_names(choice=0):
    """Returns list of ratios names"""
    letter = 'E' if choice == 0 else 'I'
    result = []
    for low in range(1, 7):
        for high in range(low + 1, 7):
            result.append(f'${letter}_{high}/{letter}_{low}$')
    return result


def get_repr(obj, *args):
    """Method returns string representation of the object."""
    result = f'{obj.__class__.__name__}('
    for arg in args:
        result += f'{arg}={obj.__getattribute__(arg)!r}, '
    return f'{result.rstrip(", ")})'


class UTF8File:
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


def get_json_object(file_name: str):
    """Returns object from JSON file"""
    with UTF8File(str(DATA_DIR / file_name)) as file:
        obj = load(file)
    return obj
