"""The module contains some common functions that used in this project."""
from datetime import datetime
from numpy import transpose, zeros
from .tabular_information import ACCEPTABLE_RARE_EARTHS


def get_sign(value):
    """Returns minus, if argument is negative,
    else it returns plus."""
    return '-' if value < 0 else '+'


def get_value_with_sign(value):
    """Returns float number as a string with sign plus or minus."""
    return f'{get_sign(value)}' + f'{abs(value): .3f}'.lstrip(' ')


def get_new_if_old_is_none(old_value, new_value):
    """Returns new_value if old_value is None, else it returns old_value."""
    return new_value if (old_value is None) else old_value


def value_to_write(value, text_separator):
    """Returns float number as a string
    with tabulation symbol or newline one in the end."""
    return f'{value:10.5f}{text_separator}'


def check_input(choice):
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


def create_table(*arrays):
    """Return transposed array of arrays."""
    data = []
    for array in arrays:
        data.append(array)
    return transpose(data)


def user_input():
    """Returns data inputted by user."""
    crystal = input('Input the name of crystal (e.g. "YNi2"): ')
    rare_earth = check_input('rare')
    w_parameter = check_input('w')
    x_parameter = check_input('x')
    return crystal, rare_earth, w_parameter, x_parameter


def get_empty_matrix(size, dimension=2):
    """Returns 1D or 2D array filled by zeros."""
    empty_matrix = None
    if dimension == 2:
        empty_matrix = zeros((size, size), dtype='float64')
    if dimension == 1:
        empty_matrix = zeros(size, dtype='float64')
    return empty_matrix


def get_time_of_execution(function):
    """Prints time of function's execution."""
    def wrapper(*args, **kwargs):
        start_time = datetime.now()
        function(*args, **kwargs)
        finish_time = datetime.now()
        print(f'Saving time: {finish_time - start_time}')
    return wrapper


def saving_file(function):
    """Prints information about saved file."""
    def wrapper(file, object_to_save, *args, **kwargs):
        print('\n',
              f'***Saving {object_to_save}***:',
              f'Saving file "{file}"...',
              'It may take some time...', sep='\n')
        function(file, object_to_save, *args, **kwargs)
        print(f'File "{file}" is saved.')
    return wrapper


if __name__ == '__main__':
    pass
