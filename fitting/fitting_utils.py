"""
This module contains fitting Lorentz function to experimental data.

"""


import sys

import numpy as np
from numpy.random import rand
from scipy.optimize import curve_fit

from common.constants import PM, INFINITY, Data
from common.physics import multi_lorentzian, multi_gaussian
from common.utils import OpenedFile
from scripts.plot_objects import CustomPlot


def get_data_from_file(
        file_name: str,
):
    """Returns three arrays (x, y, error) from file"""
    experimental_data = {
        'x': [],
        'y': [],
        'errors': [],
    }
    with OpenedFile(file_name) as file:
        for line in file:
            row = [float(value) for value in line.rstrip('\n').split('\t')]
            for index, key in enumerate(experimental_data.keys()):
                experimental_data[key].append(row[index])
    return experimental_data


def filtered_data(
        data: dict,
        min_value=-INFINITY,
        max_value=INFINITY,
):
    """Returns filtered data"""
    rows = [
        (i, j, k)
        for (i, j, k) in zip(*data.values())
        if min_value <= i <= max_value
    ]
    result = {}
    for index, key in enumerate(data.keys()):
        result[key] = np.array([value[index] for value in rows])
    return result


def print_peak_parameters(
        function,
        values,
        errors=None,
):
    """Prints parameters of peak"""
    result = f'Function name: {function.__name__}\n'
    for index, value in enumerate(values):
        result += f'{index}: {value:.3f}'
        if errors.any() and len(errors) == len(values):
            result += f' {PM} {errors[index]:.3f}'
        result += ';\n'
    print(result)


def fitting(
        function,
        data: dict,
        parameters,
        min_value: float,
        max_value: float,
):
    """Returns parameters of function fitted to data with one peak"""
    data = filtered_data(data, min_value, max_value)
    p_opt, p_cov = curve_fit(
        f=function,
        xdata=data['x'],
        ydata=data['y'],
        p0=parameters
    )
    p_err = np.sqrt(np.diag(p_cov))
    return p_opt, p_err


def multi_peak_fitting(
        function_name: str,
        data: dict,
        parameters,
        min_value: float,
        max_value: float,
):
    """Returns parameters of function fitted to data
    with several peaks (lorentzian or gaussian)"""
    function = (multi_lorentzian
                if function_name.lower() == 'lorentz'
                else multi_gaussian)
    data = filtered_data(data)
    if (len(parameters) - 1) % 3:
        print('The parameters number does not match peaks number!')
        sys.exit(1)
    else:
        p_opt, p_err = fitting(
            function=function,
            data=data,
            parameters=parameters,
            min_value=min_value,
            max_value=max_value,
        )
        return p_opt, p_err


if __name__ == '__main__':
    DATA = {
        'x': np.array([0.01 * i for i in range(-300, 500)]),
    }
    PEAK_1 = (0, 0.1, 100)
    PEAK_2 = (0.5, 0.15, 20)
    PEAK_3 = (4, 0.15, 20)
    DATA['y'] = multi_lorentzian(DATA['x'], 0.5, *PEAK_1, *PEAK_2, *PEAK_3)
    DATA['y'] += rand(len(DATA['x']))
    DATA['errors'] = DATA['y'] * 0.01
    START_PARAMETERS = (0.2, *PEAK_1, *PEAK_2, *PEAK_3)
    P_OPT, P_ERR = multi_peak_fitting(
        function_name='lorentz',
        data=DATA,
        parameters=START_PARAMETERS,
        min_value=-2,
        max_value=3,
    )
    print_peak_parameters(multi_lorentzian, P_OPT, P_ERR)
    with CustomPlot(
            data=Data(
                x=DATA['x'],
                y_set={
                    'exp': DATA['y'],
                    'fit': multi_lorentzian(DATA['x'], *P_OPT),
                },
                legend={
                    'exp': 'experiment',
                    'fit': 'fit',
                },
                errors={
                    'exp': DATA['errors'],
                    'fit': None,
                },
            )
    ) as test_plot:
        test_plot.set_labels(xlabel='x_test', ylabel='y_test', title='test')
        # test_plot.set_limits(x_min=0, x_max=5, y_min=0, y_max=100)
        # test_plot.set_locators()
        test_plot.make_plot(mode='plot')
