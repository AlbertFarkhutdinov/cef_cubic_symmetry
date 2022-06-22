"""This module contains fitting Lorentz function to experimental data."""


import os
import sys

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

from scripts.plot_objects import CustomPlot
from common.constants import PM, INFINITY, DATA_PATHS, Data
from common.physics import gaussian, multi_lorentzian, multi_gaussian


def get_data_from_file(file_name: str) -> pd.DataFrame:
    """Returns three arrays (x, y, error) from file"""
    return pd.read_csv(file_name, sep='\t', names=['x', 'y', 'errors'])


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


def multi_lorentzian_with_gauss(
        arg: float,
        center: float,
        width: float,
        amplitude: float,
        *parameters,
):
    """Returns value of multi_peak function for lorentzian."""
    return (gaussian(arg, center, width, amplitude) +
            multi_lorentzian(arg, *parameters))


def simple_fitting(data: dict):
    """Simple fitting"""
    data = filtered_data(data)
    start_width = 0.1
    peak_0 = (0, start_width, 130)
    peak_1 = (0.2, start_width, 1.7)
    peak_2 = (1.5, start_width, 0.2)
    p_opt, p_err = fitting(
        multi_lorentzian_with_gauss,
        data,
        parameters=(*peak_0, 0, *peak_1, *peak_2),
        min_value=-2,
        max_value=2,
    )
    print_peak_parameters(
        multi_lorentzian_with_gauss,
        p_opt,
        p_err,
    )
    fitted = multi_lorentzian_with_gauss(data['x'], *p_opt)
    return fitted


if __name__ == '__main__':
    EXPERIMENTAL_DATA = get_data_from_file(
        os.path.join(
            DATA_PATHS['experiment'],
            'PSI_Tb_YNi2_3meV_1.6K.dat',
        )
    )
    EXPERIMENTAL_DATA = filtered_data(
        EXPERIMENTAL_DATA,
        min_value=-2,
        max_value=2,
    )
    FITTED_DATA = simple_fitting(EXPERIMENTAL_DATA)
    with CustomPlot(
            data=Data(
                x=EXPERIMENTAL_DATA['x'],
                y_set={
                    'exp': EXPERIMENTAL_DATA['y'],
                    'fit': FITTED_DATA,
                },
                legend={
                    'exp': 'exp',
                    'fit': 'fit',
                },
                errors=None,
            )
    ) as plot:
        plot.set_labels(
            xlabel='x_test',
            ylabel='y_test',
            title='test',
        )
        # plot.set_locators()
        plot.make_plot()


if __name__ == '__main__':
    DATA = {
        'x': np.array([0.01 * i for i in range(-300, 500)]),
    }
    PEAK_1 = (0, 0.1, 100)
    PEAK_2 = (0.5, 0.15, 20)
    PEAK_3 = (4, 0.15, 20)
    DATA['y'] = multi_lorentzian(DATA['x'], 0.5, *PEAK_1, *PEAK_2, *PEAK_3)
    DATA['y'] += np.random.rand(len(DATA['x']))
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
