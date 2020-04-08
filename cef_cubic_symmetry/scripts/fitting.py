"""This module contains fitting Lorentz function to experimental data"""
from os.path import join
# from numpy.random import random
import numpy as np
from scipy.optimize import curve_fit
from scripts.plot_objects import CustomPlot
from scripts.common.constants import DATA_PATHS, Data, PM
from scripts.common.physics import lorentz, gauss
from scripts.common.utils import OpenedFile


def get_data_from_file(file_name):
    """Returns three arrays (x, y, error) from file"""
    experimental_data = {
        'x': [],
        'y': [],
        'errors': [],
    }
    with OpenedFile(file_name) as file:
        for line in file:
            row = line.rstrip('\n').split('\t')
            row = [float(value) for value in row]
            experimental_data['x'].append(row[0])
            experimental_data['y'].append(row[1])
            experimental_data['errors'].append(row[2])
    return experimental_data


def filtered_data(data, min_value=-float('inf'), max_value=float('inf')):
    """Returns filtered data"""
    data = [(i, j, k) for (i, j, k) in zip(data['x'], data['y'], data['errors'])
            if min_value <= i <= max_value]
    return {
        'x': np.array([value[0] for value in data]),
        'y': np.array([value[1] for value in data]),
        'errors': np.array([value[2] for value in data]),
    }


def lorentzian(arg, center, width, amplitude):
    """Lorentz function for one peak"""
    return amplitude * lorentz(arg, center, width)


def gaussian(arg, center, width, amplitude):
    """Gauss function for one peak"""
    return amplitude * gauss(arg, center, width)


def one_peak_fitting(function, data, parameters, min_value, max_value):
    """function for one peak fitting"""
    data = filtered_data(data, min_value, max_value)
    p_opt, p_cov = curve_fit(function, data['x'], data['y'], p0=parameters)
    p_err = np.sqrt(np.diag(p_cov))
    print(
        f'Type: {function.__name__}\n',
        f'Center ={p_opt[0]: .3f} {PM}{p_err[0]: .3f};',
        f'Width ={p_opt[1]: .3f} {PM}{p_err[1]: .3f};',
        f'Amplitude ={p_opt[2]: .3f} {PM}{p_err[2]: .3f};',
        sep=' '
    )
    parameters = zip(p_opt, p_err)
    return parameters


def multi_peak(function, arg, parameters):
    """Function for many peaks"""
    background = parameters[0]
    peaks_parameters = parameters[1:]
    peaks = []
    for index in range(0, len(peaks_parameters), 3):
        peaks.append(function(arg, *peaks_parameters[index: index + 3]))
    return background + sum(peaks)


def multi_lorentzian(arg, parameters):
    """Lorentz function for many peaks"""
    return multi_peak(lorentzian, arg, parameters)


def multi_gaussian(arg, parameters):
    """Lorentz function for many peaks"""
    return multi_peak(gaussian, arg, parameters)


def get_difference(parameters, data, p_opt):
    """Returns difference between multi_lorentzian and data"""
    return data['y'] - multi_lorentzian(data['x'], parameters) - gaussian(data['x'], *p_opt)


def simple_fitting(data):
    """Simple fitting"""
    data = filtered_data(data)
    start_width = 0.1
    parameters_0 = one_peak_fitting(gaussian,
                                    data,
                                    parameters=[0, start_width, 130],
                                    min_value=-2,
                                    max_value=0.1)
    peak_0 = gaussian(data['x'], *[i for (i, j) in parameters_0])
    data['y'] -= peak_0
    result = peak_0
    for start_values in ([0.2, start_width, 1.7],
                         [1.5, start_width, 0.2],
                         ):
        optimized = one_peak_fitting(lorentzian,
                                     data,
                                     parameters=start_values,
                                     min_value=0,
                                     max_value=2)
        peak = lorentzian(data['x'], *[i for (i, j) in optimized])
        data['y'] -= peak
        result += peak
    return result, data['y']


# def complicated_fitting(data):
#     """Complicated fitting"""
#     data = filtered_data(data)
#     start_width = 0.1
#     parameters_0 = one_peak_fitting(gaussian,
#                                     data,
#                                     parameters=[0, start_width, 130],
#                                     min_value=-2,
#                                     max_value=0.1)
#     p_opt = [i for (i, j) in parameters_0]
#     data = filtered_data(data, min_value=0.1)
#     x_array = data['x']
#     y_array = data['y']
#     start_parameters = [min(y_array)]
#     fitted_data = get_difference(start_parameters, data, p_opt)
#     counter = 0
#     optimized_parameters = start_parameters[:]
#
#     while max(fitted_data) - min(fitted_data) > 0.1:
#         counter += 1
#         print(f'\nCounter = {counter}')
#         maximum = np.argmax(fitted_data)
#         y_max = fitted_data[maximum]
#         center = x_array[maximum]
#         start_parameters += [center, start_width, (y_max - min(y_array))]
#         try:
#             optimized_parameters, _ = leastsq(
#                 func=get_difference,
#                 x0=np.array(start_parameters),
#                 args=(x_array, y_array, p_opt)
#             )
#         except TypeError:
#             print('Type error!')
#             return x_array, y_array, multi_lorentzian(x_array, optimized_parameters)
#         print(f'Optimized parameters:\nBackground = {optimized_parameters[0]: .3f}')
#         for i in range(counter):
#             center = optimized_parameters[1 + 3 * i]
#             width = optimized_parameters[2 + 3 * i]
#             amplitude = optimized_parameters[3 + 3 * i]
#             print(
#                 f'Center_{i + 1} = {center: .3f}',
#                 f'Width_{i + 1} = {width: .3f}',
#                 f'Amplitude_{i + 1} = {amplitude: .3f}',
#                 sep='\t'
#             )
#         fitted_data = get_difference(optimized_parameters, data, p_opt)
#     print('Done!')
#     return x_array, y_array, multi_lorentzian(x_array, optimized_parameters)
#
#
# def fitting_test():
#     """test fitting possibility"""
#     data = {
#         'x': np.array([0.01 * i for i in range(2000)]),
#     }
#     data['y'] = multi_lorentzian(data['x'], (0, 8, 0.6, 15, 10, 0.5, 18, 12, 9, 0.7))
#     noise_array = random(len(data['x'])) / 10
#     data['y'] += noise_array
#     x_array, y_array, fit = complicated_fitting(data)
#     with CustomPlot(data=Data(x=x_array,
#                               y_set={
#                                   'lorentz': y_array,
#                                   'fit': fit,
#                               },
#                               legend={
#                                   'lorentz': 'lorentz',
#                                   'fit': 'fit',
#                               },
#                               errors=None,
#                               )
#                     ) as test_plot:
#         # test_plot.set_limits(x_min=2, x_max=8, y_min=0, y_max=5000)
#         test_plot.set_labels(x_label='x_test', y_label='y_test', title='test')
#         # test_plot.set_locators()
#         test_plot.make_plot()


if __name__ == '__main__':
    # fitting_test()
    EXPERIMENTAL_DATA = get_data_from_file(join(DATA_PATHS['experiment'],
                                                'PSI_Tb_YNi2_3meV_1.6K.txt'))
    EXPERIMENTAL_DATA = filtered_data(EXPERIMENTAL_DATA,
                                      min_value=-2,
                                      max_value=2)
    FITTED_DATA, REST_DATA = simple_fitting(EXPERIMENTAL_DATA)
    with CustomPlot(data=Data(x=EXPERIMENTAL_DATA['x'],
                              y_set={
                                  'exp': EXPERIMENTAL_DATA['y'],
                                  'fit': FITTED_DATA,
                                  'rest': REST_DATA,
                              },
                              legend={
                                  'exp': 'exp',
                                  'fit': 'fit',
                                  'rest': 'rest',
                              },
                              errors=None,
                              )
                    ) as plot:
        # plot.set_limits(x_min=2, x_max=8, y_min=0, y_max=5000)
        plot.set_labels(x_label='x_test', y_label='y_test', title='test')
        # plot.set_locators()
        plot.make_plot()
