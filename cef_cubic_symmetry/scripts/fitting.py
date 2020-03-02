"""This module contains fitting Lorentz function to experimental data"""
from collections import namedtuple
from os.path import join
from numpy.random import random
from numpy import array, diag, sqrt
from scipy.optimize import curve_fit
import scripts.get_graph as gg
from scripts.common.constants import PATH_TO_EXPERIMENTAL_DATAFILES, ENCODING
from scripts.common.physics import lorentz

Peak = namedtuple('Peak', ['center', 'amplitude', 'width'])


def lorentz_1(arg, center, wid, amp):
    """Lorentz function for one peak"""
    return amp * lorentz(arg, center, wid)


def lorentz_2(arg, center1, wid1, amp1, center2, wid2, amp2):
    """Lorentz function for two peaks"""
    return (lorentz_1(arg, center1, wid1, amp1) +
            lorentz_1(arg, center2, wid2, amp2))


def lorentz_3(arg, center1, wid1, amp1, center2, wid2, amp2, center3, wid3, amp3):
    """Lorentz function for three peaks"""
    return (lorentz_1(arg, center1, wid1, amp1) +
            lorentz_1(arg, center2, wid2, amp2) +
            lorentz_1(arg, center3, wid3, amp3))


def get_data_from_file(file_name):
    """Returns three arrays (x, y, error) from file"""
    experimental_data = [[], [], []]
    with open(file_name, 'r', encoding=ENCODING) as file:
        lines = file.readlines()
        for line in lines:
            line = line.rstrip('\n').split('\t')
            for column, value in enumerate(line):
                experimental_data[column].append(float(value))

    for i, column in enumerate(experimental_data):
        experimental_data[i] = array(column)
    return tuple(experimental_data)


def fitting(function, x_array, y_array, parameters):
    """function for fitting"""
    p_opt, p_cov = curve_fit(function, x_array, y_array, p0=parameters)
    p_err = sqrt(diag(p_cov))
    return p_opt, p_err


def fitting_test():
    """test fitting possibility"""
    x_array = array([0.01*i for i in range(2000)])
    y_array = lorentz_3(x_array, 8, 0.6, 15, 10, 0.5, 18, 12, 9, 0.7)
    noise_array = random(len(x_array)) / 10
    p_opt, _ = fitting(lorentz_3, x_array, y_array,
                       [8, 0.6, 15, 10, 0.5, 18, 12, 9, 0.7])
    gg.get_plot(data=gg.Data(x=x_array,
                             y_set={'lorentz': y_array + noise_array,
                                    'fit': lorentz_3(x_array, *p_opt),
                                    },
                             legend={'lorentz': 'lorentz',
                                     'fit': 'fit',
                                     },
                             )
                )


if __name__ == '__main__':
    # fitting_test()
    X_DATA, Y_DATA, ERROR = get_data_from_file(join(PATH_TO_EXPERIMENTAL_DATAFILES,
                                                    'PSI_Tb_YNi2_3meV_1.6K.txt'))
    OPT, ERR = fitting(lorentz_1, X_DATA, Y_DATA, [0, 1248, 0.05])
    print(OPT, ERR, sep='\n')
    gg.get_plot(data=gg.Data(x=X_DATA,
                             y_set={'lorentz': Y_DATA,
                                    'fit': lorentz_1(X_DATA, *OPT)},
                             legend={'lorentz': 'lorentz', 'fit': 'fit'}))
