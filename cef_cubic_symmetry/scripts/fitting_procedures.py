"""This module contains fitting Lorentz function to experimental data"""
from os.path import join
from scripts.plot_objects import CustomPlot
from scripts.common.constants import DATA_PATHS, Data
from scripts.common.physics import gaussian, multi_lorentzian
import scripts.common.fitting_utils as fu


def multi_lorentzian_with_gauss(arg: float, center: float, width: float,
                                amplitude: float, *parameters):
    """Returns value of multi_peak function for lorentzian"""
    return (gaussian(arg, center, width, amplitude) +
            multi_lorentzian(arg, *parameters))


def simple_fitting(data: dict):
    """Simple fitting"""
    data = fu.filtered_data(data)
    start_width = 0.1
    peak_0 = (0, start_width, 130)
    peak_1 = (0.2, start_width, 1.7)
    peak_2 = (1.5, start_width, 0.2)
    p_opt, p_err = fu.fitting(
        multi_lorentzian_with_gauss,
        data,
        parameters=(*peak_0, 0, *peak_1, *peak_2),
        min_value=-2,
        max_value=2,
    )
    fu.print_peak_parameters(
        multi_lorentzian_with_gauss,
        p_opt,
        p_err,
    )
    fitted = multi_lorentzian_with_gauss(
        data['x'],
        *p_opt,
    )
    return fitted


if __name__ == '__main__':
    EXPERIMENTAL_DATA = fu.get_data_from_file(
        join(
            DATA_PATHS['experiment'],
            'PSI_Tb_YNi2_3meV_1.6K.dat',
        )
    )
    EXPERIMENTAL_DATA = fu.filtered_data(
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
        # plot.set_limits(
        #     x_min=2,
        #     x_max=8,
        #     y_min=0,
        #     y_max=5000,
        # )
        plot.set_labels(
            x_label='x_test',
            y_label='y_test',
            title='test',
        )
        # plot.set_locators()
        plot.make_plot()
