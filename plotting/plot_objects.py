"""The module contains class for graphs plotting."""


from typing import Optional

import matplotlib.pyplot as plt
from cycler import cycler

from common import constants as con, utils as ut


class CustomPlot:
    """Description of Plot object"""

    @staticmethod
    def _set_plot_parameters():
        """Setting of rcParams"""
        plt.rcParams.update(ut.get_json_object('plot_parameters.json'))
        prop_cycle = 'axes.prop_cycle'
        plt.rcParams[prop_cycle] = (
                cycler(color=plt.rcParams[prop_cycle]['color'])
                + cycler(linestyle=plt.rcParams[prop_cycle]['linestyle'])
        )
        plt.rcParams['figure.figsize'] = [i / 2.54 for i in (10, 10)]
        tick_parameters = {
            'direction': 'in',
            'major.pad': 3,
            'major.size': 6,
            'major.width': 2,
            'minor.size': 3,
            'minor.width': 2,
        }
        for tick in ('xtick', 'ytick'):
            for _key, _value in tick_parameters.items():
                plt.rcParams[f'{tick}.{_key}'] = _value

    def __init__(self, dpi: int = 300):
        self.dpi = dpi
        self.limits = {
            _key: None
            for _key
            in ('x_min', 'x_max', 'y_min', 'y_max')
        }
        self.fig = None
        self._ax = None

    def __enter__(self):
        """Method for entrance to context manager"""
        self._set_plot_parameters()
        self.fig, self._ax = plt.subplots(dpi=self.dpi)
        return self

    def set_labels(
            self,
            x_label: Optional[str] = None,
            y_label: Optional[str] = None,
            title: Optional[str] = None,
    ):
        """Sets labels of axis and plot"""
        self._ax.set_xlabel = x_label
        self._ax.set_ylabel = y_label
        self._ax.set_title = title

    def set_limits(
            self,
            x_min: Optional[float] = None,
            x_max: Optional[float] = None,
            y_min: Optional[float] = None,
            y_max: Optional[float] = None,
    ):
        """Sets limits of x and y intervals"""
        y_set = self.data.y_set.values()
        _limits = {
            'x_min': (x_min, min(self.data.x)),
            'x_max': (x_max, max(self.data.x)),
            'y_min': (y_min, min((min(value) for value in y_set))),
            'y_max': (y_max, max((max(value) for value in y_set))),
        }
        for _key, _value in _limits.items():
            self.limits[_key] = _value[0] or _value[1]
        if self._ax:
            for axis in ('x', 'y'):
                axis_limits = {
                    f'{axis}{lim}': self.limits[f'{axis}_{lim}']
                    for lim in ('min', 'max')
                }
                self._ax.__getattribute__(f'set_{axis}lim')(**axis_limits)

    def set_locators(
            self,
            x_major=None,
            x_minor=None,
            y_major=None,
            y_minor=None,
    ):
        """Sets major and minor ticks for plot"""
        majors = (
            x_major or (self.limits['x_max'] - self.limits['x_min']) // 5,
            y_major or (self.limits['y_max'] - self.limits['y_min']) // 5,
        )
        minors = (x_minor or majors[0] / 5, y_minor or majors[1] / 5)
        if self._ax:
            for i, axis in enumerate((self._ax.xaxis, self._ax.yaxis)):
                axis.set_major_locator(plt.MultipleLocator(majors[i]))
                axis.set_minor_locator(plt.MultipleLocator(minors[i]))

    def make_plot(
            self,
            mode='plot',
            text: con.Text = None,
    ):
        """Draws the plot at specified mode"""
        if self.fig and self._ax:
            functions = {
                key: self._ax.__getattribute__(key)
                for key in ('plot', 'scatter', 'errorbar')
            }
            for key, y_data in self.data.y_set.items():
                args = (self.data.x, y_data)
                kwargs = {
                    'label': self.data.legend[key]
                }
                if mode == 'errorbar':
                    kwargs['yerr'] = self.data.errors[key]
                    kwargs['fmt'] = 'o'
                    kwargs['elinewidth'] = 1
                functions[mode](*args, **kwargs)
            self._ax.legend()
            if text:
                plt.text(x=text.x, y=text.y, s=text.string)

    def save_or_show(
            self,
            filename=None,
            form=None,
    ):
        """Saves or shows the plot"""
        if self.fig and self._ax:
            if form:
                self.fig.savefig(f'{filename}.{form}')
                print(f'Graph {filename}.{form} is saved.')
            else:
                plt.show()

    def save_in_two_forms(
            self,
            filename: str,
            form_1='png',
            form_2='eps',
    ):
        """Saves or shows the plot"""
        self.save_or_show(filename=filename, form=form_1)
        self.save_or_show(filename=filename, form=form_2)

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Method for exit from context manager"""
        if self.fig and self._ax:
            plt.close('all')

    def __repr__(self):
        """Method returns string representation of the Plot object."""
        return ut.get_repr(self, 'data', 'dpi')
