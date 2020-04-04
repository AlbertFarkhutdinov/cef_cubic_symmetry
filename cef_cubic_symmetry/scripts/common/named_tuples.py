"""The module contains named tuples that used in this project."""

from collections import namedtuple


Element = namedtuple('Element', ['name',
                                 'total_momentum_ground',
                                 'lande_factor',
                                 'matrix_size',
                                 'f_6',
                                 'radial_integrals',
                                 'stevens_factors'])
Material = namedtuple('Material', ['rare_earth', 'crystal'])
CrossPoint = namedtuple('CrossPoint', ['rare_earth', 'w', 'x', 'ratio_name', 'difference'])
Data = namedtuple('Data', ['x', 'y_set', 'errors', 'legend'])
Text = namedtuple('Text', ['x', 'y', 'string'])
Scale = namedtuple('Scale', ['limits', 'locators'])
