"""The module contains classes for sample, crystal and RE ion."""


from typing import Union

import pandas as pd
from pretty_repr import RepresentableObject

from common.constants import DATA_DIR


class REIon(RepresentableObject):
    """
    Class with rare earth (RE) ion information.

    Parameters
    ----------
    identifier : str or int
        Identifier of RE ion.
        May be element symbol (Ce, Pr, etc.) or number of 4f-electrons.

    Attributes
    ----------
    info : pd.Series
        Information about RE ion.

        Attribute includes fields:
         - number_of_f_electrons - number of 4f-electrons;
         - element - RE element full name (cerium, praseodymium, etc.);
         - symbol - RE element symbol (Ce, Pr, etc.);
         - total_momentum_ground - total momentum in the ground state;
         - lande_factor - Lande g-factor;
         # TODO describe fields
         - f_6 - ;
         - radial_integral_2 - ;
         - radial_integral_4 - ;
         - radial_integral_6 - ;
         - stevens_factor_2 - ;
         - stevens_factor_4 - ;
         - stevens_factor_6 - ;

    """

    def __init__(
            self,
            identifier: Union[str, int],
    ) -> None:
        """Initialize self. See help(type(self)) for accurate signature."""
        self.identifier = identifier
        data = pd.read_csv(DATA_DIR / 'rare_earths_properties.csv')
        self.info = data.loc[
            data.symbol == self.identifier.capitalize()
            if isinstance(self.identifier, str)
            else data.number_of_f_electrons == self.identifier
        ].reset_index(drop=True).T[0].rename(None)

    @property
    def excluded_attributes_for_repr(self) -> set[str]:
        return {'info'}

    @property
    def matrix_size(self):
        return 2 * self.info.total_momentum_ground + 1

    @property
    def squared_momentum(self):
        __total_momentum_ground = self.info.total_momentum_ground
        return __total_momentum_ground * (__total_momentum_ground + 1)


class Crystal(RepresentableObject):
    """
    Class with crystal information.

    Parameters
    ----------
    name : str
        The name of crystal.

    """

    def __init__(self, name: str = '?'):
        """Initialize self. See help(type(self)) for accurate signature."""
        self.name = name


class Sample(RepresentableObject):
    """
    Class with sample information.

    Parameters
    ----------
    crystal : str or Crystal
        Information about crystal.
    rare_earth : str or RareEarth
        Information about impurity RE ion.

    """

    def __init__(
            self,
            crystal: Union[str, Crystal],
            rare_earth: Union[str, REIon],
    ):
        """Initialize self. See help(type(self)) for accurate signature."""
        self.crystal = (
            Crystal(name=crystal) if isinstance(crystal, str) else crystal
        )
        self.rare_earth = (
            REIon(identifier=rare_earth)
            if isinstance(rare_earth, str)
            else rare_earth
        )
