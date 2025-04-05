"""The module contains classes for sample, crystal and RE ion."""


from typing import Union

from pretty_repr import RepresentableObject

from core.re_ion import REIon


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

    def __str__(self):
        output = [
            self.crystal.name,
            f'Rare-earth ion: {self.rare_earth.info.symbol};',
            f'Number of 4f-electrons = '
            f'{self.rare_earth.info.number_of_f_electrons};',
            f'J = {self.rare_earth.info.total_momentum_ground};',
        ]
        return '\n'.join(output)
