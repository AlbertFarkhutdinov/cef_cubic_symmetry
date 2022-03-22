from common.tabular_information import Element
from common.tabular_information import RARE_EARTHS, RARE_EARTHS_NAMES


class Material:

    def __init__(
            self,
            rare_earth_symbol: Element,
            crystal: str = '?',
    ):
        self.rare_earth = RARE_EARTHS[
            RARE_EARTHS_NAMES.index(rare_earth_symbol)
        ]
        self.crystal = crystal
