import os
from typing import Optional

import pandas as pd

from common.constants import DATAFILES_DIR


class RareEarth:

    def __init__(
            self,
            symbol: Optional[str],
            number_of_f_electrons: Optional[int],
    ):
        data = pd.read_csv(
            os.path.join(DATAFILES_DIR, 'rare_earths_properties.csv')
        )
        if symbol:
            self.info = data.loc[data.symbol == symbol, :].reset_index().T[0]
            self.info.name = None
        # if number_of_f_electrons:
        #     self.info = data.loc[data.number_of_f_electrons == number_of_f_electrons, :].reset_index().T[0]
        #     self.info.name = None
