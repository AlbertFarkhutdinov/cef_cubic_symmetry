from core import Sample


class CalculationRunner:

    def __init__(
            self,
            rare_earth_ion: str,
            crystal: str,
    ) -> None:
        self.sample = Sample(crystal=crystal, rare_earth=rare_earth_ion)
