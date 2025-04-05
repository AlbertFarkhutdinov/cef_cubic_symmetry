from dataclasses import dataclass


@dataclass
class BParameters:
    """CEF parameters in B-notation."""
    b20: float = 0
    b40: float = 0
    b60: float = 0
    b22: float = 0
    b42: float = 0
    b62: float = 0
    b43: float = 0
    b63: float = 0
    b44: float = 0
    b64: float = 0
    b66: float = 0


@dataclass
class LLWParameters:
    """CEF parameters in LLW-notation."""
    w_: float
    x_: float
    f_6: float

    def get_b_parameters(self, f_6: float) -> BParameters:
        """Return CEF parameters in B-notation."""
        f_4 = 60
        b_parameters = BParameters(
            b40=self.w_ * self.x_ / f_4,
            b60=self.w_ * (1 - abs(self.x_)) / f_6
        )
        b_parameters.b44 = 5 * b_parameters.b40
        b_parameters.b64 = -21 * b_parameters.b60
        return b_parameters
