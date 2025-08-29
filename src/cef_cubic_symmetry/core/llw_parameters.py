from dataclasses import dataclass

from cef_cubic_symmetry.core.b_parameters import BParameters

F4 = 60


@dataclass
class LLWParameters:
    """CEF parameters in LLW-notation."""

    w_: float
    x_: float
    f6: float

    def get_b_parameters(self, f6: float) -> BParameters:
        """Return CEF parameters in B-notation."""
        f4 = F4
        b_parameters = BParameters(
            b40=self.w_ * self.x_ / f4,
            b60=self.w_ * (1 - abs(self.x_)) / f6,
        )
        b_parameters.b44 = 5 * b_parameters.b40
        b_parameters.b64 = -21 * b_parameters.b60
        return b_parameters
