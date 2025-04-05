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
