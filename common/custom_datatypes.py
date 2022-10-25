from dataclasses import dataclass


@dataclass
class Sample:
    rare_earth: str
    crystal: str


@dataclass
class CEFParameters:
    b20: float
    b40: float
    b60: float
    b22: float
    b42: float
    b62: float
    b43: float
    b63: float
    b44: float
    b64: float
    b66: float


@dataclass
class LLWParameters:
    w_: float
    x_: float


@dataclass
class MagnetField:
    x_: float
    z_: float


@dataclass
class REIon:
    number_of_f_electrons: int
    name: str
    total_momentum_ground: float
    lande_factor: float
    matrix_size: int
    f_6: float
    radial_integrals: float
    stevens_factors: float
