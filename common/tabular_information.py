"""
The module contains some tabular data that used in this project.

"""


from fractions import Fraction

from common.constants import Element


F4 = 60
BOHR_MAGNETON = 5.788382e-2  # meV/T

CERIUM = Element(
    number_of_f_electrons=1,
    name='Ce',
    total_momentum_ground=2.5,
    lande_factor=Fraction(6, 7),
    matrix_size=6,
    f_6=0,
    radial_integrals={
        '2': 0.3666,
        '4': 0.3108,
        '6': 0.5119,
    },
    stevens_factors={
        '2': Fraction(-2, 35),
        '4': Fraction(2, 315),
        '6': 0,
    },
)

PRASEODYMIUM = Element(
    number_of_f_electrons=2,
    name='Pr',
    total_momentum_ground=4,
    lande_factor=Fraction(4, 5),
    matrix_size=9,
    f_6=1260,
    radial_integrals={
        '2': 0.3380,
        '4': 0.2670,
        '6': 0.4150,
    },
    stevens_factors={
        '2': Fraction(-52, 2475),
        '4': Fraction(-4, 5445),
        '6': Fraction(272, 4459455),
    }
)

NEODYMIUM = Element(
    number_of_f_electrons=3,
    name='Nd',
    total_momentum_ground=4.5,
    lande_factor=Fraction(8, 11),
    matrix_size=10,
    f_6=2520,
    radial_integrals={
        '2': 0.3120,
        '4': 0.2015,
        '6': 0.3300,
    },
    stevens_factors={
        '2': Fraction(-7, 1089),
        '4': Fraction(-136, 467181),
        '6': Fraction(-1615, 42513471),
    }
)

PROMETHIUM = Element(
    number_of_f_electrons=4,
    name='Pm',
    total_momentum_ground=4,
    lande_factor=Fraction(3, 5),
    matrix_size=9,
    f_6=1260,
    radial_integrals={
        '2': 0.2917,
        '4': 0.1488,
        '6': 0.2787,
    },
    stevens_factors={
        '2': Fraction(14, 1815),
        '4': Fraction(952, 2335905),
        '6': Fraction(2584, 3864861),
    }
)

SAMARIUM = Element(
    number_of_f_electrons=5,
    name='Sm',
    total_momentum_ground=2.5,
    lande_factor=Fraction(2, 7),
    matrix_size=6,
    f_6=0,
    radial_integrals={
        '2': 0.2728,
        '4': 0.1772,
        '6': 0.2317,
    },
    stevens_factors={
        '2': Fraction(13, 315),
        '4': Fraction(26, 10395),
        '6': 0,
    }
)

EUROPIUM = Element(
    number_of_f_electrons=6,
    name='Eu',
    total_momentum_ground=0,
    lande_factor=0,
    matrix_size=1,
    f_6=0,
    radial_integrals={
        '2': 0.2569,
        '4': 0.1584,
        '6': 0.1985,
    },
    stevens_factors={
        '2': 0,
        '4': 0,
        '6': 0,
    }
)

GADOLINIUM = Element(
    number_of_f_electrons=7,
    name='Gd',
    total_momentum_ground=3.5,
    lande_factor=2,
    matrix_size=8,
    f_6=1260,
    radial_integrals={
        '2': 0.2428,
        '4': 0.1427,
        '6': 0.1720,
    },
    stevens_factors={
        '2': 0,
        '4': 0,
        '6': 0,
    }
)

TERBIUM = Element(
    number_of_f_electrons=8,
    name='Tb',
    total_momentum_ground=6,
    lande_factor=Fraction(3, 2),
    matrix_size=13,
    f_6=7560,
    radial_integrals={
        '2': 0.2302,
        '4': 0.1295,
        '6': 0.1505,
    },
    stevens_factors={
        '2': Fraction(-1, 99),
        '4': Fraction(2, 16335),
        '6': Fraction(-1, 891891),
    }
)

DYSPROSIUM = Element(
    number_of_f_electrons=9,
    name='Dy',
    total_momentum_ground=7.5,
    lande_factor=Fraction(4, 3),
    matrix_size=16,
    f_6=13860,
    radial_integrals={
        '2': 0.2188,
        '4': 0.1180,
        '6': 0.1328,
    },
    stevens_factors={
        '2': Fraction(-2, 315),
        '4': Fraction(-8, 135135),
        '6': Fraction(4, 3864861),
    }
)

HOLMIUM = Element(
    number_of_f_electrons=10,
    name='Ho',
    total_momentum_ground=8,
    lande_factor=Fraction(5, 4),
    matrix_size=17,
    f_6=13860,
    radial_integrals={
        '2': 0.2085,
        '4': 0.1081,
        '6': 0.1810,
    },
    stevens_factors={
        '2': Fraction(-1, 450),
        '4': Fraction(-1, 30030),
        '6': Fraction(-5, 3864861),
    }
)

ERBIUM = Element(
    number_of_f_electrons=11,
    name='Er',
    total_momentum_ground=7.5,
    lande_factor=Fraction(6, 5),
    matrix_size=16,
    f_6=13860,
    radial_integrals={
        '2': 0.1991,
        '4': 0.0996,
        '6': 0.1058,
    },
    stevens_factors={
        '2': Fraction(4, 1575),
        '4': Fraction(2, 45045),
        '6': Fraction(8, 3864861),
    }
)

THULIUM = Element(
    number_of_f_electrons=12,
    name='Tm',
    total_momentum_ground=6,
    lande_factor=Fraction(7, 6),
    matrix_size=13,
    f_6=7560,
    radial_integrals={
        '2': 0.1905,
        '4': 0.0921,
        '6': 0.0953,
    },
    stevens_factors={
        '2': Fraction(1, 99),
        '4': Fraction(8, 49005),
        '6': Fraction(-5, 891891),
    }
)

YTTERBIUM = Element(
    number_of_f_electrons=13,
    name='Yb',
    total_momentum_ground=3.5,
    lande_factor=Fraction(8, 7),
    matrix_size=8,
    f_6=1260,
    radial_integrals={
        '2': 0.1826,
        '4': 0.0854,
        '6': 0.0863,
    },
    stevens_factors={
        '2': Fraction(2, 63),
        '4': Fraction(-2, 1155),
        '6': Fraction(4, 27027),
    }
)

RARE_EARTHS = (
    CERIUM, PRASEODYMIUM, NEODYMIUM, PROMETHIUM, SAMARIUM, EUROPIUM,
    GADOLINIUM, TERBIUM, DYSPROSIUM, HOLMIUM, ERBIUM, THULIUM, YTTERBIUM
)

ACCEPTABLE_RARE_EARTHS = tuple(
    element.name
    for element in RARE_EARTHS
    if element.f_6 != 0
)
RARE_EARTHS_NAMES = [element.name for element in RARE_EARTHS]
