"""
This module is used for checking the project code for compliance with PEP8.

"""

from pylint_af import PyLinter


if __name__ == '__main__':
    PyLinter().check()
