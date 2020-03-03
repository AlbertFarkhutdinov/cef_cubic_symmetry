"""This module contains unittests for module utils.py"""
import unittest
from sys import path as sys_path
from os import getcwd
from os.path import join
sys_path.append(join(getcwd(), '..'))
from scripts.common import path_utils
from scripts.common import constants


class PathUtilsTests(unittest.TestCase):
    """Class with tests for path_utils.py"""
    def test_get_paths_without_parameters(self):
        """Tests for function get_paths without parameters"""
        self.assertEqual(join(constants.DATAFILES_DIR,
                              'energies\\YNi2_Tm\\energy_YNi2_Tm.dat'),
                         path_utils.get_paths(constants.PATH_TO_ENERGY_DATAFILES,
                                              'energy',
                                              material={'crystal': 'YNi2', 'rare_earth': 'Tm'}))

    def test_get_paths_with_w(self):
        """Tests for function get_paths with parameter W"""
        self.assertEqual(join(constants.DATAFILES_DIR,
                              'energies\\YNi2_Tm\\energy_YNi2_Tm_w+1.000.dat'),
                         path_utils.get_paths(constants.PATH_TO_ENERGY_DATAFILES,
                                              'energy',
                                              material={'crystal': 'YNi2', 'rare_earth': 'Tm'},
                                              parameters={'w': 1}))

    def test_get_paths_with_x(self):
        """Tests for function get_paths with parameter x"""
        self.assertEqual(join(constants.DATAFILES_DIR,
                              'energies\\YNi2_Tm\\energy_YNi2_Tm_x-1.000.dat'),
                         path_utils.get_paths(constants.PATH_TO_ENERGY_DATAFILES,
                                              'energy',
                                              material={'crystal': 'YNi2', 'rare_earth': 'Tm'},
                                              parameters={'x': -1}))


if __name__ == '__main__':
    unittest.main()
