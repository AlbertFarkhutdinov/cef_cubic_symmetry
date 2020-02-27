"""This module contains unittests for module utils.py"""
import unittest
from os.path import join
from scripts.common import path_utils
from scripts.common import constants


class PathUtilsTests(unittest.TestCase):
    """Class with tests for path_utils.py"""
    def test_get_paths(self):
        """Tests for function get_paths"""
        material = {'crystal': 'YNi2', 'rare_earth': 'Tm'}
        self.assertEqual(join(constants.DATAFILES_DIR,
                              'energies\\YNi2_Tm\\energy_YNi2_Tm.dat'),
                         path_utils.get_paths(constants.PATH_TO_ENERGY_DATAFILES,
                                              'energy', 'dat', material))
        self.assertEqual(join(constants.DATAFILES_DIR,
                              'energies\\YNi2_Tm\\energy_YNi2_Tm_w+1.000.dat'),
                         path_utils.get_paths(constants.PATH_TO_ENERGY_DATAFILES,
                                              'energy', 'dat', material, {'w': 1}))
        self.assertEqual(join(constants.DATAFILES_DIR,
                              'energies\\YNi2_Tm\\energy_YNi2_Tm_w-1.000.dat'),
                         path_utils.get_paths(constants.PATH_TO_ENERGY_DATAFILES,
                                              'energy', 'dat', material, {'w': -1}))
        self.assertEqual(join(constants.DATAFILES_DIR,
                              'energies\\YNi2_Tm\\energy_YNi2_Tm_x+1.000.dat'),
                         path_utils.get_paths(constants.PATH_TO_ENERGY_DATAFILES,
                                              'energy', 'dat', material, {'x': 1}))
        self.assertEqual(join(constants.DATAFILES_DIR,
                              'energies\\YNi2_Tm\\energy_YNi2_Tm_x-1.000.dat'),
                         path_utils.get_paths(constants.PATH_TO_ENERGY_DATAFILES,
                                              'energy', 'dat', material, {'x': -1}))


if __name__ == '__main__':
    unittest.main()
