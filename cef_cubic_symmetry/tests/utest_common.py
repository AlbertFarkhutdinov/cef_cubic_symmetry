"""This module contains unittests for module common.py"""
import unittest
import os
from cef_object_scripts import common


class CommonTests(unittest.TestCase):
    """Class with tests for common.py"""
    def test_get_sign(self):
        """Tests for function get_sign"""
        self.assertEqual('+', common.get_sign(1))
        self.assertEqual('+', common.get_sign(0))
        self.assertEqual('-', common.get_sign(-1))

    def test_value_to_write(self):
        """Tests for function value_to_write"""
        self.assertEqual('  3.14159\t', common.value_to_write(3.1415926535, '\t'))
        self.assertEqual(' -3.14159\n', common.value_to_write(-3.1415926535, '\n'))

    def test_get_temperature(self):
        """Tests for function get_temperature"""
        self.assertEqual(0, common.get_temperature(0, 10))
        self.assertEqual(10, common.get_temperature(10, 0))
        self.assertEqual(0, common.get_temperature(None, 0))
        self.assertEqual(0, common.get_temperature(0, None))
        self.assertEqual(None, common.get_temperature(None, None))

    def test_get_paths(self):
        """Tests for function get_paths"""
        self.assertEqual(os.path.join(common.BASE_DIR,
                                      'energy_datafiles\\YNi2_Tm\\energy_YNi2_Tm.dat'),
                         common.get_paths(common.PATH_TO_ENERGY_DATAFILES,
                                          'energy', 'dat', 'YNi2', 'Tm'))
        self.assertEqual(os.path.join(common.BASE_DIR,
                                      'energy_datafiles\\YNi2_Tm\\energy_YNi2_Tm_w+1.000.dat'),
                         common.get_paths(common.PATH_TO_ENERGY_DATAFILES,
                                          'energy', 'dat', 'YNi2', 'Tm', w_parameter=1))
        self.assertEqual(os.path.join(common.BASE_DIR,
                                      'energy_datafiles\\YNi2_Tm\\energy_YNi2_Tm_w-1.000.dat'),
                         common.get_paths(common.PATH_TO_ENERGY_DATAFILES,
                                          'energy', 'dat', 'YNi2', 'Tm', w_parameter=-1))
        self.assertEqual(os.path.join(common.BASE_DIR,
                                      'energy_datafiles\\YNi2_Tm\\energy_YNi2_Tm_x+1.000.dat'),
                         common.get_paths(common.PATH_TO_ENERGY_DATAFILES,
                                          'energy', 'dat', 'YNi2', 'Tm', x_parameter=1))
        self.assertEqual(os.path.join(common.BASE_DIR,
                                      'energy_datafiles\\YNi2_Tm\\energy_YNi2_Tm_x-1.000.dat'),
                         common.get_paths(common.PATH_TO_ENERGY_DATAFILES,
                                          'energy', 'dat', 'YNi2', 'Tm', x_parameter=-1))


if __name__ == '__main__':
    unittest.main()
