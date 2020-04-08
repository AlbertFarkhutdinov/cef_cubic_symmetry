"""This module contains unittests for module utils.py"""
import unittest
from sys import path as sys_path
from os import getcwd
from os.path import join
sys_path.append(join(getcwd(), '..'))
from scripts.common import utils


class UtilsTests(unittest.TestCase):
    """Class with tests for utils.py"""
    def test_get_sign(self):
        """Tests for function get_sign"""
        self.assertEqual('+', utils.get_sign(1))
        self.assertEqual('+', utils.get_sign(0))
        self.assertEqual('-', utils.get_sign(-1))

    def test_get_value_with_sign(self):
        """Tests for function get_value_with_sign"""
        self.assertEqual('+1.000', utils.get_value_with_sign(1))
        self.assertEqual('+0.346', utils.get_value_with_sign(0.3456))
        self.assertEqual('-1.753', utils.get_value_with_sign(-1.75321))

    def test_get_default(self):
        """Tests for function get_default"""
        self.assertEqual(0, utils.get_default(0, 10))
        self.assertEqual(10, utils.get_default(10, 0))
        self.assertEqual(0, utils.get_default(None, 0))
        self.assertEqual(0, utils.get_default(0, None))
        self.assertEqual(None, utils.get_default(None, None))


if __name__ == '__main__':
    unittest.main()
