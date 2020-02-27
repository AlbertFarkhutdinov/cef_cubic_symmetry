"""This module contains unittests for module utils.py"""
import unittest
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

    def test_get_new_if_old_is_none(self):
        """Tests for function get_new_if_old_is_none"""
        self.assertEqual(0, utils.get_new_if_old_is_none(0, 10))
        self.assertEqual(10, utils.get_new_if_old_is_none(10, 0))
        self.assertEqual(0, utils.get_new_if_old_is_none(None, 0))
        self.assertEqual(0, utils.get_new_if_old_is_none(0, None))
        self.assertEqual(None, utils.get_new_if_old_is_none(None, None))

    def test_value_to_write(self):
        """Tests for function value_to_write"""
        self.assertEqual('   3.14159\t', utils.value_to_write(3.1415926535, '\t'))
        self.assertEqual('  -3.14159\n', utils.value_to_write(-3.1415926535, '\n'))


if __name__ == '__main__':
    unittest.main()
