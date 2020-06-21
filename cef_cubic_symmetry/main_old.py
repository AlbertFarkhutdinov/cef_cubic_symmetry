"""The module contains the console manager for work with project."""
from scripts.common.constants import BASE_DIR, Material
from scripts.common.utils import check_input
from scripts.cubic_cef_object import Cubic


def main():
    """Main procedure"""
    print(f'Working directory: {BASE_DIR}\n')
    while True:
        command = input('\nInput command (e.g. "help"): ')
        print()
        kwargs = {
            'material': None,
            'llw_parameters': {},
        }
        if command not in ('help', 'exit'):
            crystal = input('Input the name of crystal (e.g. "YNi2"): ')
            rare_earth = check_input('rare')
            kwargs['material'] = Material(
                crystal=crystal,
                rare_earth=rare_earth
            ),
            kwargs['llw_parameters'] = {
                'w': check_input('w'),
            }
        if command != 'save_peak_dat':
            kwargs['llw_parameters']['x'] = check_input('x')

        methods = {
            'get_object': None,
            'get_one_dot': 'get_one_dot',
            'load_data': 'load_data',
            'save_peak_dat': 'save_peak_dat',
            'save_parameters': 'save_to_file',
            'save_spectra': 'save_spectra_with_one_temperature',
            'save_sus': 'save_susceptibility',
        }
        if command in methods:
            cef_object = Cubic(**kwargs)
            kwargs = {}
            if command == 'save_peak_dat':
                kwargs['number_of_intervals'] = check_input('intervals')
            elif command == 'save_spectra':
                kwargs['gamma'] = 0.16
                kwargs['temperature'] = 3
            if command is not None:
                cef_object.__getattribute__(command)(**kwargs)
            if command in ('get_object', 'load_data'):
                print(cef_object)
        elif command == 'help':
            print(
                '\nList of available commands:',
                'get_object - prints CEF-object for values inputted by user.',
                'get_one_dot - prints energy transfers for values inputted by user.',
                'save_parameters - saves parameters of CEF to '
                '.cgf-file for values inputted by user.',
                'load_data - loads CEF-parameters from .cfg-file for values inputted by user.',
                'save_peak_dat - saves dependence of energy transfer '
                'on x-parameter of CEF to .dat-file for values '
                'inputted by user.',
                'save_spectra - saves neutron inelastic scattering spectra '
                'to .dat-file for values inputted by user.',
                'save_sus - saves magnetic susceptibilities to '
                '.dat-file for values inputted by user.',
                sep='\n'
            )
        elif command == 'exit':
            break
        else:
            print("If you don't know, which commands can be inputted, enter 'help'")


if __name__ == '__main__':
    main()
