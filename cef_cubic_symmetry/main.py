from scripts.common.constants import BASE_DIR, Material
from scripts.common.utils import user_input, check_input
from scripts.cubic_cef_object import Cubic
print(f'Working directory: {BASE_DIR}\n')
while True:
    command = input('\nInput command (e.g. "help"): ')
    print()
    if command == 'get_one_dot':
        crystal, rare_earth, w, x = user_input()
        Cubic(
            material=Material(
                crystal=crystal,
                rare_earth=rare_earth
            ),
            llw_parameters={
                'w': w,
                'x': x,
            }
        ).get_one_dot()
    elif command == 'get_object':
        crystal, rare_earth, w, x = user_input()
        print(
            Cubic(
                material=Material(
                    crystal=crystal,
                    rare_earth=rare_earth
                ),
                llw_parameters={
                    'w': w,
                    'x': x,
                }
            )
        )
    elif command == 'load_data':
        crystal, rare_earth, w, x = user_input()
        cef_object = Cubic(
                material=Material(
                    crystal=crystal,
                    rare_earth=rare_earth
                ),
                llw_parameters={
                    'w': w,
                    'x': x,
                }
            )
        cef_object.load_data()
        print(cef_object)
    elif command == 'save_energy_dat':
        crystal = input('Input the name of crystal (e.g. "YNi2"): ')
        rare_earth = check_input('rare')
        w = check_input('w')
        number_of_intervals = check_input('intervals')
        Cubic(
            material=Material(
                crystal=crystal,
                rare_earth=rare_earth
            ),
            llw_parameters={
                'w': w,
            }
        ).save_energy_dat(number_of_intervals)
    elif command == 'save_parameters':
        crystal, rare_earth, w, x = user_input()
        Cubic(
            material=Material(
                crystal=crystal,
                rare_earth=rare_earth
            ),
            llw_parameters={
                'w': w,
                'x': x,
            }
        ).save_to_file()
    elif command == 'save_spectra':
        crystal, rare_earth, w, x = user_input()
        Cubic(
            material=Material(
                crystal=crystal,
                rare_earth=rare_earth
            ),
            llw_parameters={
                'w': w,
                'x': x,
            }
        ).save_spectra_with_one_temperature(gamma=0.16, temperature=3)
    elif command == 'save_sus':
        crystal, rare_earth, w, x = user_input()
        Cubic(
            material=Material(
                crystal=crystal,
                rare_earth=rare_earth
            ),
            llw_parameters={
                'w': w,
                'x': x,
            }
        ).save_susceptibility()
    elif command == 'help':
        print('\nList of available commands:')
        print('get_object - prints CEF-object for values inputted by user.')
        print('get_one_dot - prints energy transfers for values inputted by user.')
        print('save_parameters - saves parameters of CEF to .cgf-file for values inputted by user.')
        print('load_data - loads CEF-parameters from .cfg-file for values inputted by user.')
        print('save_energy_dat - saves dependence of energy transfer on x-parameter of CEF to .dat-file for values '
              'inputted by user.')
        print('save_spectra - saves neutron inelastic scattering spectra to .dat-file for values inputted by user.')
        print('save_sus - saves magnetic susceptibilities to .dat-file for values inputted by user.')
    elif command == 'exit':
        break
    else:
        print("If you don't know, which commands can be inputted, enter 'help'")
