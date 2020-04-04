from scripts.common.constants import BASE_DIR
from scripts.common.utils import user_input, check_input
from scripts.common.named_tuples import Material
from scripts.get_cef_object import CF
from scripts import get_results
print(f'Working directory: {BASE_DIR}\n')
while True:
    command = input('\nInput command (e.g. "help"): ')
    if command == 'get_one_dot':
        print()
        crystal, rare_earth, w, x = user_input()
        get_results.get_one_dot(material=Material(crystal=crystal,
                                                  rare_earth=rare_earth),
                                parameters={'w': w,
                                            'x': x})
    elif command == 'get_object':
        print()
        crystal, rare_earth, w, x = user_input()
        print(CF(Material(crystal=crystal,
                          rare_earth=rare_earth)).get_llw_cubic_object(parameters={'w': w,
                                                                                   'x': x}))
    elif command == 'load_data':
        print()
        crystal, rare_earth, w, x = user_input()
        cef_object = get_results.load_data(material=Material(crystal=crystal,
                                                             rare_earth=rare_earth))
        print(cef_object)
    elif command == 'save_energy_dat':
        print()
        crystal = input('Input the name of crystal (e.g. "YNi2"): ')
        rare_earth = check_input('rare')
        w = check_input('w')
        number_of_intervals = check_input('intervals')
        get_results.save_energy_dat(crystal, rare_earth, w, number_of_intervals)
    elif command == 'save_parameters':
        print()
        crystal, rare_earth, w, x = user_input()
        get_results.save_parameters(crystal, rare_earth, w, x)
    elif command == 'save_spectra':
        print()
        crystal, rare_earth, w, x = user_input()
        get_results.save_parameters(crystal, rare_earth, w, x)
    elif command == 'save_sus':
        print()
        crystal, rare_earth, w, x = user_input()
        get_results.save_susceptibility(crystal, rare_earth, w, x)
    elif command == 'help':
        print()
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
        print()
        print("If you don't know, which commands can be inputted, enter 'help'")
