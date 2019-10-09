from cef_object_scripts import core
from cef_object_scripts import get_results

# save_info('Start')
print('Hello!')
crystal = input('Input the name of crystal: ')
while True:
    command = input('\nInput command: ')
    if command == 'get_one_dot':
        rare_earth = core.check_input('rare')
        w = core.check_input('w')
        x = core.check_input('x')
        get_results.get_one_dot(w, x, rare_earth)
    elif command == 'load_data':
        rare_earth = core.check_input('rare')
        w = core.check_input('w')
        x = core.check_input('x')
        cef_object = get_results.load_data(crystal, rare_earth, w, x)
        print(cef_object)
    elif command == 'save_energy_dat':
        rare_earth = core.check_input('rare')
        w = core.check_input('w')
        number_of_intervals = core.check_input('intervals')
        get_results.save_energy_dat(crystal, rare_earth, w, number_of_intervals)
    elif command == 'save_parameters':
        rare_earth = core.check_input('rare')
        w = core.check_input('w')
        x = core.check_input('x')
        get_results.save_parameters(crystal, rare_earth, w, x)
    elif command == 'save_spectra':
        rare_earth = core.check_input('rare')
        w = core.check_input('w')
        x = core.check_input('x')
        get_results.save_parameters(crystal, rare_earth, w, x)
    elif command == 'save_sus':
        rare_earth = core.check_input('rare')
        w = core.check_input('w')
        x = core.check_input('x')
        get_results.save_susceptibility(crystal, rare_earth, w, x)
    elif command == 'help':
        print('Help:')
        print('get_one_dot - command, which ')
        print('load_data - command, which ')
        print('save_energy_dat - command, which ')
        print('save_parameters - command, which ')
        print('save_spectra - command, which ')
        print('save_susceptibility - command, which ')
    elif command == 'exit':
        break
    else:
        print("If you don't know, which commands can be inputted, enter 'help'")

# save_info('Finish')
