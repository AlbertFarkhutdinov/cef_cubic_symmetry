from cef_object_scripts import get_results

# save_info('Start')
print('Hello!')
while True:
    command = input('\nInput command: ')
    if command == 'get_one_dot':
        get_results.get_one_dot()
    elif command == 'load_data':
        pass
    elif command == 'save_energy_dat':
        pass
    elif command == 'save_parameters':
        pass
    elif command == 'save_spectra':
        pass
    elif command == 'save_susceptibility':
        pass
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
