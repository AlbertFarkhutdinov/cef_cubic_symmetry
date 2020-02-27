"""The module contains the procedure of CEF parameters defining with spectra saving."""
from os.path import join
from .common.constants import PATH_TO_EXPERIMENTAL_DATAFILES
from . import get_results as gr
from . import get_graph as gg
from .fitting import get_data_from_file

MATERIAL = {'crystal': 'YNi2', 'rare_earth': 'Er'}
PARAMETERS = {'w': 1}
EXPERIMENTAL_ENERGY_1 = 0.6
EXPERIMENTAL_ENERGY_2 = 1.3
EXPERIMENTAL_RATIO = EXPERIMENTAL_ENERGY_2 / EXPERIMENTAL_ENERGY_1
TEMPERATURE_1 = 1.4
TEMPERATURE_2 = 15
DATA_1 = get_data_from_file(join(PATH_TO_EXPERIMENTAL_DATAFILES,
                                 f'PSI_{MATERIAL["rare_earth"]}_YNi2_3meV_{TEMPERATURE_1}K.txt'))
DATA_2 = get_data_from_file(join(PATH_TO_EXPERIMENTAL_DATAFILES,
                                 f'PSI_{MATERIAL["rare_earth"]}_YNi2_3meV_{TEMPERATURE_2}K.txt'))
LEGEND = {TEMPERATURE_1: f'{TEMPERATURE_1} K',
          TEMPERATURE_2: f'{TEMPERATURE_2} K',
          'diff': f'{TEMPERATURE_1} K - {TEMPERATURE_2} K'}

gg.get_energy_transfer_plot(material=MATERIAL,
                            parameters={'w': 0, 'x': 0},
                            temperatures=(TEMPERATURE_1, TEMPERATURE_2),
                            data=gg.Data(x=DATA_1[0],
                                         y_set={TEMPERATURE_1: DATA_1[1],
                                                TEMPERATURE_2: DATA_2[1],
                                                'diff': DATA_1[1] - DATA_2[1],
                                                },
                                         legend=LEGEND),
                            scale=gg.Scale(limits=gg.Limits(x_min=0,
                                                            x_max=2.2,
                                                            y_min=0,
                                                            y_max=20),
                                           locators=gg.Locators(x_major=1,
                                                                x_minor=0.2,
                                                                y_major=5,
                                                                y_minor=1)
                                           ),
                            )


CROSSES = gr.find_cross(EXPERIMENTAL_RATIO, material=MATERIAL, parameters=PARAMETERS)
RECALCULATED_CROSSES = gr.recalculation(CROSSES, EXPERIMENTAL_ENERGY_1, material=MATERIAL)
print('Cross points:')
for point in RECALCULATED_CROSSES:
    print(f'w = {point.w : 6.3f};\tx = {point.x : .3f};')

for point in RECALCULATED_CROSSES:
    for temperature in (TEMPERATURE_1, TEMPERATURE_2):
        gr.save_spectra_with_one_temperature(material=MATERIAL,
                                             parameters={'w': point.w,
                                                         'x': point.x,
                                                         'gamma': 0.16,
                                                         'T': temperature})

    energies, intensities = gr.save_spectra_with_two_temperatures(material=MATERIAL,
                                                                  parameters={'w': point.w,
                                                                              'x': point.x,
                                                                              'gamma': 0.16},
                                                                  temperature_1=TEMPERATURE_1,
                                                                  temperature_2=TEMPERATURE_2)

    gg.get_energy_transfer_plot(material=MATERIAL,
                                parameters={'w': point.w, 'x': point.x},
                                temperatures=(TEMPERATURE_1, TEMPERATURE_2),
                                data=gg.Data(x=energies,
                                             y_set=intensities,
                                             legend=LEGEND),
                                scale=gg.Scale(limits=gg.Limits(x_min=0,
                                                                x_max=2.2,
                                                                y_min=0,
                                                                y_max=4000),
                                               locators=gg.Locators(x_major=1,
                                                                    x_minor=0.2,
                                                                    y_major=500,
                                                                    y_minor=100)
                                               ),
                                )
    print('Graph saved')
