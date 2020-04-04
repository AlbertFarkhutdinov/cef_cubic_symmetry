"""The module contains the procedure of CEF parameters defining with spectra saving."""
import os
from scripts.common.constants import PATH_TO_EXPERIMENTAL_DATAFILES
from scripts.common.named_tuples import Material, Data, Scale
from scripts import get_results as gr
from scripts import get_graph as gg
from scripts.fitting import get_data_from_file


SPECTROMETER = 'PSI'
INITIAL_ENERGY = 3
MATERIAL = Material(crystal='YNi2', rare_earth='Er')

EXPERIMENTAL_ENERGIES = (0.6, 1.3)
EXPERIMENTAL_RATIO = EXPERIMENTAL_ENERGIES[1] / EXPERIMENTAL_ENERGIES[0]

gg.get_llw_ratios_plot(material=MATERIAL,
                       experimental_value=EXPERIMENTAL_RATIO,
                       y_limits={'max_value': 3, 'min_value': 1},
                       y_ticks={'y_major': 0.5, 'y_minor': 0.1})


TEMPERATURES = (1.4, 15)
DATA = []
for temperature in TEMPERATURES:
    DATA.append(get_data_from_file(os.path.join(PATH_TO_EXPERIMENTAL_DATAFILES,
                                                '_'.join([SPECTROMETER,
                                                          MATERIAL.rare_earth,
                                                          MATERIAL.crystal,
                                                          f'{INITIAL_ENERGY}meV',
                                                          f'{temperature}K.txt'])
                                                )
                                   )
                )

gg.get_spectrum_experiment(material=MATERIAL,
                           temperatures=TEMPERATURES,
                           data=tuple(DATA),
                           scale=Scale(
                               limits={
                                   'x_min': 0,
                                   'x_max': 2.2,
                                   'y_min': 0,
                                   'y_max': 20,
                               },
                               locators={
                                   'x_major': 1,
                                   'x_minor': 0.2,
                                   'y_major': 5,
                                   'y_minor': 1,
                               }
                           ),
                           )


PARAMETERS = {'w': 1}
CROSSES = gr.find_cross(EXPERIMENTAL_RATIO, material=MATERIAL, parameters=PARAMETERS)
RECALCULATED_CROSSES = gr.recalculation(CROSSES, EXPERIMENTAL_ENERGIES[0], material=MATERIAL)
print('Cross points:')
for point in RECALCULATED_CROSSES:
    print(f'w = {point.w : 6.3f};\tx = {point.x : .3f};')

for point in RECALCULATED_CROSSES:
    for temperature in TEMPERATURES:
        gr.save_spectra_with_one_temperature(material=MATERIAL,
                                             parameters={'w': point.w,
                                                         'x': point.x,
                                                         'gamma': 0.16,
                                                         'T': temperature})

    energies, intensities = gr.save_spectra_with_two_temperatures(material=MATERIAL,
                                                                  parameters={'w': point.w,
                                                                              'x': point.x,
                                                                              'gamma': 0.16},
                                                                  temperature_1=TEMPERATURES[0],
                                                                  temperature_2=TEMPERATURES[1])

    gg.get_spectrum_theory(material=MATERIAL,
                           parameters={'w': point.w, 'x': point.x},
                           temperatures=TEMPERATURES,
                           data=Data(x=energies,
                                     y_set=intensities,
                                     errors=None,
                                     legend={TEMPERATURES[0]: f'{TEMPERATURES[0]} K',
                                             TEMPERATURES[1]: f'{TEMPERATURES[1]} K',
                                             'diff': f'{TEMPERATURES[0]} K - {TEMPERATURES[1]} K'}),
                           scale=Scale(
                               limits={
                                   'x_min': 0,
                                   'x_max': 2.2,
                                   'y_min': 0,
                                   'y_max': 4000,
                               },
                               locators={
                                   'x_major': 1,
                                   'x_minor': 0.2,
                                   'y_major': 500,
                                   'y_minor': 100,
                               }
                           ),
                           )
