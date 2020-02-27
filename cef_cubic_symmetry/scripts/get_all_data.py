"""The module contains the procedure of CEF parameters defining with spectra saving."""
from . import get_results as gr
from . import get_graph as gg

MATERIAL = {'crystal': 'YNi2', 'rare_earth': 'Tm'}
PARAMETERS = {'w': 1}
EXPERIMENTAL_ENERGY_1 = 0.45606
EXPERIMENTAL_ENERGY_2 = 2.09349
EXPERIMENTAL_RATIO = EXPERIMENTAL_ENERGY_2 / EXPERIMENTAL_ENERGY_1
TEMPERATURES = (3, 5, 25)

CROSSES = gr.find_cross(EXPERIMENTAL_RATIO, material=MATERIAL, parameters=PARAMETERS)
RECALCULATED_CROSSES = gr.recalculation(CROSSES, EXPERIMENTAL_ENERGY_1, material=MATERIAL)
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
    for temperature in (3, 5):
        energies, intensities = gr.save_spectra_with_two_temperatures(material=MATERIAL,
                                                                      parameters={'w': point.w,
                                                                                  'x': point.x,
                                                                                  'gamma': 0.16},
                                                                      temperature_1=temperature,
                                                                      temperature_2=25)

        colors = {temperature: 'black', 25: 'red', 'diff': 'green'}
        legend = {temperature: f'{temperature} K', 25: '25 K', 'diff': f'{temperature} - 25 K'}
        gg.get_energy_transfer_plot(material=MATERIAL,
                                    parameters={'w': point.w, 'x': point.x},
                                    temperatures=(temperature, 25),
                                    data=gg.Data(x=energies,
                                                 y_set=intensities,
                                                 legend=legend),
                                    scale=gg.Scale(limits=gg.Limits(x_min=0,
                                                                    x_max=4,
                                                                    y_min=0,
                                                                    y_max=4000),
                                                   locators=gg.Locators(x_major=1,
                                                                        x_minor=0.2,
                                                                        y_major=500,
                                                                        y_minor=100)
                                                   ),
                                    )
        print('Graph saved')
