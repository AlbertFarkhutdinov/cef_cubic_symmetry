"""The module contains the procedure of CEF parameters defining with spectra saving."""
from cef_object_scripts import get_results as gr
from cef_object_scripts import get_graph as gg

CRYSTAL_NAME = 'YNi2'
RARE_EARTH_NAME = 'Tm'
W = 1
EXPERIMENTAL_ENERGY_1 = 0.45606
EXPERIMENTAL_ENERGY_2 = 2.09349
EXPERIMENTAL_RATIO = EXPERIMENTAL_ENERGY_2 / EXPERIMENTAL_ENERGY_1
TEMPERATURES = (3, 5, 25)

CROSSES = gr.find_cross(EXPERIMENTAL_RATIO, CRYSTAL_NAME, RARE_EARTH_NAME, W)
RECALCULATED_CROSSES = gr.recalculation(CROSSES, EXPERIMENTAL_ENERGY_1,
                                        CRYSTAL_NAME, RARE_EARTH_NAME)
print('Cross points:')
for point in RECALCULATED_CROSSES:
    print(f'w = {point.w : 6.3f};\tx = {point.x : .3f};')

for point in RECALCULATED_CROSSES:
    for temperature in TEMPERATURES:
        gr.save_spectra_with_one_temperature(CRYSTAL_NAME,
                                             RARE_EARTH_NAME,
                                             point.w, point.x,
                                             temperature,
                                             gamma=0.16)
    for temperature in (3, 5):
        energies, intensities = gr.save_spectra_with_two_temperatures(crystal=CRYSTAL_NAME,
                                                                      rare_earth=RARE_EARTH_NAME,
                                                                      w_parameter=point.w,
                                                                      x_parameter=point.x,
                                                                      temperature_1=temperature,
                                                                      temperature_2=25)

        colors = {temperature: 'black', 25: 'red', 'diff': 'green'}
        legend = {temperature: f'{temperature} K', 25: '25 K', 'diff': f'{temperature} - 25 K'}
        gg.get_energy_transfer_plot(crystal=CRYSTAL_NAME,
                                    rare_earth=RARE_EARTH_NAME,
                                    w_parameter=point.w, x_parameter=point.x,
                                    x_min=0, x_max=4, y_max=4000,
                                    text_x=2.2, text_y=2300,
                                    temperature_1=temperature, temperature_2=25,
                                    x_data=energies, y_data_set=intensities,
                                    x_major_locator=1, y_major_locator=500,
                                    color_set=colors, y_legend_set=legend, mode='png')
        print('Graph saved')
