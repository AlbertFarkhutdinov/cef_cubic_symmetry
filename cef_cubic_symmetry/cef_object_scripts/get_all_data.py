from cef_object_scripts import get_results, get_graph

crystal_name = 'YNi2'
rare_earth_name = 'Tm'
w = 1
experimental_energy_1 = 0.45606
experimental_energy_2 = 2.09349
experimental_ratio = experimental_energy_2 / experimental_energy_1
temperatures = [3, 5, 25]

crosses = get_results.find_cross(experimental_ratio, crystal_name, rare_earth_name, w)
recalculated_crosses = get_results.recalculation(crosses, experimental_energy_1, crystal_name, rare_earth_name)
print('Cross points:')
for point in recalculated_crosses:
    print(f'w = {point.w : 6.3f};\tx = {point.x : .3f};')

for point in recalculated_crosses:
    for temperature in temperatures:
        get_results.save_spectra_with_one_temperature(crystal_name, rare_earth_name, point.w, point.x, temperature,
                                                      gamma=0.16)
    for temperature in (3, 5):
        energies, intensities = get_results.save_spectra_with_two_temperatures(crystal=crystal_name,
                                                                               rare_earth=rare_earth_name,
                                                                               w=point.w,
                                                                               x=point.x,
                                                                               temperature_1=temperature,
                                                                               temperature_2=25)

        colors = {temperature: 'black', 25: 'red', 'diff': 'green'}
        legend = {temperature: f'{temperature} K', 25: '25 K', 'diff': f'{temperature} - 25 K'}
        get_graph.get_energy_transfer_plot(crystal=crystal_name, rare_earth=rare_earth_name,
                                           w=point.w, x=point.x,
                                           temperature_1=temperature, temperature_2=25,
                                           x_data=energies, y_data_set=intensities,
                                           color_set=colors, y_legend_set=legend, mode='png')
        print('Graph saved')
