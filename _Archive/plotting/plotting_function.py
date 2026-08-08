def plotting_temperature(data_x, data_y, units):
    plt.plot(data_x, data_y, linestyle='-', marker='o', lw=0.5)
    plt.grid(True, lw=0.1, color='grey')
    temp, length = ['^{o}F', 'ft']
    if units.lower() == 'si':
        temp, length = ['^{o}C', 'm']
    plt.ylabel(fr'Temperature ${temp}$')
    plt.xlabel(fr'Length ${length}$')
