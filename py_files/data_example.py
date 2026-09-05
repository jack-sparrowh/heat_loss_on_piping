import numpy as np
import pandas as pd

def user_table_example():

    # USER TABLE EXAMPLE
    user_array_1 = np.array([
        [4, 3.8, 0.1, 500] for i in range(4)
    ])
    
    user_array_2 = np.array([
        [4, 3.8, 0.1, 10] for i in range(200)
    ])
    
    user_array_3 = np.array([
        [4, 3.8, 0.1, 2000]
    ])
    
    user_array_imperial = np.array([
        [4, 3.8, 0.5, 831],
        [3, 2.8, 0.5, 2727],
        [2, 1.8, 0.5, 207]
    ])

    user_array_si = np.array([
        [4, 3.8, 0.5, 253],
        [3, 2.8, 0.5, 831],
        [2, 1.8, 0.5, 63]
    ])

    return [user_array_1, user_array_2, user_array_3, user_array_imperial, user_array_si]

def properties_table_example():

    # PROPERTIES TABLE EXAMPLE
    properties_table_imperial = pd.Series(
        [1800, 68, 60, 9.7, 0.025, 24155, 1, 50, 1], 
        index=[
            'inlet_temperature', 
            'ambient_temperature', 
            'h_fluid_to_pipe', 
            'k_pipe', 
            'k_insulation', 
            'mass_flow_rate', 
            'specific_heat_capacity', 
            'solar_flux', 
            'emissivity'
        ])

    properties_table_si = pd.Series(
        [982, 20, 340, 16.8, 0.043, 10942, 4.2, 157.7, 1], 
        index=[
            'inlet_temperature', 
            'ambient_temperature', 
            'h_fluid_to_pipe', 
            'k_pipe', 
            'k_insulation', 
            'mass_flow_rate', 
            'specific_heat_capacity', 
            'solar_flux', 
            'emissivity'
        ])

    return [properties_table_imperial, properties_table_si]