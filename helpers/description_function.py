def notes(desc_file):
    description = desc_file['description'] + '\n'
    parameters = '\nParameters\n----------\n'
    params_group = desc_file['parameters']
    for param in params_group:
        current_param = params_group[param]
        parameters = parameters + param + ': ' + current_param[0] + '\n\t' + current_param[1] + f' [{current_param[2]}]' + '\n'
    returned = '\nReturns\n-------\n' + desc_file['returns'][0] + ': ' + desc_file['returns'][1] + '\n\t' + desc_file['returns'][2]
    full_note = description + parameters + returned
    return full_note

temp_func = {
    'description': 'Calculate fluid temperature at provided length and combined heat transfer coefficient.',
    'parameters': {
        'length': ['float, array_like', 'Length at which temperature is to be calculated.'],
        'h_combined_surface': ['float, array_like', 'Combined heat transfer coefficient at the insulation-air interface.'],
        'inlet_temperature': ['float', 'Temperature of fluid at length = 0.'],
        'ambient_temperature': ['float', 'Ambient temperature around piping.'],
        'densit': ['float', 'Density of a fluid.'],
        'volumetric_flow_rat': ['float', 'Volumetric flow rate of fluid.'],
        'specific_heat_capacity': ['float', 'Specific heat capacity of fluid.']
    },
    'returns': [
        'T_at_length', 'float, array_like', 'Temperature calculated at provided length.'
    ]
}

cp_func = {
    'description': 'Calculate specific heat capacity of air at reference temperature, which is film temperature of air around piping.\nModel is taken from Cengel, and is fitted to a range 273K - 1800K. Below 273K cp is equal to cp at 273K and above 1800K cp is equal to cp at 1800K.',
    'parameters': {
        'ref_temperature': ['float, array_like', 'Reference temperature at insulation-air interface.']
    },
    'returns': [
        '_cp', 'float, array_like', 'Specific heat capacity.'
    ]
}

u_func = {
    'description': 'Calculate kinematic viscosity of air at reference temperature, which is film temperature of air around piping.\nModel is taken from ???.',
    'parameters': {
        'ref_temperature': ['float, array_like', 'Reference temperature at insulation-air interface.', 'K']
    },
    'returns': [
        '_u', 'float, array_like', 'Kinematic viscosity [Pa*s]'
    ]
}

Ra_func = {
    'description': 'Calculate Rayleigh number of air at film temperature for horizontal cylindrical body. It is assummed that Pr=0.71 for air in all temperature range.',
    'parameters': {
        'surface_temperature': ['float, array_like', 'Insulation surface temperature', 'K'],
        'ambient_temperature': ['float', 'Ambient temperature around piping', 'K'],
        'length': ['float, array_like', 'Length to determine outer diameter of piping', 'ft/m - units_handler determines units']
    },
    'returns': [
        '_Ra', 'float, array_like', 'Rayleigh number.'
    ]
}

h_func = {
    'description': 'Calculate combined heat transfer coefficient at insulation-air interface.\nRadiation heat transfer coeff. is in its typical form.\nConvective heat transfer is based on Churchill and Chu correlation formula.',
    'parameters': {
        'length': ['float, array_like', 'Length at which combined heat transfer coeff. is to be determined', 'ft'],
        'surface_temperature': ['float, array_like', 'Insulation surface temperature', 'F'],
        'ambient_temperature': ['float', 'Ambient temperature around piping', 'F'],
        'insulation_thickness': ['float', 'Thickness of insulation', 'inch'],
        'emissivity': ['float', 'Emissivity of insulation material or material covering insulaiton.', '']
    },
    'returns': [
        '_h_comb_vec', 'array_like', 'Array with heat transfer coefficients of radiative and convective effects.'
    ]
}

Q_th_p_func = {
    'description': 'Calculate total heat transfer rate per unit length through whole piping cross section.',
    'parameters': {
        'length': ['float, array_like', 'Length at which total heat transfer rate is to be determined', 'ft'],
        'h_combined_surface': ['float, array_like, array_like', 'Combined heat transfer coefficient', 'BTU/ft2.h.F'],
        'ambient_temperature': ['float', 'Ambient temperature around piping', 'F']
    },
    'returns': [
        '_Q_through', 'float, array_like', 'Total heat transfer rate per unit length.'
    ]
}

Q_to_p_func = {
    'description': 'Calculate total heat transfer rate per unit length transferred from the fluid to insulation.',
    'parameters': {
        'length': ['float, array_like', 'Length at which total heat transfer rate is to be determined', 'ft'],
        'h_combined_surface': ['float, array_like, array_like', 'Combined heat transfer coefficient', 'BTU/ft2.h.F'],
        'surface_temperature': ['float, array_like', 'Insulation surface temperature', 'F']
    },
    'returns': [
        '_Q_through', 'float, array_like', 'Total heat transfer rate per unit length from fluid to insulation.'
    ]
}

Q_out_func = {
    'description': 'Calculate total heat transfer rate per unit length transferred from the insulation surface to the air surrounding piping.',
    'parameters': {
        'length': ['float, array_like', 'Length at which total heat transfer rate is to be determined', 'ft'],
        'h_combined_surface': ['float, array_like, array_like', 'Combined heat transfer coefficient', 'BTU/ft2.h.F'],
        'surface_temperature': ['float, array_like', 'Insulation surface temperature', 'F'],
        'ambient_temperature': ['float', 'Ambient temperature around piping', 'F']
    },
    'returns': [
        '_Q_through', 'float, array_like', 'Total heat transfer rate per unit length from insulation to air.'
    ]
}

loss_func = {
    'description': 'Calculates squared difference between total heat transfer and heat transfer to the outer insulation surface at given surface temperature and length.',
    'parameters': {
        'x':['list, array_like', 'Surface temperature', 'F'],
        'length': ['float, array_like', 'Length at which total heat transfer rate is to be determined', 'ft']
    },
    'returns': [
        'loss', 'float, array_like', 'Squared difference.'
    ]
}

calc_tab_func = {
    'description': 'Prepare a semi-final dataframe containing values used by final function to create user-redable table.\nCreates a pandas dataframe with cumulative length, surface temperature minimizing loss function,\nheat transfer coefficients for both radiative heat transfer and convective and calculated temperature of fluid at specified length.',
    'parameters': {
        'mode':['str', 'Mode of calculation. Either "conservative" or "exact"', ''], 
        'length_series':['array_like', 'Series or array of cummulative lengths', ''], 
        'inlet_temperature': ['float', 'Temperature of fluid at length = 0', 'F'], 
        'ambient_temperature': ['float', 'Ambient temperature around piping', 'F']
    },
    'returns': [
        '_semi_final_table', 'pd.DataFrame', 'Semi-final table.'
    ]
}

usr_prepare_table_func = {
    'description': 'Take Nx4 datafame and make it Mx4 where M is the amount of rows completely filled by user. Excel treats empty cells as 0, and we don\'t expect user to provide 0 at any point.',
    'parameters': {
        'user_table':['pd.DataFrame', 'User table', '']
    },
    'returns': [
        'prep_user_table', 'pd.DataFrame', 'Prepared User Table.'
    ]
}

usr_expand_table_func = {
    'description': 'Extend prepared user table by cumulative length, and handle length units.',
    'parameters': {
        'user_table':['pd.DataFrame', 'User table', ''],
        'units':['str', 'Either "SI" or "Imperial"', '']
    },
    'returns': [
        'ext_prep_user_table', 'pd.DataFrame', 'Extended User Table.'
    ]
}

prop_prepare_table_func = {
    'description': 'Simple handling of properties table. Take 1x9 dataframe and add names to indicies and convert to series for ease of calling.',
    'parameters': {
        'properties_table':['pd.DataFrame', 'Properties table', '']
    },
    'returns': [
        'prep_properties_table', 'pd.DataFrame', 'Prepared properties table.'
    ]
}

units_handler_func = {
    'description': 'Converts all SI units to Impertial if units is SI, otherwise does nothing.',
    'parameters': {
        'properties_table':['pd.DataFrame', 'Properties table', ''],
        'units':['str', 'Either "SI" or "Imperial"', '']
    },
    'returns': [
        'units_handled_prop_table', 'pd.DataFrame', 'Prepared properties table, with handled units.'
    ]
}

user_units_func = {
    'description': 'Prints out all necessary units based on choosen metric system.',
    'parameters': {
        'units':['str', 'Either "SI" or "Imperial"', '']
    },
    'returns': [
        'units_list', 'list', 'List of units.'
    ]
}

segment_print_func = {
    'description': 'Prints out units of length based on choosen metric system.',
    'parameters': {
        'units':['str', 'Either "SI" or "Imperial"', '']
    },
    'returns': [
        'string', 'str', 'String of length units.'
    ]
}

final_table_func = {
    'description': 'Creates final table before units conversion. It calculates all selected heats to show the difference between heat to the insulation and from the insulation as a sum of radiation and convection. It also attaches simple to folow column names.',
    'parameters': {
        'mode':['str', 'Mode of calculation. Either "conservative" or "exact"', '']
    },
    'returns': [
        '_final_table', 'pd.DataFrame', 'Final table ready for unit conversion.'
    ]
}


final_table_units_handler_func = {
    'description': 'Converts final table to correct units.',
    'parameters': {
        'mode':['str', 'Mode of calculation. Either "conservative" or "exact"', ''],
        'units':['str', 'Either "SI" or "Imperial"', '']
    },
    'returns': [
        '_uh_final_table', 'pd.DataFrame', 'Final table.'
    ]
}
