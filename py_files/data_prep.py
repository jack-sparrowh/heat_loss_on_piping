import numpy as np
import pandas as pd

def user_table_prepare(user_table, units):
    
    '''
    Take Nx4 datafame and make it Mx4 where M is the amount of rows completely filled by user. 
    Excel treats empty cells as 0, and we don't expect user to provide whole row filled with 0 at any point.

    Parameters
    ----------
    user_table: pd.DataFrame
        User table []
    units: str
        Either "SI" or "Imperial" []
    
    Returns
    -------
    prep_user_table: pd.DataFrame
        Prepared User Table.
    '''
    
    user_table.columns = [
        'pipe_outer_diameter', 
        'pipe_inner_diameter', 
        'insulation_thickness', 
        'segment_length'
    ]
    prep_user_table = user_table[user_table.notnull().all(1)].copy()

    if units.lower() == 'si':
        prep_user_table.loc[:, 'segment_length'] = prep_user_table.loc[:, 'segment_length'] * 1 / 0.3048 # 1 ft = 0.3048 m
    
    return prep_user_table
    
# def usr_expand_table(user_table, units):
#     
#     '''
#     Extend prepared user table by cumulative length, and handle length units.
# 
#     Parameters
#     ----------
#     user_table: pd.DataFrame
#         User table []
#     units: str
#         Either "SI" or "Imperial" []
#     
#     Returns
#     -------
#     ext_prep_user_table: pd.DataFrame
#         Extended User Table.
#     '''
#     
#     prep_user_table = usr_prepare_table(user_table, units)
#     
#     if units.lower() == 'si':
#         prep_user_table.loc[:, 'segment_length'] = prep_user_table.loc[:, 'segment_length'] * 1 / 0.3048 # 1 ft = 0.3048 m
# 
#     ext_prep_user_table = prep_user_table.copy()
#     return ext_prep_user_table
    
def prop_table_prepare(properties_table):
    
    '''
    Simple handling of properties table. Take 1x9 dataframe and add names to indicies and convert to series for ease of calling.

    Parameters
    ----------
    properties_table: pd.DataFrame
        Properties table []
    
    Returns
    -------
    prep_properties_table: pd.DataFrame
        Prepared properties table. 
    '''
    
    prep_properties_table = properties_table.copy().iloc[:, 0]
    prep_properties_table.name = 'Values'
    prep_properties_table.index = [
        'inlet_temperature',
		'ambient_temperature',
        'h_fluid_to_pipe',
        'k_pipe',
        'k_insulation',
        'mass_flow_rate',
        'specific_heat_capacity',
        'solar_flux',
        'emissivity'
    ]
    return prep_properties_table

def prop_table_units_handler(properties_table, units):
    
    '''
    Converts all SI units to Impertial if units is SI, otherwise does nothing.

    Parameters
    ----------
    properties_table: pd.DataFrame
        Properties table []
    units: str
        Either "SI" or "Imperial" []
    
    Returns
    -------
    units_handled_prop_table: pd.DataFrame
        Prepared properties table, with handled units.
    '''
    
    _uh_prop_table = prop_table_prepare(properties_table)
    
    if units.lower() == 'si':
        _uh_prop_table['inlet_temperature'] = 1.8 * _uh_prop_table['inlet_temperature'] + 32
        _uh_prop_table['ambient_temperature'] = 1.8 * _uh_prop_table['ambient_temperature'] + 32
        _uh_prop_table['h_fluid_to_pipe'] = _uh_prop_table['h_fluid_to_pipe'] * 1 / 5.67826 # 1 BTU/ft2.h.F = 5.67826 W/m2.K
        _uh_prop_table['k_pipe'] = _uh_prop_table['k_pipe'] * 1 / 1.73073 # 1 BTU/ft.h.F = 1.73073 W/m.K
        _uh_prop_table['k_insulation'] = _uh_prop_table['k_insulation'] * 1 / 1.73073 # 1 BTU/ft.h.F = 1.73073 W/m.K
        _uh_prop_table['mass_flow_rate'] = _uh_prop_table['mass_flow_rate'] * 1 / 0.453592 # 1 lb/h = 0.453592 kg/h
        _uh_prop_table['specific_heat_capacity'] = _uh_prop_table['specific_heat_capacity'] * 1 / 4.1868 # 1 BTU/lb.F = 4.1868 kJ/kg.K
        _uh_prop_table['solar_flux'] = _uh_prop_table['solar_flux'] * 1 / 3.154105 # 1 BTU/ft2.h = 3.154105 W/m2
    
    return _uh_prop_table

def user_units(units):
    
    '''
    Prints out all necessary units based on choosen metric system.

    Parameters
    ----------
    units: str
        Either "SI" or "Imperial" []
    
    Returns
    -------
    units_list: list
        List of units.
    '''
    
    units_list = [
        'F',
        'F',
        'BTU/ft2.h.F',
        'BTU/ft.h.F',
        'BTU/ft.h.F',
        'lbm/h',
        'BTU/lbm.F',
        'BTU/ft2.h',
        ''
    ]
    
    if units.lower() == 'si':
        units_list = [
            'C',
            'C',
            'W/m2.K',
            'W/m.K',
            'W/m.K',
            'kg/h',
            'kJ/kg.K',
            'W/m2',
            ''
        ]
    
    return units_list
    
def segment_length_print(units):
    
    '''
    Prints out units of length based on choosen metric system.

    Parameters
    ----------
    units: str
        Either "SI" or "Imperial" []
    
    Returns
    -------
    string: str
        String of length units.
    '''
    length = 'ft'
    if units.lower() == 'si':
        length = 'm'
    
    return f'Segment length [{length}]'