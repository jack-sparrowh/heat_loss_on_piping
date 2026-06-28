def usr_prepare_table(user_table):

    '''
    Take Nx4 datafame and make it Mx4 where M is the amount of rows completely filled by user. 
    Excel treats empty cells as 0, and we don't expect user to provide 0 at any point.

    Parameters
    ----------
    user_table: pd.DataFrame
        User table []

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
    return prep_user_table

def usr_expand_table(user_table, units):

    '''
    Extend prepared user table by cumulative length, and handle length units.

    Parameters
    ----------
    user_table: pd.DataFrame
        User table []
    units: str
        Either "SI" or "Imperial" []

    Returns
    -------
    ext_prep_user_table: pd.DataFrame
        Extended User Table.
    '''

    prep_user_table = usr_prepare_table(user_table)

    # add first row with segment_length of 0
    prep_user_table = pd.concat(
        [
            prep_user_table.iloc[0, :].to_frame().T,
            prep_user_table
        ]
    ).reset_index(drop=True)

    if units.lower() == 'si':
        prep_user_table.loc[:, 'segment_length'] = prep_user_table.loc[:, 'segment_length'] * 1 / 0.3048 # 1 ft = 0.3048 m

    prep_user_table.loc[0, 'segment_length'] = 1e-6
    prep_user_table['cumulative_length'] = prep_user_table.loc[:, 'segment_length'].cumsum()
    ext_prep_user_table = prep_user_table.copy()
    return ext_prep_user_table

def prop_prepare_table(properties_table):

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
        'volumetric_flow_rate',
        'density',
        'specific_heat_capacity',
        'h_fluid_to_pipe',
        'k_pipe',
        'k_insulation',
        'emissivity'
    ]
    return prep_properties_table

def units_handler(properties_table, units):

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

    units_handled_prop_table = prop_prepare_table(properties_table)

    if units.lower() == 'si':
        units_handled_prop_table['inlet_temperature'] = 1.8 * units_handled_prop_table['inlet_temperature'] + 32
        units_handled_prop_table['ambient_temperature'] = 1.8 * units_handled_prop_table['ambient_temperature'] + 32
        units_handled_prop_table['density'] = units_handled_prop_table['density'] * 1 / 16.0185 # 1 lb/ft3 = 16.0185 kg/m3
        units_handled_prop_table['specific_heat_capacity'] = units_handled_prop_table['specific_heat_capacity'] * 1 / 4.1868 # 1 BTU/lb.F = 4.1868 kJ/kg.K
        units_handled_prop_table['h_fluid_to_pipe'] = units_handled_prop_table['h_fluid_to_pipe'] * 1 / 5.67826 # 1 BTU/ft2.h.F = 5.67826 W/m2.K
        units_handled_prop_table['k_pipe'] = units_handled_prop_table['k_pipe'] * 1 / 1.73073 # 1 BTU/ft.h.F = 1.73073 W/m.K
        units_handled_prop_table['k_insulation'] = units_handled_prop_table['k_insulation'] * 1 / 1.73073 # 1 BTU/ft.h.F = 1.73073 W/m.K

    return units_handled_prop_table

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
        'LPM',
        'lbm/ft3',
        'BTU/lbm.F',
        'BTU/ft2.h.F',
        'BTU/ft.h.F',
        'BTU/ft.h.F',
        ''
    ]

    if units.lower() == 'si':
        units_list = [
            'C',
            'C',
            'LPM',
            'kg/m3',
            'kJ/kg.K',
            'W/m2.K',
            'W/m.K',
            'W/m.K',
            ''
        ]

    return units_list

def segment_print(units):

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
