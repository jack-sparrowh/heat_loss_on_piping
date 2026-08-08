import numpy as np
import pandas as pd
from scipy.optimize import minimize

def calc_mode_table(
    mode, 
    length_series, 
    inlet_temperature, 
    ambient_temperature, 
    T_kwargs={},
    H_kwargs={}
    ):

    '''
    Prepare a semi-final dataframe containing values used by final function to create user-redable table.
    Creates a pandas dataframe with cumulative length, surface temperature minimizing loss function,
    heat transfer coefficients for both radiative heat transfer and convective and calculated
    temperature of fluid at specified length.

    Parameters
    ----------
    mode: str
        Mode of calculation. Either "conservative" or "exact" []
    length_series: array_like
        Series or array of cummulative lengths []
    inlet_temperature: float
        Temperature of fluid at length = 0 [F]
    ambient_temperature: float
        Ambient temperature around piping [F]

    Returns
    -------
    _semi_final_table: pd.DataFrame
        Semi-final table.
    '''

    if mode.lower() == 'conservative':
        # calculate assuming h_comb -> inft
        _surf_temp = np.nan * np.ones(length_series.size)
        _h_comb = np.nan * np.ones((length_series.size, 2))
        _T = T_partial(length_series.values, 10000, **T_kwargs)

    elif mode.lower() == 'exact':
        # calculate by minimizing loss at each length given by user
        _surf_temp_at_length_minim = np.array([[
            _length, 
            minimize(
                loss_function, 
                [( inlet_temperature + ambient_temperature ) / 2], 
                args=(_length)
            ).x[0]
        ] for _length in length_series])
        _h_comb = h_combined_partial(_surf_temp_at_length_minim[:, 0], _surf_temp_at_length_minim[:, 1])
        _surf_temp = _surf_temp_at_length_minim[:, 1]
        _T = T_partial(length_series.values, _h_comb.sum(1), **T_kwargs)

    _semi_final_table = pd.DataFrame(
        np.c_[length_series, _surf_temp, _h_comb[:, 0], _h_comb[:, 1], _T],
        columns=['cumulative_length', 'surface_temperature_minim', 'h_rad', 'h_conv', 'temperature_at_length']
    )

    return _semi_final_table

def final_table(
    mode,
    **kwargs
    ):

    '''
    Creates final table before units conversion. It calculates all selected heats to show the difference 
    between heat to the insulation and from the insulation as a sum of radiation and convection. 
    It also attaches simple to folow column names.

    Parameters
    ----------
    mode: str
        Mode of calculation. Either "conservative" or "exact" []

    Returns
    -------
    _final_table: pd.DataFrame
        Final table ready for unit conversion.
    '''

    _semi_table = calc_mode_table_partial(mode, **kwargs)

    _semi_table.loc[:, 'segment_length'] = _semi_table['cumulative_length'].diff().fillna(0)

    if mode.lower() == 'conservative':
        _semi_table.loc[:, 'heat_to_ins'] = np.array([
            Q_through_pipe_partial(
                _semi_table.loc[i, 'cumulative_length'],
                10000
            ) for i in range(_semi_table.shape[0])
        ])
        _semi_table.loc[:, 'surface_temperature_minim'] = properties_table.ambient_temperature

    elif mode.lower() == 'exact':
        _semi_table.loc[:, 'heat_to_ins'] = np.array([
            Q_through_pipe_partial(
                _semi_table.loc[i, 'cumulative_length'],
                _semi_table.loc[i, ['h_rad', 'h_conv']].sum()
            ) for i in range(_semi_table.shape[0])
        ])

    _semi_table.loc[:, 'heat_by_convection'] = np.array([
        Q_out_partial(
            _semi_table.loc[i, 'cumulative_length'],
            _semi_table.loc[i, 'h_conv'],
            _semi_table.loc[i, 'surface_temperature_minim']
        ) for i in range(_semi_table.shape[0])
    ])

    _semi_table.loc[:, 'heat_by_radiation'] = np.array([
        Q_out_partial(
            _semi_table.loc[i, 'cumulative_length'],
            _semi_table.loc[i, 'h_rad'], 
            _semi_table.loc[i, 'surface_temperature_minim']
        ) for i in range(_semi_table.shape[0])
    ])

    _semi_table.loc[:, 'heat_sum'] = _semi_table.loc[:, 'heat_by_convection'] + _semi_table.loc[:, 'heat_by_radiation']

    _semi_table.loc[:, 'diff'] = np.abs(_semi_table.loc[:, 'heat_to_ins'] - _semi_table.loc[:, 'heat_sum'])

    columns_to_display = [
        'segment_length',
        'cumulative_length',
        'temperature_at_length',
        'surface_temperature_minim',
        'heat_to_ins',
        'heat_by_radiation',
        'heat_by_convection',
        'heat_sum',
        'diff'
    ]

    _final_table = _semi_table.loc[:, columns_to_display]

    return _final_table

def final_table_units_handler(mode, units):

    '''
    Converts final table to correct units.

    Parameters
    ----------
    mode: str
        Mode of calculation. Either "conservative" or "exact" []
    units: str
        Either "SI" or "Imperial" []

    Returns
    -------
    _uh_final_table: pd.DataFrame
        Final table.
    '''

    _uh_final_table = final_table(mode)

    rename_columns = lambda length, temp, heat: [
        f'Segment Length [{length}]',
        f'Cumulative Length [{length}]',
        f'Temperature [{temp}]',
        f'Surface Temperature [{temp}]',
        f'Heat Transfer to Insulation [{heat}]',
        f'Heat Transfer by Radiation [{heat}]',
        f'Heat Transfer by Convection [{heat}]',
        f'Heat Transfer Total [{heat}]',
        f'Difference [{heat}]'
    ]

    length, temp, heat = ['ft', 'F', 'BTU/h.ft']

    if units.lower() == 'si':
        length, temp, heat = ['m', 'C', 'W/m']
        _uh_final_table['segment_length'] = _uh_final_table['segment_length'] * 0.3048 # 1 ft = 0.3048 m
        _uh_final_table['cumulative_length'] = _uh_final_table['cumulative_length'] * 0.3048 # 1 ft = 0.3048 m
        _uh_final_table['temperature_at_length'] = (_uh_final_table['temperature_at_length'] - 32) / 1.8
        _uh_final_table['surface_temperature_minim'] = (_uh_final_table['surface_temperature_minim'] - 32) / 1.8
        _uh_final_table['heat_to_ins'] = _uh_final_table['heat_to_ins'] * 0.293071 / 0.3048 # 1 BTU/h.ft = 0.293071 / 0.3048 W/m
        _uh_final_table['heat_by_radiation'] = _uh_final_table['heat_by_radiation'] * 0.293071 / 0.3048 # 1 BTU/h.ft = 0.293071 / 0.3048 W/m
        _uh_final_table['heat_by_convection'] = _uh_final_table['heat_by_convection'] * 0.293071 / 0.3048 # 1 BTU/h.ft = 0.293071 / 0.3048 W/m
        _uh_final_table['heat_sum'] = _uh_final_table['heat_sum'] * 0.293071 / 0.3048 # 1 BTU/h.ft = 0.293071 / 0.3048 W/m
        _uh_final_table['diff'] = _uh_final_table['diff'] * 0.293071 / 0.3048 # 1 BTU/h.ft = 0.293071 / 0.3048 W/m

    _uh_final_table.columns = rename_columns(length, temp, heat)

    return _uh_final_table.round(3)
