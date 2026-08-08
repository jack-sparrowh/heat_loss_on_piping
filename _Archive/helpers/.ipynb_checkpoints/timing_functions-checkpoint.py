import pandas as pd
import numpy as np
import time
from functools import wraps, partial
from scipy.optimize import minimize

# number of iterations are a function c(x) = 298x + 217
# time follows t(x) = 0.40438+0.41398x+0.001353x^2
# or 0.249+0.457x-0.0000438x^{2}

L = [
    [4, 3.8, 0.5, 876],
    [3, 2.8, 0.5, 2727],
    [2, 1.8, 0.5, 208],
    [2, 1.8, 0.5, 208]
    ]
L = np.array(L)
L = np.c_[L, L[:, -1].cumsum()]
L = pd.DataFrame(L, columns='pipe_outer_diameter pipe_inner_diameter insulation_thickness segment_length cumulative_length'.split(' '))

# prepare
user_table = pd.concat([L.iloc[0, :].to_frame().T, L]).reset_index(drop=True)
user_table.iloc[0, -2:] = 1e-6

properties_table = pd.Series([87, 21, 60, 9.7, 0.025, 190, 60, 1, 0.1], index=['inlet_temperature', 'ambient_temperature', 'h_fluid_to_pipe', 'k_pipe', 'k_insulation', 'volumetric_flow_rate', 'density', 'specific_heat_capacity', 'emissivity'])

c=0
func_times = {}
t = 0
def timer(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        global c, t
        global func_times
        start = time.perf_counter()
        result = func(*args, **kwargs)
        end = time.perf_counter()
        if func.__name__ not in func_times.keys():
            func_times[func.__name__] = np.array([0, end - start])
        else:
            func_times[func.__name__] += np.array([1, end - start])
        c += 1
        t += end - start
        #print(f'{func.__name__} took {end - start:.4f} seconds')
        return result
    return wrapper

# setting it to 1 shifts diameter of the last segment
# lets say we have 2 segments
# 4" 100ft
# 3" 100ft
# if we set haviside to 0 we have:
# from 0 to 100 ft 4" piping
# and from 100+0.00000...1 to 200ft 3"
# setting haviside to 1 shifts diameter to fist segment as 3", which doesnt make any difference as this is just arbitraty which diameter at 100ft we take.
# it does not change results of calculations ,it changes parameters for different segments but the end result is the same 
# it just seems weird for folks w/o mathematical background as it makes a big jump on the segments between diameter change
# so we'll go with 0 here to avoid explaining something that doesn't make any difference
H = lambda x: np.heaviside(x, 0)

def K(
    D_pipe_outer, 
    D_pipe_inner, 
    ins_thickness,
    h_combined_surface,
    h_fluid_to_pipe, 
    k_pipe, 
    k_insulation
    ):

    '''
    Calculate combined thermal conductivity.

    Parameters
    ----------
    D_pipe_outer: float, array_like
        Outer pipe diameter [inch].
    D_pipe_inner: float, array_like
        Inner pipe diameter [inch].
    ins_thickness: float, array_like
        Radius of insulation thickness [inch].
    h_combined_surface: float, array_like
        Combined surface heat transfer coeff (convection + radiation effect) [BTU/ft2.h.F]
    h_fluid_to_pipe: float, array_like
        Heat transfer coeff of heat transfer from the fluid to piping wall [BTU/ft2.h.F]
        If not known, can be set to values above 10 only if the flow is not laminar.
    k_pipe: float, array_like
        Thermal conductivity of piping [BTU/ft.h.F]
    k_insulation: float, array_like
        Thermal conductivity of insulation [BTU/ft.h.F]

    Returns
    -------
    _K: float, array_like
        Combined thermal conductivity [BTU/ft.h.F]
    '''

    _R = (12 * 2) / (h_fluid_to_pipe * D_pipe_inner) +\
         np.log( D_pipe_outer / D_pipe_inner ) / k_pipe +\
         np.log( (D_pipe_outer + 2*ins_thickness) / D_pipe_outer ) / k_insulation +\
         (12 * 2) / (h_combined_surface * (D_pipe_outer + 2*ins_thickness))
    
    _K = 1 / _R
    
    return _K

def heaviside_sum(
    length, 
    h_combined_surface, 
    prep_user_table, 
    mode='integral', 
    func=None, 
    **kwargs
    ):

    '''
    Calculate heaviside summation or integral depending on mode.
    Because calculations can be made for segments of piping of different
    inner and outer diameter as well as insulation thickness and length,
    given heaviside summation approach one can create a continuous function
    of piecewise continuous functions as a "train function".

    Parameters
    ----------
    length: float, array_like
        Length at which properties are calculated [ft]
    h_combined_surface: float, array_like
        Combined surface heat transfer coeff (convection + radiation effect) [BTU/ft2.h.F]
    prep_user_table: pd.DataFrame
        Prepared user table.
    mode: str
        Either "sum" or "integral"
    func: default None
        Function used inside heaviside summation, which must always take 2 variables

    Returns
    -------
    _sum: float, array_like
        Heaviside sum or integral
    '''
    
    
    try:
        # this calculation requires that even if h_combined_surface is of 0th dim
        # it still must be of prep_user_table shape, thus by multiplying by
        # np.ones(dim of prep_user_table) we make sure the calculation continues
        _hv = pd.Series(h_combined_surface * np.ones(prep_user_table.shape[0]))
    except ValueError:
        # if h_combined_surface is not shape 1 or prep_user_table shape
        # then it must still be an ndarray, but of wrong shape
        # to not stop the flow, we set it to first value and let the user know
        _hv = pd.Series(h_combined_surface[0] * np.ones(prep_user_table.shape[0]))
        print("heat transfer coeff vector set to its first value")
        print(f"provided vector size of {h_combined_surface.size} is different from user table of size {prep_user_table.shape[0]}")
    
    if func is None:
        func = lambda table, h_vec: K_partial(
            table.pipe_outer_diameter,
            table.pipe_inner_diameter,
            table.insulation_thickness,
            h_vec,
            **kwargs
        )

    _put = prep_user_table.copy()
    
    # this must be close to 0 but not 0. if it is, then 0 is provided to K which returns NaN and then summation gives NaN.
    _shift_hv = _hv.shift().fillna(1E-9)
    _shift_put = _put.shift().fillna(1E-9) 

    # must be changed to a vector, leaving it as series does not trigger broadcasting and gives ValueError
    _len_diff = length - _shift_put.loc[:, 'cumulative_length'].values.reshape(-1, 1) 
    _func_diff = (func(_put, _hv) - func(_shift_put, _shift_hv)).values.reshape(-1, 1)
    
    _sum = H(_len_diff) * _func_diff
    if mode == 'integral':
        _sum = _sum * _len_diff
    return _sum.sum(0)

def T_l(
    length, 
    h_combined_surface, 
    inlet_temperature, 
    ambient_temperature, 
    density,
    volumetric_flow_rate,
    specific_heat_capacity, 
    **kwargs
    ):

    '''
    Calculate fluid temperature at provided length and combined heat transfer coefficient.

    Parameters
    ----------
    length: float, array_like
        Length at which temperature is to be calculated [ft]
    h_combined_surface: float, array_like
        Combined heat transfer coefficient at the insulation-air interface [BTU/ft2.h.F]
    inlet_temperature: float
        Temperature of fluid at length = 0 [F]
    ambient_temperature: float
        Ambient temperature around piping [F]
    densit: float
        Density of a fluid [lbm/ft3]
    volumetric_flow_rat: float
        Volumetric flow rate of fluid [LPM]
    specific_heat_capacity: float
        Specific heat capacity of fluid [BTU/lbm.R]
    
    Returns
    -------
    T_at_length: float, array_like
        Temperature calculated at provided length [F]
    '''
    
    # 1LPM = 0.0353147 ft3/min * 60 min / 1 h
    exponential_part = np.exp(
        -2 * np.pi / ( density * volumetric_flow_rate * specific_heat_capacity * 0.0353147 * 60 ) *\
        heaviside_sum_partial(length, h_combined_surface, **kwargs)
    )
    
    T_at_length = ambient_temperature + ( inlet_temperature - ambient_temperature ) * exponential_part
    
    return T_at_length

def cp(ref_temperature):
    
    '''
    Calculate specific heat capacity of air at reference temperature, which is film temperature of air around piping.
    Model is taken from Cengel, and is fitted to a range 273K - 1800K. Below 273K cp is equal to cp at 273K and above 
    1800K cp is equal to cp at 1800K.

    Parameters
    ----------
    ref_temperature: float, array_like
        Reference temperature at insulation-air interface [K]
    
    Returns
    -------
    _cp: float, array_like
        Specific heat capacity [kJ/kg.K]
    '''
    
    # all constants here come from "Thermodynamics: An Engineering Approach" by Y. Cengel
    _cp_f = lambda ref_temperature: (28.11 + 0.1967e-2 * ref_temperature + 0.4802e-5 * ref_temperature ** 2 - 1.966e-9 * ref_temperature ** 3) / 28.96
    _cp = _cp_f(ref_temperature)
    _cp[ref_temperature < 273] = _cp_f(273)
    _cp[ref_temperature > 1800] = _cp_f(1800)
    return _cp
    
def u(ref_temperature):
    
    '''
    Calculate kinematic viscosity of air at reference temperature, which is film temperature of air around piping.
    
    Parameters
    ----------
    ref_temperature: float, array_like
        Reference temperature at insulation-air interface [K]
    
    Returns
    -------
    _u: float, array_like
        Kinematic viscosity [Pa*s]
    '''
    
    # Sutherland's Law
    _u = 1.458e-6 * ref_temperature ** 1.5 / (ref_temperature + 110.4)
    return _u

def Ra(
    surface_temperature,
    ambient_temperature,
    length
    ):
        
    '''
    Calculate Rayleigh number of air at film temperature for horizontal cylindrical body. 
    It is assummed that Pr=0.72 for air in all temperature range.

    Parameters
    ----------
    surface_temperature: float, array_like
        Insulation surface temperature [K]
    ambient_temperature: float
        Ambient temperature around piping [K]
    length: float, array_like
        Length to determine outer diameter of piping [ft/m - units_handler takes care of units]
    
    Returns
    -------
    _Ra: float, array_like
        Rayleigh number.
    '''
    
    film_temperature = (surface_temperature + ambient_temperature) / 2
    
    # Pr = 0.72 -> mean value for air on range from 220K to 2200K
    # g = 9.81 -> gravitational acceleration
    # P / R = 101.325 / 0.287 -> atmospheric pressure by gas constant for air
    # 0.0254 -> conversion from inch to meter
    _Ra = 0.72 * 9.81 * (101.325 / 0.287) ** 2 *\
        np.abs(surface_temperature - ambient_temperature) / (film_temperature ** 3 * u(film_temperature) ** 2) *\
        (pipe_diameter_at_length(length) * 0.0254) ** 3
    return _Ra

def h_combined(
    length,
    surface_temperature,
    ambient_temperature,
    emissivity
    ):
    
    '''
    Calculate combined heat transfer coefficient at insulation-air interface.
    Radiation heat transfer coeff. is in its typical form.
    Convective heat transfer is based on Churchill and Chu correlation formula.
    
    Parameters
    ----------
    length: float, array_like
        Length at which combined heat transfer coeff. is to be determined [ft]
    surface_temperature: float, array_like
        Insulation surface temperature [F]
    ambient_temperature: float
        Ambient temperature around piping [F]
    emissivity: float
        Emissivity of insulation material or material covering insulaiton. []
    
    Returns
    -------
    _h_comb_vec: array_like
        Array with heat transfer coefficients of radiative and convective effects.
    '''
    
    surface_temperature_R = surface_temperature + 459.67
    ambient_temperature_R = ambient_temperature + 459.67
    
    # 1.71E-9 -> Stefan-Boltzmann constant
    h_rad = emissivity * 1.71E-9 *\
            (surface_temperature_R ** 2 + ambient_temperature_R ** 2) *\
            (surface_temperature_R + ambient_temperature_R)
    
    surface_temperature_K = surface_temperature_R * 5/9
    ambient_temperature_K = ambient_temperature_R * 5/9
    film_temperature_K = (surface_temperature_K + ambient_temperature_K) / 2
    
    # 0.1762280394 -> conversion from W/m.K to BTU/ft.h.F
    # 1000 -> conversion from kW to W
    # Pr = 0.72 -> mean value for air on range from 220K to 2200K
    # 0.6 + 0.387 / 1.20326 -> constants from Churchill and Chu correlation assuming Pr=0.72
    # 0.0254 -> conversion from inch to meters
    h_conv = 0.1762280394 * 1000 / 0.72 *\
             (0.6 + 0.387 / 1.20326 * Ra(surface_temperature_K, ambient_temperature_K, length) ** (1/6)) ** 2 *\
             u(film_temperature_K) * cp(film_temperature_K) / (pipe_diameter_at_length(length) * 0.0254)
    
    _h_comb_vec = np.c_[h_rad, h_conv]
    
    return _h_comb_vec

def Q_through_pipe(
    length, 
    h_combined_surface, 
    ambient_temperature,
    T_kwargs={},
    H_kwargs={}
    ):
    
    '''
    Calculate total heat transfer rate per unit length through whole piping cross section.

    Parameters
    ----------
    length: float, array_like
        Length at which total heat transfer rate is to be determined [ft]
    h_combined_surface: float, array_like, array_like
        Combined heat transfer coefficient [BTU/ft2.h.F]
    ambient_temperature: float
        Ambient temperature around piping [F]
    
    Returns
    -------
    _Q_through: float, array_like
        Total heat transfer rate per unit length.
    '''
    
    _Q_through = 2 * np.pi * ( T_partial(length, h_combined_surface, **T_kwargs) - ambient_temperature ) *\
            heaviside_sum_partial(length, h_combined_surface, mode='sum', **H_kwargs)
            
    return _Q_through

def Q_to_pipe(
    length,
    h_combined_surface,
    surface_temperature,
    T_kwargs={},
    H_kwargs={}
    ):
    
    '''
    Calculate total heat transfer rate per unit length transferred from the fluid to insulation.

    Parameters
    ----------
    length: float, array_like
        Length at which total heat transfer rate is to be determined [ft]
    h_combined_surface: float, array_like, array_like
        Combined heat transfer coefficient [BTU/ft2.h.F]
    surface_temperature: float, array_like
        Insulation surface temperature [F]
    
    Returns
    -------
    _Q_through: float, array_like
        Total heat transfer rate per unit length from fluid to insulation.
    '''
    
    _Q_to = 2 * np.pi * ( T_partial(length, h_combined_surface, **T_kwargs) - surface_temperature ) *\
            heaviside_sum_partial(length, 10000, mode='sum', **H_kwargs)
            
    return _Q_to
    
def Q_out(
    length, 
    h_combined_surface, 
    surface_temperature,
    ambient_temperature
    ):
    
    '''
    Calculate total heat transfer rate per unit length transferred from the insulation surface to the air surrounding piping.

    Parameters
    ----------
    length: float, array_like
        Length at which total heat transfer rate is to be determined [ft]
    h_combined_surface: float, array_like, array_like
        Combined heat transfer coefficient [BTU/ft2.h.F]
    surface_temperature: float, array_like
        Insulation surface temperature [F]
    ambient_temperature: float
        Ambient temperature around piping [F]
    
    Returns
    -------
    _Q_through: float, array_like
        Total heat transfer rate per unit length from insulation to air.
    '''

    _Q_o = np.pi * pipe_diameter_at_length(length) / 12 * h_combined_surface * ( surface_temperature - ambient_temperature )
    
    return _Q_o

def loss_function(
    x,
    length=None,
    T_kwargs={},
    H_kwargs={}
    ):
    
    '''
    Calculates squared difference between total heat transfer and heat transfer to 
    the outer insulation surface at given surface temperature and length.

    Parameters
    ----------
    x: list, array_like
        Surface temperature [F]
    length: float, array_like
        Length at which total heat transfer rate is to be determined [ft]
    
    Returns
    -------
    loss: float, array_like
        Squared difference.
    '''
    
    surface_temperature = x
    loss = ( 
            Q_to_pipe(length, h_combined_partial(length, surface_temperature).sum(1), surface_temperature, **T_kwargs, **H_kwargs) -\
            Q_through_pipe_partial(length, h_combined_partial(length, surface_temperature).sum(1), **T_kwargs, **H_kwargs) 
    ) ** 2
    return loss

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

# LOADING MAINFRAME FUNCTIONS
K_partial = partial(K, h_fluid_to_pipe=properties_table.h_fluid_to_pipe, k_pipe=properties_table.k_pipe, k_insulation=properties_table.k_insulation)
heaviside_sum_partial = partial(heaviside_sum, prep_user_table=user_table)
T_partial = partial(T_l, inlet_temperature=properties_table.inlet_temperature, ambient_temperature=properties_table.ambient_temperature, density=properties_table.density, volumetric_flow_rate=properties_table.volumetric_flow_rate, specific_heat_capacity=properties_table.specific_heat_capacity)

# LOADING LOSS FUNCTION HELPERS
h_combined_partial = partial(h_combined, ambient_temperature=properties_table.ambient_temperature, emissivity=properties_table.emissivity)
Q_through_pipe_partial = partial(Q_through_pipe, ambient_temperature=properties_table.ambient_temperature)
Q_out_partial = partial(Q_out, ambient_temperature=properties_table.ambient_temperature)

# LOADING MAIN CALCULATION FRAMEWORK
calc_mode_table_partial = partial(calc_mode_table, length_series=user_table.cumulative_length, inlet_temperature=properties_table.inlet_temperature, ambient_temperature=properties_table.ambient_temperature)

# DIAMETER HELPERS
pipe_diameter_at_length = partial(heaviside_sum, h_combined_surface=10000, prep_user_table=user_table, mode='sum', func=lambda table, _: table.pipe_outer_diameter + 2 * table.insulation_thickness)

# adding new to be able to run T_l with ambiguously many h values
pipe_outer_diameter_at_length = partial(heaviside_sum, h_combined_surface=10000, prep_user_table=user_table, mode='sum', func=lambda table, _: table.pipe_outer_diameter)
pipe_inner_diameter_at_length = partial(heaviside_sum, h_combined_surface=10000, prep_user_table=user_table, mode='sum', func=lambda table, _: table.pipe_inner_diameter)
ins_thickness_at_length = partial(heaviside_sum, h_combined_surface=10000, prep_user_table=user_table, mode='sum', func=lambda table, _: table.insulation_thickness)