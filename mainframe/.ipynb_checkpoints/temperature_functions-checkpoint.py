import numpy as np
import pandas as pd

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
