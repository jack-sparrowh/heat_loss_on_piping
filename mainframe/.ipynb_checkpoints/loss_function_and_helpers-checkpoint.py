import numpy as np
import pandas as pd

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
