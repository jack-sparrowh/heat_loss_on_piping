import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.integrate import quad

class Air:

    '''
    Describes air properties around piping. 
    Provides all necessary functions for heat transfer coefficient determination.
    
    Parameters
    ----------
    outer_diameter: float, array-like
    	Piping outer diameter [inch]
    insulation_thickness: float, array-like
    	Insulation thickness [inch]
    ambient_temperature: float, array-like
    	Ambient temperature [K]
    emissivity: float
    	Piping surface emissivity []
    
    Returns
    -------
    Air_class: 
    	Air class.
    '''
    
    def __init__(
        self,
        outer_diameter,
        insulation_thickness,
        ambient_temperature,
        emissivity
    ):

        self.total_diameter = outer_diameter + 2 * insulation_thickness
        self.ambient_temperature = ambient_temperature
        self.emissivity = emissivity
        
    
    def cp(
        self,
        ref_temperature
    ):

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
        
        # going with np.array below assures that even when ref_temp is int all calculations carry on
        _cp = np.array(_cp_f(ref_temperature))
        _cp[ref_temperature < 273] = _cp_f(273)
        _cp[ref_temperature > 1800] = _cp_f(1800)
        return _cp

    def u(
        self,
        ref_temperature
    ):
    
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
        self,
        ambient_temperature,
        surface_temperature,
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
    
        Returns
        -------
        _Ra: float, array_like
            Rayleigh number.
        '''

        film_temperature = (surface_temperature + ambient_temperature) / 2
    
        # Pr = 0.72 -> mean value for air on range from 220K to 2200K
        # g = 9.81 -> gravitational acceleration
        # P / R = 101.325 / 0.287 -> atmospheric pressure over gas constant for air
        # 0.0254 -> conversion from inch to meter
        _Ra = 0.72 * 9.81 * (101.325 / 0.287) ** 2 *\
            np.abs(surface_temperature - ambient_temperature) / (film_temperature ** 3 * self.u(film_temperature) ** 2) *\
            (self.total_diameter * 0.0254) ** 3
        
        return _Ra
    
    def h_combined(
        self,
        surface_temperature
    ):
    
        '''
        Calculate combined heat transfer coefficient at insulation-air interface.
        Radiation heat transfer coeff. is in its typical form.
        Convective heat transfer is based on Churchill and Chu correlation formula.
    
        Parameters
        ----------
        surface_temperature: float, array_like
            Insulation surface temperature [F]
    
        Returns
        -------
        _h_comb_vec: array_like
            Array with heat transfer coefficients of radiative and convective effects.
        '''
    
        surface_temperature_R = surface_temperature + 459.67
        ambient_temperature_R = self.ambient_temperature + 459.67
    
        # 1.71E-9 -> Stefan-Boltzmann constant
        h_rad = self.emissivity * 1.71E-9 *\
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
                 (0.6 + 0.387 / 1.20326 * self.Ra(surface_temperature_K, ambient_temperature_K) ** (1/6)) ** 2 *\
                 self.u(film_temperature_K) * self.cp(film_temperature_K) / (self.total_diameter * 0.0254)
    
        _h_comb_vec = np.c_[h_rad, h_conv]
    
        return _h_comb_vec

class Piping_Segment:

    '''
    Describes one piping section defined by user.

    Parameters
    ----------
    outer_diameter: float
        Outer diameter [inch]
    inner_diameter: float
        Inner diameter [inch]
    insulation_thickness: float
        Insulation thickness [inch]
    segment_length: float
        Length of pipng segment [ft]
    start_length: float
        Initial cumulative length. If there is more than one segment, then start_length_i = end_segment_(i-1) [F]
    inlet_temperature: float
        Temperature at start length. If there is more than one segment, then inlet_temperature_i = final_temperature_(i-1) [F]
    ambient_temperature: float
        Ambient temperature [F]
    h_fluid_to_pipe: float
        Heat transfer coeff from fluid to wall of piping [BTU/ft2.h.F]
    k_pipe: float
        Thermal conductivity of piping material [BTU/ft.h.F]
    k_insulation: float
        Thermal conductivity of insulation material [BTU/ft.h.F]
    mass_flow_rate: float
        Mass flow rate [lbm/h]
    specific_heat_capacity: float
        Specific heat capacity of liquid [BTU/lbm.R]
    h_combined_func: function
        Function for calculating combined heat transfer coeff. Must take one argument only []
    
    Returns
    -------
    Piping_Segment_class: 
        Piping segment class.
    '''
    
    def __init__(
        self,
        outer_diameter,
        inner_diameter,
        insulation_thickness,
        segment_length,
        start_length,
        inlet_temperature,
        ambient_temperature,
        h_fluid_to_pipe, 
        k_pipe, 
        k_insulation,
        mass_flow_rate,
        specific_heat_capacity,
        solar_flux,
        h_combined_func
    ):

        # piping geometry
        self.outer_diameter = outer_diameter
        self.inner_diameter = inner_diameter
        self.insulation_thickness = insulation_thickness
        self.segment_length = segment_length
        self.start_length = start_length
        self.total_diameter = self.outer_diameter + 2 * self.insulation_thickness
        self.end_length = self.start_length + self.segment_length
        
        # piping thermal properties
        self.h_fluid_to_pipe = h_fluid_to_pipe
        self.k_pipe = k_pipe
        self.k_insulation = k_insulation

        # fluid properties
        self.inlet_temperature = inlet_temperature
        self.ambient_temperature = ambient_temperature
        self.mass_flow_rate = mass_flow_rate
        self.specific_heat_capacity = specific_heat_capacity

        # add heat rate per segment length, that will be populated during solution for total piping
        self.heat_rate_psl_total = 0
        self.heat_rate_psl_surface = 0 

        # heat transfer coefficient functoin
        self.h_comb_vec = lambda T_surface: h_combined_func(T_surface)
        self.h_comb_func = lambda T_surface: h_combined_func(T_surface).sum(1)

        # define derivative function
        self.derivative = lambda f, x: ( f(x+0.00001) - f(x-0.00001) ) / 0.00002

        # find inlet surface temperature
        self.find_inlet_surface_temperature()
        
        # set normal bounds for find_surface_temperature
        # we use sorted() for when inlet temperature is lower than ambient temperature
        self.bounds=sorted([self.ambient_temperature, self.inlet_surface_temperature])

        # solar flux
        self.solar_flux = solar_flux
        if solar_flux > 0:
            self.solar_radiation_helper()

    def K(
        self,
        T_surface
    ):
    
        '''
        Calculate combined thermal conductivity.
    
        Parameters
        ----------
        T_surface: float, array_like
            Insulation surface temperature [F]
    
        Returns
        -------
        _K: float, array_like
            Combined thermal conductivity [BTU/ft.h.F]
        '''

        _R = (12 * 2) / (self.h_fluid_to_pipe * self.inner_diameter) +\
             np.log( self.outer_diameter / self.inner_diameter ) / self.k_pipe +\
             np.log( self.total_diameter / self.outer_diameter ) / self.k_insulation +\
             (12 * 2) / (self.h_comb_func(T_surface) * self.total_diameter)
    
        self.K_combined = 1 / _R

        return self.K_combined


    def ratio_of_resistances(
        self,
        T_surface
    ):

        '''
        Calculate ratio of thermal resistances. The returned ratio finds R_fa/R_sa where
        R_fa = 1/(2*pi*K)
        R_sa = 24/(2*pi*D*h).

        Parameters
        ----------
        T_surface: float
        	Insulation surface temperature [F]
        
        Returns
        -------
        _res_ratio: float
        	Calculated resistance ratio.
        '''

        _res_ratio = self.total_diameter * self.h_comb_func(T_surface) / (24 * self.K(T_surface))

        return _res_ratio
    
    def find_fluid_temperature(
        self, 
        T_surface
    ):
    
        '''
        Calculate fluid temperature at provided surface temperature.
    
        Parameters
        ----------
        T_surface: float, array_like
            Insulation surface temperature [F]
    
        Returns
        -------
        _T_fluid: float, array_like
            Temperature calculated at provided surface temperature [F]
        '''

        _T_fluid = self.ambient_temperature + (T_surface - self.ambient_temperature) *\
                    self.ratio_of_resistances(T_surface)
    
        self.fluid_temperature = _T_fluid

        return self.fluid_temperature

    def find_inlet_surface_temperature(
        self
    ):
        '''
        Calculate inlet surface temperature.

        Function minimizez squared difference of:
        R_fa/R_sa*(Ts-Ta)-(Tf-Ta)

        Parameters
        ----------
        
        Returns
        -------
        self.inlet_surface_temperature: float
        	Inlet surface temperature.
        '''

        _T_surf_loss_func = lambda T_surface: (
            (T_surface - self.ambient_temperature) * self.ratio_of_resistances(T_surface) -\
            (self.inlet_temperature - self.ambient_temperature)
        ) ** 2

        # we use sorted() for when inlet temperature is lower than ambient temperature
        self.minim_inlet_surface_temperature = minimize_scalar(
            _T_surf_loss_func,
            bounds=sorted([self.ambient_temperature, self.inlet_temperature])
        )

        # make sure to extract float, so scipy can accept it
        self.inlet_surface_temperature = self.minim_inlet_surface_temperature.x[0]

    def solar_radiation_helper(
        self,
        threshold=0.01
    ):
        
        # calculate constant temperature solar flux
        self.const_temp_sol_flux = self.h_comb_func(self.inlet_surface_temperature) \
                                    * (self.inlet_surface_temperature - self.ambient_temperature)

        # find how close provided solar flux is to constant temp solar flux
        _closeness = self.solar_flux / self.const_temp_sol_flux - 1

        _loss = lambda T_s: ( self.h_comb_func(T_s) * (T_s - self.ambient_temperature) - self.solar_flux ) ** 2

        if _closeness >= 0:
            
            upper_bound = minimize_scalar(_loss).x[0]
            self.bounds = sorted([self.ambient_temperature, upper_bound])

            if _closeness <= threshold:
                self.solar_flux = (1 + threshold) * self.const_temp_sol_flux

        else:

            lower_bound = minimize_scalar(_loss).x[0]
            self.bounds = sorted([lower_bound, self.inlet_surface_temperature])

            if _closeness >= -threshold:
                self.solar_flux = (1 - threshold) * self.const_temp_sol_flux

    def integrand(
        self,
        T_surface
    ):
        '''
        Find value of integrand, given surface temperature.

        Parameters
        ----------
        T_surface: float
        	Insulation surface temperature [F]
        
        Returns
        -------
        _integrand: float
        	Calculated integrand.
        '''

        ln_res_ratio = lambda T_surface: np.log(self.ratio_of_resistances(T_surface))

        _I = 1 / self.K(T_surface) * ( 1 / (T_surface - self.ambient_temperature) + self.derivative(ln_res_ratio, T_surface) )

        _integrand = _I / ( 1 - self.solar_flux / (self.h_comb_func(T_surface) * (T_surface - self.ambient_temperature)) )

        # make sure to extract float, so scipy can accept it
        return _integrand[0]

    def integral_loss_function(
        self,
        T_surface,
        length
    ):

        '''
        Calculates squared sum of integral over surface temperature and integral over distance.

        Parameters
        ----------
        T_surface: float
        	Insulation surface temperature [F]
        length: float
        	Distance at which to find surface temperature [ft]
        
        Returns
        -------
        loss_function: float
        	Calculated loss.
        '''
        
        _lambda = self.mass_flow_rate * self.specific_heat_capacity / (2 * np.pi)
        integral = quad(self.integrand, self.inlet_surface_temperature, T_surface, limit=500)
        loss = (integral[0] + length / _lambda) ** 2

        return loss

    def find_surface_temperature(
        self,
        length
    ):

        '''
        Finds surface temperature, at specified length by minimizing integral_loss_funciton.

        Parameters
        ----------
        length: float
        	Distance at which to find surface temperature [ft]
        
        Returns
        -------
        self.surface_temperature: float
        	Surface temperature.
        '''

        self.minim_surface_temperature = minimize_scalar(
            self.integral_loss_function,
            bounds=self.bounds,
            args=(length)
        )

        self.surface_temperature = self.minim_surface_temperature.x
        
        return self.surface_temperature
        
    def calculate_heat_per_length(
        self, 
        k, 
        LMTD
    ):

        '''
        Calculates heat transfer rate per unit length of piping.
        The form of equation might suggests that only conductivity can be passed as an argument,
        but one can utilize thermal resistance network approach to connect conductivity to any
        other reqiured value. As an example, user can pass:
        (Ts-Ta) / (Tf-Ta) * D_t/2 * h
        to find how heat transfer changes with changing heat transfer coeff. on surface of piping.
        
        Parameters
        ----------
        k: float, array-like
        	Pseudo-conductivity. Must have a unit of conductivity [BTU/ft.h.F]
        LMTD: float, array-like
        	Log-mean temperature difference [F]
        
        Returns
        -------
        _hpl: float, array-like
        	Heat transfer per length.
        '''

        _hpl = 2 * np.pi * k * LMTD
        
        return _hpl
        
    def temperature_profile(self, n=2):

        '''
        Find fluid temperature over a whole segment.
        Parameter n controls to how many parts a given segment is to be split to.
        Calculations are carried out by filling separate lists, and then they are combined
        to one pandas dataframe. At the end, this dataframe is saved into class as "temp_prof".
        
        Parameters
        ----------
        n: int
        	Number of parts to split the given segment to. Default n=2 []
        
        Returns
        -------
        self.temp_prof: pd.DataFrame
        	Temperature profile.
        '''

        # initialize at x=0
        temp_profile = [self.inlet_temperature]
        surf_temp_profile = [self.inlet_surface_temperature]
        h_comb_profile = [self.h_comb_vec(self.inlet_surface_temperature).flatten()]
        segment_length = [0]

        # run through rest n-1 segments
        for _l in np.linspace(self.segment_length / (n - 1), self.segment_length, n - 1):
 
            _T_surf_l = self.find_surface_temperature(_l)
            _h_comb_vec_l = self.h_comb_vec(_T_surf_l).flatten() 
            _T_l = self.find_fluid_temperature(_T_surf_l)

            temp_profile.append(_T_l[0])
            surf_temp_profile.append(_T_surf_l)
            h_comb_profile.append(_h_comb_vec_l)
            segment_length.append(_l)

        _result_df = pd.DataFrame(
            np.c_[segment_length, temp_profile, surf_temp_profile, h_comb_profile],
            columns=["length", "temp", "surf_temp", "h_rad", "h_conv"]
        )

        self.temp_prof = _result_df

        return self.temp_prof