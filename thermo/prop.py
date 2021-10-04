# Previous method of one python module for substance is not desireable, and not very flexible. Instead, a class is made, and is iniatilized with a CoolProp fluid name
# The intention is to pass around the fluid-properties constantly, to where-ever it is needed

from CoolProp.CoolProp import PropsSI, PhaseSI
import physical_constants

class FluidProperties:
    def __init__(self, fluid):
        # CoolProp works with various fluids, but experience has taught me that sometimes, unexpected results can occur.
        # These unexpected results only appeared with mixtures so far, but precautions are taken anyway.
        verified_fluids = ["HEOS::Water", "water" ,"nitrogen"]
        # Documentation says "HEOS::Water" is the most accurate, but it is also slow.
        # "IF97::Water" is faster but slightly less accurate

        if fluid  not in verified_fluids:
            print("WARNING: you are using a fluid that has not been verified.")

        self.fluid = fluid

    def get_enthalpy(self, T, p):
        """Returns SPECIFIC enthalpy for fluider through CoolProp library
        
        Arguments:
            T {K} -- Temperature
            P {Pa} -- Pressure
        
        Returns:
            h {J/kg}
        """
        h = PropsSI('H','T',T,'P',p,self.fluid)
        return h

    def get_density(self, T, p):
        """Returns density of fluid through CoolProp library

        Args:
            T (K): Temperature
            p (Pa): Pressure

        Returns:
            rho (kg/m^3): Density
        """

        rho = PropsSI('DMASS','T', T, 'P', p, self.fluid)
        return rho

    def get_liquid_density_at_psat(self, p_sat):
        """Returns density of fluid as saturated liquid

        Args:
            p_sat (Pa): Saturation Pressure

        Returns:
            rho_sat_l (kg/m^3): Density of saturated liquid
        """

        return PropsSI('DMASS','Q', 0, 'P', p_sat, self.fluid)

    def get_vapour_density_at_psat(self, p_sat):
        """Returns density of fluid as saturated liquid

        Args:
            p_sat (Pa): Saturation Pressure

        Returns:
            rho_sat_g (kg/m^3): Density of saturated gas
        """

        return PropsSI('DMASS','Q', 1, 'P', p_sat, self.fluid)

    def get_pressure(self, T, rho):
        """Returns pressure of fluid through CoolProp library

        Args:
            T (K): Temperature
            rho (kg/m^3): Density

        Returns:
            P {Pa}: Pressure
        """
        return PropsSI('P', 'T', T, 'DMASS', rho, self.fluid)

    def get_viscosity(self, T, p):
        """Returns viscosity of fluid through CoolProp Library

        Args:
            T (K): Temperature
            p (Pa): Pressure

        Returns:
            mu (Pa*s): Viscosity
        """
        viscosity = PropsSI('viscosity', 'T', T, 'P', p, self.fluid)
        return viscosity
    
    def get_liquid_saturation_viscosity(self, p_sat):
        """Returns viscosity at saturation point p_sat in liquid phase

        Args:
            p_sat (Pa): Saturation pressure
        
        Returns:
            mu (Pa*s): Viscosity
        """
        return  PropsSI('viscosity','Q', 0, 'P', p_sat, self.fluid)
    
    def get_gas_saturation_viscosity(self, p_sat):
        """Returns viscosity at saturation point p_sat in gas phase

        Args:
            p_sat (Pa): Saturation pressure
        Returns:
            mu (Pa*s): Viscosity
        """
        return  PropsSI('viscosity','Q', 1, 'P', p_sat, self.fluid)



    def get_specific_heat_ratio(self, T, p):
        """Returns ratio of specific heats (gamma) through CoolProp library

        Args:
            T (K): Temperature
            p (Pa): Pressure

        Returns:
            cp/cv (gamma) (-): Specific heat ratio (gamma)
        """

        cp = PropsSI('CPMASS', 'T', T, 'P', p, self.fluid)
        cv = PropsSI('CVMASS', 'T', T, 'P', p, self.fluid)

        return cp/cv

    def get_cp(self, T, p):
        """Return specific heat under constant pressure

        Args:
            T (K): Temperature
            p (Pa): Pressure

        Returns:
            cp (J/(kg*K)): Specific heat under constant pressure
        """
        return PropsSI('CPMASS', 'T', T, 'P', p, self.fluid)

    def get_liquid_saturation_cp(self, p_sat):
        """Return specific heat under constant pressure at liquid saturation point

        Args:
            p_sat (Pa): Saturation pressure

        Returns:
            cp (J/(kg*K)): Specific heat under constant pressure
        """
        return PropsSI('CPMASS', 'Q', 0, 'P', p_sat, self.fluid)

    def get_gas_saturation_cp(self, p_sat):
        """Return specific heat under constant pressure at gas saturation point

        Args:
            p-sat (Pa): Saturation pressure

        Returns:
            cp (J/(kg*K)): Specific heat under constant pressure
        """
        return PropsSI('CPMASS', 'Q', 1, 'P', p_sat, self.fluid)

    def get_thermal_conductivity(self, T, p):
        """Return thermal conductivity

        Args:
            T (K): Temperature
            p (Pa): Pressure

        Returns:
            kappa (W/(m*K)): Thermal conductivity of fluid
        """
        return PropsSI("CONDUCTIVITY", 'T', T, 'P', p, self.fluid)

    def get_liquid_saturation_conductivity(self, p_sat):
        """Returns conductivity at saturation point p_sat in liquid phase

        Args:
            p_sat (Pa): Saturation pressure

        Returns:
            kappa (W/(m*K)): Thermal conductivity of fluid
        """
        return  PropsSI('CONDUCTIVITY','Q', 0, 'P', p_sat, self.fluid)

    def get_gas_saturation_conductivity(self, p_sat):
        """Returns conductivity at saturation point p_sat in gas phase

        Args:
            p_sat (Pa): Saturation pressure

        Returns:
            kappa (W/(m*K)): Thermal conductivity of fluid
        """
        return  PropsSI('CONDUCTIVITY','Q', 1, 'P', p_sat, self.fluid)

    def get_phase(self, T,p):
        """Returns the phase of fluid through CoolProp library

        Arguments:
            T {K} -- Temperature
            P {Pa} -- Pressure
        
        Returns:
            Returns phase according to CoolProp
            ('liquid'/'gas'/'supercritical_gas'/'supercritical'/'supercritical_liquid')
        """ 

        return PhaseSI('T',T,'P',p,self.fluid)

    def get_specific_gas_constant(self):
        """Returns the specific gas constant for the fluid of choice

        Returns:
            R (J/kg): Specific gas constant (in terms of mass)
        """

        molar_mass = PropsSI("MOLAR_MASS", self.fluid) # Returned in kg/mol
        gas_constant = PropsSI("GAS_CONSTANT", self.fluid) # Returned in J/(mol*K)

        return gas_constant/molar_mass

    # def get_surface_tension(self, T, p):
    #     """Returns surface tension at specified state

    #     Args:
    #         T (K): Temperature
    #         p (Pa): Pressure

    #     Returns:
    #         gamma(N/m): Surface tension
    #     """

    #     return PropsSI("SURFACE_TENSION", "T", T, "P", p, self.fluid)

    def get_surface_tension_at_psat(self, p_sat):
        """Returns surface tension at saturation temperature in liquid state.
        Helper function for get_Bond_number

        Args:
            p_sat (Pa): Saturation Pressure

        Returns:
            gamma(N/m): Surface tension of saturated liquid
        """

        return PropsSI("SURFACE_TENSION", "Q", 0,  "P", p_sat, self.fluid)

    ### FLOW SIMILARITY PARAMETERS BELOW
    # Strictly speaking not part of the thermodynamics, but it is useful to put them here
    
    def get_Reynolds_from_velocity(self, T, p, L_ref, u):
        """Returns the Reynolds number based on given thermodynamic state, chosen reference length and flow velocity

        Args:
            T (K): Temperature
            p (Pa): Pressure
            L_ref (m): Reference length (such as channel length, hydraulic diamater, etc..)
            u (m/s): Flow velocity NOTE: velocity MUST be at same reference state, as density can change from location to location

        Returns:
            Re (-): Reynolds number
        """
        viscosity = self.get_viscosity(T=T, p=p) # [Pa*s] mu - dynamic viscosity
        rho = self.get_density(T=T, p=p) # [kg/m^3] rho - density

        return rho*u*L_ref/viscosity # [-] Reynolds number

    def get_Reynolds_from_mass_flow(self, T, p, L_ref, m_dot, A):
        """Returns the Reynolds number based on given thermodynamic state, chosen reference length and flow velocity

        Args:
            T (K): Temperature
            p (Pa): Pressure
            L_ref (m): Reference length (such as channel length, hydraulic diamater, etc..)
            m_dot (kg/s): Flow velocity NOTE: velocity MUST be at same reference state, as density can change from location to location
            A (m^2): Area through which fluid flows

        Returns:
            Re (-): Reynolds number
        """

        viscosity = self.get_viscosity(T=T, p=p) # [Pa*s]
        rho = self.get_density(T=T, p=p) # [kg/m^3] density
        u = m_dot/(rho*A)

        return rho*u*L_ref/viscosity # [-]


    def get_Prandtl(self, T, p):
        """Return the Prandtly number based on given thermodynamic state

        Args:
            T (K): Temperature
            p (Pa): Pressure

        Returns:
            Pr (-): Prandtl number - ratio of momentum diffusivity/thermal diffusivity
        """
        cp = self.get_cp(T=T, p=p) # [J/kg]
        viscosity = self.get_viscosity(T=T, p=p) # [Pa*s]
        #NOTE: might throw an error if close to saturation point
        kappa = self.get_thermal_conductivity(T=T, p=p) # [W/(m*K)]

        return cp*viscosity/kappa # [-] Prandtl number

    def get_saturation_Prandtl_liquid(self, p_sat):
        """Return the Prandtl number of a saturated liquid at p_sat

        Args:
            p_sat (Pa): Saturation pressure
        """
        cp = self.get_liquid_saturation_cp(p_sat = p_sat) # [J/kg]
        viscosity = self.get_liquid_saturation_viscosity(p_sat=p_sat) # [Pa*s]
        kappa = self.get_liquid_saturation_conductivity(p_sat=p_sat) # [W/(m*K)]

        return cp*viscosity/kappa # [-] Prandtl number

    def get_saturation_Prandtl_gas(self, p_sat):
        """Return the Prandtl number of a saturated gas at p_sat

        Args:
            p_sat (Pa): Saturation pressure
        """
        cp = self.get_gas_saturation_cp(p_sat = p_sat) # [J/kg]
        viscosity = self.get_gas_saturation_viscosity(p_sat=p_sat) # [Pa*s]
        kappa = self.get_gas_saturation_conductivity(p_sat=p_sat) # [W/(m*K)]

        return cp*viscosity/kappa # [-] Prandtl number
        

    def get_Bond_number(self, p_sat, L_ref):
        """ Return the Bond/Eotvos number at the saturation point.
        Most literature seems to give no specific reference state, 
        so this state seems to make the most sense to me, as it concerns boiling fluids

        Args:
            p_sat (Pa): Saturation Pressure
            L_ref (m): Reference Length

        Returns:
            Bo/Eo [-]: Bond/Eötvös number at saturation pressure
        """
        #T_sat = self.get_saturation_temperature(p=p_sat) # [K] Saturation temperature (required to estab)
        surface_tension = self.get_surface_tension_at_psat(p_sat=p_sat) # [N/m]
        # Density in liquid and vapour phase at saturation point
        rho_sat_l = self.get_liquid_density_at_psat(p_sat=p_sat) # [kg/m^3] 
        rho_sat_g = self.get_vapour_density_at_psat(p_sat=p_sat) # [kg/m^3]
        
        return (rho_sat_l-rho_sat_g)*physical_constants.g0*L_ref**2 / surface_tension # [-]

    def get_temperature_from_enthalpy_and_pressure(self, h, p):
        """Return the temperature based on SPECIFIC enthalpy and pressure

        Args:
            h (J/kg): SPECIFIC enthalpy
            p (Pa): Pressure

        Returns:
            T (K): Temperature
        """
        return PropsSI("T","H", h, "P", p, self.fluid)

    def get_saturation_temperature(self, p):
        """Return the saturation temperature of the fluid based on pressure

        Args:
            p (Pa): Saturation pressure

        Returns:
            T_sat (K): Saturation temperature
        """
        # The saturation temperature is found by setting the vapour quality Q somewhere between 0 and 1 (no superheating or subcooling)
        return PropsSI("T" ,"P", p, "Q", 0, self.fluid) # [K] Saturation temperature

    def get_saturation_pressure(self, T):
        """Return the saturation temperature of the fluid based on temperature

        Args:
            T (K): Saturation temperature

        Returns:
            p_sat (Pa): Saturation temperature
        """
        # The saturation pressure is found by setting the vapour quality Q somehwhere between 0 and 1 (no superheating and cooling)
        return PropsSI("P", "T", T, "Q", 0, self.fluid) # [Pa] Saturation pressure

    def get_saturation_enthalpy_gas(self, p):
        """Return the enthalpy of the fluid as a saturated gas

        Args:
            p (Pa): pressure

        Returns:
            h_sat: Saturation enthalpy at the given pressure (gas phase)
        """

        return PropsSI("H", 'P', p, 'Q', 1, self.fluid) # [J/kg] Specific enthalpy

    def get_saturation_enthalpy_liquid(self, p):
        """Return the enthalpy of the fluid as a saturated liquid

        Args:
            p (Pa): pressure

        Returns:
            h_sat: Saturation enthalpy at the given pressure (liquid phase)
        """

        return PropsSI("H", 'P', p, 'Q', 0, self.fluid) # [J/kg] Specific enthalpy

    def get_critical_temperature(self):
        """Return the temperature at the critical point

        Returns:
            T_crit (K): Critical temperature
        """
        return PropsSI("TCRIT", self.fluid) # [K] Critical temperature

    def get_critical_pressure(self):
        """Return the pressure at the critical point

        Returns:
            p_crit (Pa): Critical pressure
        """

        return PropsSI("PCRIT", self.fluid) # [Pa] Critical pressure

    def get_Mach_number(self, T, p, u):
        """Returns the Mach number by dividing provided velocity with speed of sound

        Args:
            T (K): Temperature
            p (Pa): Pressure
            u (m/s): Velocity

        Returns:
            M (-): Mach number
        """

        speed_of_sound = PropsSI("speed_of_sound", 'T', T, 'P', p, self.fluid)

        return u/speed_of_sound # [-] Return the mach number based on given velocity u