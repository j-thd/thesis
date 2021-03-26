import numpy as np
from physical_constants import stefan_boltzmann, g0
from thermo.prop import FluidProperties

def ideal_enthalpy_change(T_inlet, p_inlet, T_outlet, p_outlet, fp: FluidProperties):
    """Returns specific enthalpy change based on simple chamber inlet and outlet conditions.
    This should give the power the micro-heater must transfer in ideal conditions with no heat losses.
    In addition returns a warning if the final state is not gaseous.
    
    Arguments:
        T_inlet {K} -- Inlet temperature
        p_inlet {Pa} -- Inlet pressure
        T_outlet {K} -- Outlet temperature
        p_outlet {Pa} -- Outlet pressure
        fp {object} -- FluidProperties object
    
    Returns:
        delta_h {J/(kg*K)} -- 
    """

    h_inlet = fp.get_enthalpy(T=T_inlet, p=p_inlet)
    h_outlet = fp.get_enthalpy(T=T_outlet, p=p_outlet)

    outlet_phase = fp.get_phase(T=T_outlet,p=p_outlet)
    if not (outlet_phase == 'gas' or outlet_phase == 'supercritical_gas'):
        print("Warning: Phase at chamber exit is not gaseous but {}".format(outlet_phase))
        
    return h_outlet-h_inlet

def ideal_power_consumption(mass_flow, T_inlet, p_inlet, T_outlet, p_outlet, fp: FluidProperties):
    """Returns power consumption assuming no heat losses, based on enthalpy change and mass flow alone
    
    Arguments:
        mass_flow {kg/s} -- Mass flow
        T_inlet {K} -- Inlet temperature
        p_inlet {Pa} -- Inlet pressure
        T_outlet {K} -- Outlet temperature
        p_outlet {Pa} -- Outlet pressure
    
    Returns:
        {W} -- Required micro-heater power
    """

    return mass_flow * ideal_enthalpy_change(T_inlet=T_inlet, p_inlet=p_inlet, T_outlet=T_outlet, p_outlet=p_outlet, fp=fp)

def ideal_heater_temperature(P_mh, T_inlet, T_outlet, A, k, d):
    """Return the ideal heater temperature, assuming that all heat is conducted from distance d towards the chamber, with an uniform temperature distribution and no further heat losses

    
    Arguments:
        P_mh {W} -- Required micro-heater power
        T_inlet {K} -- Inlet temperature of chamber
        T_outlet {K} -- Outlet temperature of chamber
        A {m^2} -- Surface area of chamber wall (assumed to be same as heater surface area)
        k {W/(m*K)} -- Thermal conductivity of chip
        d {m} -- Thickness of chip i.e: distance from chamber wall to heater
    
    Returns:
        T_mh {K -- Temperature needed to provide required power to chamber
    """
    T_ref = (T_inlet + T_outlet) / 2

    return P_mh * d / (A*k) + T_ref

def radiation_loss(T, A, emmisivity):
    """Return radiation loss based on black-body radiation 
    
    Arguments:
        T {K} -- Temperature of radiator
        A {m^2} -- Area from which it is radiated
        emmisivity {[type]} -- [description]
    
    Returns:
        [type] -- [description]
    """
    #print(type(T))
    T = float(T)
    #print(type(T))
    #print("Temperature**4: {} ".format(T**4))
    #print(type(A))
    A = float(A)
    print("Area: {}".format(A))
    #print(type(A))
    #print(type(emmisivity))
    emmisivity = float(emmisivity)
    #print(type(emmisivity))
    #print(type(stefan_boltzmann))
    P_rad = emmisivity * A * stefan_boltzmann * T**4
    print("P_rad: {}".format(P_rad))
    #print(type(P_rad))
    return emmisivity * A * stefan_boltzmann * T**4

def convective_heat_flow(heat_transfer_coefficient, T_wall, T_ref, A_wall):
    """Return the heat flow due to any form of heat transfer (such as convection or conduction)
    heat_transfer_coefficient is written out fully to avoid confusion with specific enthalpy h

    NOTE: the heat transfer coefficient must be positive, and heat flowing away from T_wall when T_wall>T_ref is considered to be NEGATIVE heat flow
    Args:
        heat_transfer_coefficient (W/(m^2*K)): h -heat transfer coefficient due to arbitrary form of transfer such as conduction (kappa) or convection (h_conv)
        T_wall (K): Temperature at wall
        T_ref (K): Reference temperature (NOTE: it must be known what the reference temperature was for determining h_conv, as it is an emperical relation)
        A_wall (m^2): Area through which heat is transferred

    Returns:
        Q_dot (W): Heat flow through area A_wall 
    """

    assert(heat_transfer_coefficient>0) # Return error when heat transfer coefficient is not positive
    return -heat_transfer_coefficient*(T_wall-T_ref)*A_wall # [W] Q_dot: heat flow

def T_wall_from_heat_flow(Q_dot, heat_transfer_coefficient, T_ref, A_wall):
    """ Return the wall temperature from known heat flow, h_conv, reference temperature and chamber wall dimension

    Args:
        Q_dot (W): Required heat flux from wall to fluid (is negative if heat flows from wall to fluid)
        heat_transfer_coefficient (W/(m^2*K)): Heat transfer coefficient from arbitrary from of heat transfer, such as conduction (kappa) or convection (h_conv)
        T_ref (K): Reference temperature( NOTE: it must be known what the reference temperature was for determing the coefficent h_conv as it is an emperical relation)
        A_wall (m^2): Area through which the heat is transfffer

    Returns:
        T_wall (K): The temperature of the wall that would result in the heat transfer Q_dot
    """
    assert(heat_transfer_coefficient>0) # Return error when heat transfer coefficient is not positive. Q_dot must be negative when T_wall is T_ref after all
    assert(Q_dot < 0) # Short warning to check if Q_dot is negative, which is what is desired, if one wants to heat up the flow

    return -Q_dot/(heat_transfer_coefficient*A_wall) + T_ref # [K] Returns required wall temperature (should be inverse of function above)

def h_conv_from_Stanton(Stanton, u, T_ref, p_ref, fp: FluidProperties):
    """Return heat transfer coefficient dependent on Stanton number and thermodynamic state\
        Temperature and pressure are passed instead of T and p, as these abstract away constant computations of cp and rho\
            and it is easier to pass around the same state variables time and time again

            WARNING: REFERENCE THERMODYNAMIC STATE (T_ref, p_ref) MUST BE EQUAL TO THOSE WITH WHICH NUSSELT NUMBER WAS DETERMINED

    Args:
        Stanton (-): Stanton number: dimensionless flow characteristic
        u (m/s): flow velocity
        T_ref (K): Temperature
        p_ref (Pa): Pressure
        fp (FluidProperties): object to use to obtain properties of fluid

    Returns:
        h_conv (W/(m^2*K)): convective heat transfer coefficient based on Stanton number, flow velocity and thermodynamic state
    """
    cp = fp.get_cp(T=T_ref, p=p_ref) # [J/kg] Specific heat capacity under constant pressure
    rho = fp.get_density(T=T_ref,p=p_ref) # [kg/m^3] Fluid density

    return Stanton*rho*u*cp # [W/(m^2*K)] h_conv: Convective heat transfer coefficient

def wetter_perimeter_rectangular(w_channel, h_channel):
    """Calculate the wetted perimeter (for hydraulic diameters) of a rectangular channel

    Args:
        w_channel (m): channel width
        h_channel (m): chann height (or depth)

    Returns:
        P (m): Wetter perimeter of a rectangular channel
    """
    return 2 * (w_channel + h_channel)

def hydraulic_diameter_rectangular(w_channel, h_channel):
    """ Calculate the hydraulic diameter for a rectangular channel

    Args:
        w_channel (m): Channel width
        h_channel (m): Channel height/depth

    Returns:
        D_h (m): Hydraulic diameter of a rectangular channel
    """
    A = w_channel*h_channel # [m^2] Channel cross-sectional area`
    P = wetter_perimeter_rectangular(w_channel=w_channel, h_channel=h_channel) # [m] Wetted perimeter of channel

    return 4 * A / P # [m] 

def velocity_from_mass_flow(A, m_dot, rho):
    """ Calculate channel velocity from mass flow

    Args:
        A (m^2): Channel area
        m_dot (kg/s): Mass flow
        rho (kg/m^3): Density

    Returns:
        u (m/s): Velocity inside a channel
    """


    return m_dot / (rho* A)
    
def mass_flow(A, u, rho):
    """ Returns the mass flows based on density, flow velocity and area

    Args:
        A (m^2): Channel area through which the fluid flow
        u (m/s): Flow velocity
        rho (kg/m^3): Density of the fluid

    Returns:
        m_dot (kg/s): Mass flow through the chamber channels
    """
    return rho*A*u # [kg/s]


def required_power(m_dot, delta_h):
    """Required power to raise enthalpy of flow with mass flow m_dot with delta_h

    Args:
        m_dot (kg/s): Mass flow
        delta_h (J/kg): Specific enthalpy increase
    """

    return m_dot * delta_h # [W] Required power (Power is positive if heat flows from T_wall to T_ref)

def required_heater_area(Q_dot, h_conv, T_wall, T_ref):

    #assert(Q_dot > 0) # Check if Q_dot is positive as per convention
    #assert(T_wall > T_ref)

    return Q_dot/(h_conv * (T_wall - T_ref)) # [m^2] Required heater area 

def required_chip_area(L_channel, w_channel, w_channel_margin):
    return L_channel * (w_channel + 2*w_channel_margin)

def Reynolds_from_mass_flow(m_dot, A_channel, L_ref, mu):
    """Calculate Reynolds number from mass flow

    Args:
        m_dot (kg/s): Mass flow
        A_channel (m^2): Cross-sectional area of channel
        L_ref (m): Reference length (probably hydraulic diameter)
        mu (Pa*s): Viscosity

    Returns:
        Re (-): Reynold's number
    """
    return m_dot*L_ref/ (A_channel*mu)



def Froude_number(u, L_ref):
    """ Calculate Froude number from velocity and reference length

    Args:
        u (m/s): Flow velocity
        L_ref (m): Reference length (like hydraulic diameter)

    Returns:
        Fr (-): Froude number
    """

    return u / (g0 * L_ref)**0.5 # [-] Froude number

def delta_enthalpy_per_section(h):
    """Calculate the enthalpy difference per section, to facilite calculation of how much power must be added in a section. The first section defaults to zero.

    Args:
        h (J/kg): Enthalpy at each station

    Returns:
        delta_h: Enthalpy difference between section i and i-1 (i=0 defaults to zero)
    """
    # Heating power required in section to increase temp by dT or quality by dx. Use enthalpy difference
    # Make the first element zero to make the array conventiently precisely as long as T. Arbitrary choice as to whether last or first element is zero
    delta_h = np.hstack((0, np.diff(h))) # [J/kg]
    return delta_h