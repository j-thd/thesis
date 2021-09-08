import math
import numpy as np
import scipy.optimize
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
    #A = float(A)
    #print("Area: {}".format(A))
    #print(type(A))
    #print(type(emmisivity))
    emmisivity = float(emmisivity)
    #print(type(emmisivity))
    #print(type(stefan_boltzmann))
    #P_rad = emmisivity * A * stefan_boltzmann * T**4
    #print("P_rad: {}".format(P_rad))
    #print(type(P_rad))
    return emmisivity * A * stefan_boltzmann * T**4

def radiation_loss_numpy(T_np_array, A, emmisivity):
    # T must become iterable to convert each temperature to a Python float (to prevent overflow in numpy float objects when T**4)
    T_iter = np.nditer(T_np_array, flags=['c_index']) # [K] 
    P_rad = np.zeros_like(T_np_array)

    for T_i in T_iter:
        T = float(T_i) # [K] Convert temperature to Python Float
        P_rad[T_iter.index] = emmisivity * A * stefan_boltzmann * T**4 # [W]

    return P_rad
    

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

def wetted_perimeter_rectangular(w_channel, h_channel):
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
    P = wetted_perimeter_rectangular(w_channel=w_channel, h_channel=h_channel) # [m] Wetted perimeter of channel

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

def inlet_manifold_width(channel_amount, w_channel, w_channel_spacing, inlet_manifold_width_factor):
    """Calculate total inlet manifold width, with margins, based on preset factors related to channel width

    Args:
        channel_amount (-): Number of channels
        w_channel (m): Channel width
        w_channel_spacing (m): Distance between channels
        inlet_manifold_width_factor (-): Factor multiplied with channel width, that determines margin width

    Returns:
        inlet_manifold_width [m]: The total width of the inlet manifold
    """
    total_channel_width = channel_amount*w_channel + (channel_amount-1)*w_channel_spacing # [m]
    inlet_manifold_width_margin = inlet_manifold_width_factor * w_channel # [m]

    return total_channel_width + 2 *inlet_manifold_width_margin

def inlet_manifold_length(w_inlet_manifold, inlet_manifold_length_factor):
    """Determine inlet manifold length by multiplying it with a preset linear factor

    Args:
        w_inlet_manifold (m): Width of inlet manifold
        inlet_manifold_length_factor (-): Multiplication factor to determine inlet length (large enough to divide and settle fluid)

    Returns:
        inlet_manifold_length (m): The length of the inlet manifold
    """
    return w_inlet_manifold * inlet_manifold_length_factor # [m] Length of the inlet manifold

def total_chip_width(w_inlet_manifold, w_outer_margin, w_throat, AR_exit):
    """Calculate total chip width, based on either the total width of the inlet manifold or the nozzle exit width, depending on which is larger.
    A preset margin is added around the widest area (for a rectanguler chip)

    The outer chip width must check whether the nozzle exit or exit width is more important.

    Args:
        inlet_manifold_width
        w_outer_margin (m): Spacing around the outside of outer channels
        w_throat (m): Throat width
        AR_exit (-): Exit area ratio

    Returns:
        w_chip: (m) Total width of chip
    """


    w_nozzle_exit =w_throat*AR_exit # [m] Nozzle exit width

    return 2*w_outer_margin +  np.maximum( w_inlet_manifold , w_nozzle_exit )# [m] Total width of chip

def total_chip_length(l_inlet_manifold, l_channel, l_outlet):
    """ Return total channel length

    Args:
        l_inlet (m): Inlet length (manifold, etc.)
        l_channel (m): Heating channel length
        l_outlet (m): Exit manifold + nozzle length

    Returns:
        l_chip: (m) Total chip length
    """
    return l_inlet_manifold + l_channel + l_outlet

def outlet_length(w_channel, w_channel_spacing, channel_amount, convergent_half_angle, w_throat, divergent_half_angle, AR_exit, l_exit_manifold=0.0):
    """Calculate the outlet length (entire nozzle + exit manifold) of a channel for the purpose of determining the entire rectangular chip size

    Args:
        w_channel (m): Channel width
        w_channel_spacing (m): Spacing between channels
        channel_amount (-): Amount of parallel channels
        convergent_half_angle (rad): Half-angle of convergent nozzle part
        w_throat (m): width
        divergent_half_angle (rad): Half-angle of divergent nozzle part
        AR_exit (-): Exit area ratio
        l_exit_manifold (m, optional): Length of exit manifold, before convergent nozzle section begins. Defaults to 0.

    Returns:
        l_outlet (m): Total outlet length
    """
    # The exit manifold length determines the convergent part of the nozzle length, if a convergent angle is chosen to be fixed
    w_exit_manifold = channel_amount * w_channel + (channel_amount-1) * w_channel_spacing # [m] Total width of exit manifold (before feeding into nozzle)
    # L is the length of the virtual right-angled triangle that would appear if the nozzle was 0 micron wide.
    L_conv = 0.5*w_exit_manifold/math.tan(convergent_half_angle) # [m] Virtual nozzle length
    # Now the actual length
    l_convergent = L_conv * ( 1 - w_throat/w_exit_manifold ) # [m] Length of virtual convergent triangle up to throat 
    # The same process but for the exit
    L_div = 0.5*w_throat*AR_exit / math.tan(divergent_half_angle)
    l_divergent = L_div * ( 1 - 1/AR_exit) # [m] Length form virtual divergent triangle from the throat 

    return l_convergent + l_divergent + l_exit_manifold # [m] Total outlet length

def basic_substrate_heat_loss(T_top, kappa, emissivity, thickness, A_substrate):
    """Determine the heat loss through the substrate the chip is mounted on, by assuming it only conducts through it and radiates away
        Equations are solved by equation conduction and radiation

    Args:
        T_top (K): Temperature of Top of substrate
        kappa (W/(m*K)): Thermal conductivity of substrate
        emissivity (-): Emissivity of substrate
        thickness (m): Thickness of substrate
        A_substrate (m^2): Area of substrate (assumed to chip area connected to it, and area through which it radiates are the same)

    Returns:
        P_loss: (W) Heat loss through bottom of substrate
    """

    # T is the unknown temperature of the bottom side of the substrate. For the correct solution, the radiation and conduction are equal
    func_radiation = lambda T: emissivity*float(T)**4*stefan_boltzmann # [W/m^2] Radiative heat flux
    func_conduction = lambda T: kappa*(T_top - T) / thickness # [W/m^2] 
    func_equilibrium = lambda T: func_radiation(T) - func_conduction(T) # [W/m^2] Conservation of energy 
    
    sol = scipy.optimize.root_scalar(func_equilibrium,bracket=([240,T_top]))
    T_bottom = sol.root # [K] The root that puts heat fluxes in equilibrium
    P_loss = func_conduction(T_bottom) * A_substrate # [W] The actual power loss (could also be multiplied with radation heat flux)
    return P_loss # [W]

def area_ratio_contraction(w_inlet_manifold, w_channel, channel_amount):
    """Smaller area of channel divided by larger area of channel, which for the constant depth case is the ratio of widt

    Args:
        w_inlet_manifold (m): Width of inlet manifold (larger section)
        w_channel (m): Channel width in smaller section
        channel_amount (0): Amount of channels

    Returns:
        A_s/A_l : Area ratio of the sudden contraction
    """
    return w_channel*channel_amount/w_inlet_manifold # [-] Area ratio of smaller section A_s after contraction divided by larger section A_l before contraction