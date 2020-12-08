"""In this file the complete models from chamber inlet to nozzle outlet which require only 0-dimensional geometry inputs are considered
"""
import numpy as np
from scipy.optimize import optimize, root_scalar, minimize
from basic.IRT_corrections import hydraulic_diameter

from basic.chamber import T_wall_from_heat_flow, h_conv_from_Stanton, hydraulic_diameter_rectangular, radiation_loss, required_chip_area, required_heater_area, required_power, wetter_perimeter_rectangular, convective_heat_flow, mass_flow, velocity_from_mass_flow, ideal_enthalpy_change
from thermo.convection import Nu_DB, Nusselt_Dittus_Boelter, Stanton_from_Nusselt_and_velocity, Stanton_from_Nusselt_func_and_velocity
from thermo.prop import FluidProperties
import thermo.convection as conv
import basic.IRT_corrections as IRTc
import basic.IRT as IRT


def NIT_Rajeev(F_desired, p_inlet, h_throat, w_throat, throat_roc, AR_exit, p_back, divergence_half_angle, fp: FluidProperties, heating=True):
    pass

def engine_performance_from_F_and_T(F_desired, p_chamber, T_chamber, AR_exit, p_back, fp: FluidProperties):
    """ Returns Nozzle Inlet State based on IRT

    Args:
        F_desired (N): Desired thrust
        p_chamber (Pa): Chamber thrust/nozzle inlet pressure
        T_chamber (K): Chamber temperature
        AR_exit (-): A_exit / A_throat
        p_back (Pa): Back pressure
        fp (FluidProperties): Used to access fluid properties

    Raises:
        Exception: Raised when no solution is found for given input

    Returns:
        dict{A_throat [m^2], m_dot [kg/s]}: Throat area and mass flow
    """
    # Default range for temperature for root-finding algorithm
    A_low = 1e-20 # [K]
    A_high = 1e-6 # [K]

    R = fp.get_specific_gas_constant() # [J/(kg*K)]
    gamma = fp.get_specific_heat_ratio(T_chamber,p_chamber) # [-] Specific heat ratio
    # If x is zero, then the desired thrust is found. Gamma is changed depending on temperature as well
    x = lambda A_t: F_desired - IRT.get_engine_performance(p_chamber=p_chamber, T_chamber=T_chamber, A_throat=A_t, \
                        AR_exit=AR_exit,p_back=p_back, gamma=gamma,R=R)['thrust']

    root_result = root_scalar(x, bracket=[A_low,A_high], xtol=1e-12
    )

    if root_result.converged:
        A_throat = root_result.root # [m^2]
    else:
        raise Exception("No solution found for given temperature")

    # Now the chamber temperature is known, return the unknown parameters

    ep = IRT.get_engine_performance(p_chamber=p_chamber, T_chamber=T_chamber, A_throat=A_throat, \
                        AR_exit=AR_exit,p_back=p_back, gamma=fp.get_specific_heat_ratio(T_chamber,p_chamber),R=R)
    # Add throat area to dictionary
    ep['A_throat'] = A_throat # This is usually an input for get_engine_performance, but now added to make it one complete solution
    #ep['Isp'] = IRT.effective_ISP(thrust=F_desired,m_dot=ep['m_dot']) # [s] Effective specific impulse

    return ep

def DB_IRT(F_desired, p_inlet, T_inlet, T_chamber, L_channel, w_channel, h_channel, w_throat, h_throat, AR_exit, p_back, fp: FluidProperties, heating=True):


    inlet_state = engine_performance_from_F_and_T(F_desired=F_desired, p_chamber=p_inlet, T_chamber=T_chamber, AR_exit=AR_exit, p_back=p_back, fp=fp)
    A_throat = inlet_state['A_throat'] # [m^2] Throat area necessary for thrust F, given other conditions
    m_dot = inlet_state['m_dot'] # [kg/s] Mass flow determined necessary for thrust F, given other conditions

    # Report on flow conditions and achieved throat area
    print("Throat area: {:3.4f} um^2".format(A_throat*1e12))
    print("Mass flow: {:3.4f} mg/s".format(m_dot*1e6))

    # Calculate resultant throat width
    w_throat = A_throat / h_throat # [m] Throat width

    # Calculate what the heat flux is for the given conditions
    A_channel = w_channel * h_channel # [m^2] Channel cross-section
    rho_inlet = fp.get_density(T=T_inlet, p=p_inlet) # [kg/m^3]
    u_inlet = velocity_from_mass_flow(m_dot=m_dot, rho=rho_inlet, A=A_channel) # [m/s]
    print("u_inlet: {:3.9f} m/s".format(u_inlet))
    wetted_perimeter = wetter_perimeter_rectangular(w_channel = w_channel, h_channel=h_channel) # [m]
    D_hydraulic = hydraulic_diameter(A=A_channel,wetted_perimeter=wetted_perimeter) # [m]

    # Now the nusselt number from Dittus Boelter can be calcuated, with the outlet temperature naturally being the chamber temperature from IRT
    # No pressure drop is assumed as well
    Nu = Nusselt_Dittus_Boelter(T_inlet=T_inlet, T_outlet=T_chamber, p =p_inlet, D_hydraulic=D_hydraulic, L_channel=L_channel, u=u_inlet, fp=fp, heating=heating,supressExceptions=True) # [-] Nusselt number
    # Reference temperature for Nu must be the same as for St
    T_bulk = (T_inlet+T_chamber)/2 # [K] Bulk temperature is the reference temperature for DB
    St = Stanton_from_Nusselt_and_velocity(Nu=Nu, T_ref=T_bulk, p_ref=p_inlet, u=u_inlet, L_ref=D_hydraulic, fp=fp) # [-] Stanton number
    print("Stanton: {:3.9f} ".format(St))
    h_conv = h_conv_from_Stanton(Stanton=St, u=u_inlet, T_ref=T_bulk, p_ref=p_inlet, fp=fp) # [W/(m^2*K)] Heat transfer coefficient for conduction

    # Now everything is known, only one wall temperature will result in the desired heat flow
    # From the change in enthalpy and mass flow the required heat flow can be found
    delta_h = ideal_enthalpy_change(T_inlet=T_inlet, p_inlet=p_inlet, T_outlet=T_chamber, p_outlet=p_inlet,fp=fp) # [J/kg] Change in enthalpy across chamber
    Q_dot = -delta_h*m_dot # [W] Required heat flow from wall to fluid (is negative because it flows AWAY from the wall to the fluid)
    A_wall = L_channel * w_channel # [m^2] Channel wall through which the heat is conducted (one-sided heating is assumed)
    T_wall = T_wall_from_heat_flow(Q_dot=Q_dot, heat_transfer_coefficient=h_conv,T_ref=T_bulk, A_wall=A_wall)
    print("Required wall temperature: {:4.3f} K".format(T_wall))

    return {    'm_dot': m_dot, 
                'A_throat': A_throat,
                'u_inlet': u_inlet,
                'T_wall': T_wall}


# def DB_Rajeev(T_inlet, T_wall, p_inlet, u, L_channel, w_channel, h_channel, w_throat, h_throat, throat_roc, AR_exit, p_back, divergence_half_angle, fp: FluidProperties, heating=True):
    
#     # First, iterations must be performed to find outlet temperature that results from the provided wall temperature. 
#     # This is because the Dittus Boelter relation works with the bulk temperature, which depends on outlet temperature.

#     # Use scipy root_scalar to converge to a result where outlet temperature matches with the bulk temperature (which must match)
#     assert( T_inlet < T_wall) # To avoid surprises when incorrect inputs are given
#     x = lambda T_guess : T_guess - DB_T_chamber_outlet(T_inlet=T_inlet,T_outlet_guess=T_guess,T_wall=T_wall,p_inlet=p_inlet, L_channel=L_channel, u=u, w_channel=w_channel, h_channel=h_channel, fp=fp)
#     root_result = root_scalar(x, bracket=[T_inlet, T_wall], xtol=0.001) # [K] Converge to outlet temperature with 0.01K precision
#     if root_result.converged:
#         T_outlet= root_result.root # [K] Result if converged
#     else:
#         raise Exception("Root finding of outlet temperature did not converge")
#     # Capture engine performance from Rajeev
#     ep = IRTc.Rajeev_complete(p_chamber=p_inlet, T_chamber=T_outlet,w_throat=w_throat,h_throat=h_throat, throat_roc=throat_roc, AR_exit=AR_exit, p_back=p_back, divergence_half_angle=divergence_half_angle, fp=fp, is_cold_flow=False)

#     # Some results to verify if models are in agreement (probably not as this stage!)
#     A_channel = w_channel * h_channel # [m^2] Channel cross-sectional area
#     rho_inlet = fp.get_density(T=T_inlet,p=p_inlet)
#     m_dot_input = mass_flow(A=A_channel,u=u, rho=rho_inlet)
#     print("Mass flow input: {:3.5f} mg/s".format(m_dot_input*1e6))
#     print("Ideal mass flow: {:3.5f} mg/s".format(ep['m_dot_ideal']*1e6))
#     print("Real mass flow {:3.5f} mg/s".format(ep['m_dot_real']*1e6))

# def DB_T_chamber_outlet(T_inlet, T_outlet_guess, T_wall, p_inlet, L_channel,u, w_channel, h_channel, fp: FluidProperties):
#     """ WARNING: do not use function in isolation, but as part of iteration that converges T_outlet_guess towards actual T_outlet
#     Intermediate function in the model uses Dittus Boelter relation . It returns the chamber outlet temperature, based on a guess of this temperature,
#     which means that the result is only correct if T_outlet_guess and T_outlet match within an acceptable margin of error. This function helps to iterate and converge to the correct value

#     Args:
#         T_inlet (K): Inlet temperature
#         T_outlet_guess (K): Guess of the outlet temperature
#         T_wall (K): Wall temperature
#         p_inlet (Pa): Temperature at inlet (no pressure drop assumed over the rest of the chamber)
#         L_channel (m): Channel length
#         u (m/s): Flow velcoity
#         w_channel (m): Channel width
#         h_channel (m): Channel height
#         fp (FluidProperties): Object to access properties of fluid with

#     Returns:
#         T_outlet (K): Outlet temperature
#     """

#     A_channel = w_channel * h_channel # [m^2] Cross-sectional diameter of chamber
#     wetted_perimeter = wetter_perimeter_rectangular(w_channel=w_channel, h_channel=h_channel) # [m] Perimeter of channel
#     D_hydr_chamber = IRTc.hydraulic_diameter(A=A_channel, wetted_perimeter=wetted_perimeter) # [m] Hydraulic diameter of channel
#     Nu = conv.Nusselt_Dittus_Boelter(T_inlet, T_outlet_guess, p=p_inlet, D_hydraulic=D_hydr_chamber, L_channel=L_channel, u=u, fp=fp,supressExceptions=True) # [-] Nusselt number of of the flow
#     # It is important to use the exact same reference temperature and pressurefor Stanton, as for the case in which Nusselt was determined.
#     # Pressure is in all cases assumed to be pressre at chamber inlet calculations
#     # Similarly, the reference length must be the same so Nu and Re are in agreement and lengths cancel out. Which is the hydraulic diameter in this case
#     T_bulk = (T_inlet+T_outlet_guess)/2 # [K] For Dittus Boelter this is the bulk temperature
#     St = conv.Stanton_from_Nusselt_and_velocity(Nu=Nu,T_ref=T_bulk, p_ref=p_inlet, u=u, L_ref=D_hydr_chamber, fp=fp) # [-]
#     h_conv = h_conv_from_Stanton(Stanton=St, u=u, T_ref=T_bulk, p_ref=p_inlet, fp=fp) # [W/(kg*m^2)]

#     # Now the enthalpy increase of the fluid due to the heat flow into the fluid must be determined
#     A_wall = L_channel*w_channel # [m^2] Chamber wall, through which the fluid is heated
#     Q_wall = convective_heat_flow(heat_transfer_coefficient=h_conv, T_wall=T_wall, T_ref=T_bulk, A_wall=A_wall) # [W] Heat flow towards wall (is negative when heat flows away from wall to the flow)
    
#     rho_inlet = fp.get_density(T=T_inlet,p=p_inlet) # [kg/m^3]
#     m_dot = mass_flow(A=A_channel,u=u,rho=rho_inlet) # [kg/s]
#     # Q_wall is negative as heat flows away from T_wall to T_bulk, so enthalpy increase is -Q_wall
#     delta_h = -Q_wall/m_dot # [J/kg] Increase in specific enthalpy of the flow at the channel outlet. 

#     print("Delta h: {:3.5f}".format(delta_h))

#     # The initial enthalpy must be known to determine the thermodynamic state at the end
#     h_inlet = fp.get_enthalpy(T_inlet,p_inlet) # [J/kg]
#     print("h_inlet {:3.5f}".format(h_inlet))
#     h_outlet = h_inlet + delta_h # [J/kg]
#     T_outlet = fp.get_temperature_from_enthalpy_and_pressure(h=h_outlet,p=p_inlet) # [K] The temperature at the outlet, based on the guess.
#     # Note, this is NOT the correct temperature unless T_outlet_guess and T_outlet match within an acceptable margin of error 
#     # Because the reference temperature T_bulk depends on T_outlet_guess
#     print (T_outlet)

#     return T_outlet # [K] Outlet temperature (based on guess of outlet temperature)

def chamber_performance_from_Nu(Nu_func, T_inlet, T_chamber, T_ref, T_wall, p_ref, m_dot, A_channel, L_ref, fp: FluidProperties):
    """ Function that calculates the power consumption and heating area for a specific chamber

    Args:
        Nu_func (-): Nusselt function, that implement the emperical relation of choice
        T_inlet (K): Chamber inlet temperature
        T_chamber (K): Chamber outlet temperature (same as T_chamber for IRT)
        T_ref (K): Reference temperature for the Nusselt relation and flow similary parameters
        T_wall (K): Wall temperature
        p_ref (Pa): Chamber pressure (no pressure drop assumed)
        m_dot (kg/s): Mass flow
        A_channel (m^2): Cross-sectional area of the chamber, through which the fluid flows
        L_ref (m): Reference length for Nusselt relation and flow similarty parameters
        fp (FluidProperties): Object from which Fluid Properties are determined

    Returns:
        dictionary with heater area, power required to heat up flow and Nusselt number
    """
    ## Make sure all parameters are calculated at the same reference state (including velocity!)
    # Pr and Re parameters needed for most Nusselt relation
    rho_ref = fp.get_density(T=T_ref, p=p_ref) # [kg/m^3] Reference density
    u_ref = velocity_from_mass_flow(A=A_channel, m_dot=m_dot, rho=rho_ref) # [m/s] Speed at reference state
    print("u_ref {} m/s".format(u_ref))
    Re_ref = fp.get_Reynolds_from_mass_flow(T=T_ref, p=p_ref, L_ref=L_ref, m_dot=m_dot, A=A_channel) # [-] Reynolds number at reference state
    Pr_ref = fp.get_Prandtl(T=T_ref, p=p_ref) # [-] Prandtl number at reference state
    # Now the Nusselt can be determined
    Nusselt = Nu_func(args={'Re':Re_ref, 'Pr': Pr_ref}) # [-] Nusselt number at given state (used for plotting purposes)
    Stanton = Stanton_from_Nusselt_func_and_velocity(Nu_func=Nu_func, m_dot=m_dot, A=A_channel, T_ref=T_ref, p_ref=p_ref, L_ref=L_ref, fp=fp) # [-] Stanton number at reference state
    h_conv = h_conv_from_Stanton(Stanton=Stanton, u=u_ref, T_ref=T_ref, p_ref=p_ref, fp=fp)
    # Now determine how much energy must be convected
    delta_h = ideal_enthalpy_change(T_inlet=T_inlet, p_inlet=p_ref, T_outlet=T_chamber, p_outlet=p_ref, fp=fp) # [J/kg]
    Q_dot = required_power(m_dot=m_dot, delta_h=delta_h) # [W] Required power to achieve delta_h
    A_heater = required_heater_area(Q_dot=Q_dot, h_conv=h_conv, T_wall=T_wall, T_ref=T_ref) # [m^2]
    assert(A_heater>0)
    # Return a dictionary with interesting values
    return {'A_heater': A_heater,
            'Q_dot': Q_dot,
            'Nusselt': Nusselt,
            'Re_ref': Re_ref,
            'Pr_ref': Pr_ref}

def total_power_single_channel(T_wall, w_channel, Nu_func, T_inlet, T_chamber, T_ref, p_ref, m_dot, h_channel, w_channel_margin, emmisivity, fp):
    """ Function that calculates the total power consumption of a specific chamber, in order to optimize the chamber

    Args:
        T_wall (K): Wall temperature
        w_channel (m): Channel width
        Nu_func (-): Nusselt function
        T_inlet (K): Chamber inlet temperature
        T_chamber (K): Chamber outlet temperature (same as T_c in IRT)
        T_ref (K): Reference temperature for the Nusselt relation and flow similary parameters
        p_ref (Pa): Reference pressure for the Nusselt relation and flow similary parameters (same as inlet pressure as no pressure drop is assumed)
        m_dot (kg/s): Mass flow
        h_channel (m): Channel height
        w_channel_margin (m): The amount of margin around the chamber for structural reasons. Important because it also radiates heat
        emmisivity (-) Emmisivity of the chamber material (for calculation of radation losses)
        fp (-  ): [description]
    """

    ## NOTE: just a test, assumptions are a bit sloppy, and this is just to test optimization approach at this stage

    A_channel = w_channel * h_channel # [m^2] Channel area
    print("A_channel: {}".format(A_channel))
    wetted_perimeter = wetter_perimeter_rectangular(w_channel=w_channel, h_channel=h_channel) # [m] Distance of channel cross-section in contact with fluid
    print("wetted_perimeter: {}".format(wetted_perimeter))
    # Reference length is hydraulic diameter
    D_hydraulic = hydraulic_diameter_rectangular(w_channel=w_channel, h_channel=h_channel) # [m] Hydraulic diameter

    cp = chamber_performance_from_Nu(Nu_func=Nu_func, T_inlet=T_inlet, T_chamber=T_chamber, T_ref=T_ref, T_wall=T_wall, p_ref=p_ref, m_dot=m_dot, A_channel=A_channel, L_ref=D_hydraulic, fp=fp)
    ## Assume single-sided heater heat chamber material through-out at T-wall, and all along the wetted perimeter heat flows equally to teh propellant
    A_heater = cp['A_heater'] # [m^2]
    assert(A_heater>0)
    L_channel = A_heater / wetted_perimeter # [m] Required length of channel to achieve required heat flow
    print("L_channel: {}".format(L_channel))
    A_chip = required_chip_area(L_channel=L_channel, w_channel=w_channel, w_channel_margin=w_channel_margin) # [m^2] Assume radiation occurs from silicon on just a single side of the heater
    print("A_microheater: {}".format(A_chip))
    P_radiation = radiation_loss(T=T_wall, A=A_chip, emmisivity=emmisivity) # [W] Radiation loss
    P_total = P_radiation + cp['Q_dot'] # [W] Total power consumption, the variable to minimize
    print("Heat loss {}".format(P_radiation/P_total))

    return P_total

def minimize_total_power_single_channel(T_wall_superheat_min, T_wall_superheat_max, w_channel_min, w_channel_max, Nu_func, T_inlet, T_chamber, T_ref, p_ref, m_dot, h_channel, w_channel_margin, emmisivity, fp):
    """Attempts to find the lowest total power consumption by varying wall superheat and channel width

    Args:
        T_wall_superheat (K): Superheat of wall in relation to T_chamber (min and max for optim. algorithm)
        w_channel (m): Channel width (for optimization algorithm)
        Nu_func (-): Nusselt function
        T_inlet (K): Chamber inlet temperature
        T_chamber (K): Chamber outlet temperature (same as T_c in IRT)
        T_ref (K): Reference temperature for the Nusselt relation and flow similary parameters
        p_ref (Pa): Reference pressure for the Nusselt relation and flow similary parameters (same as inlet pressure as no pressure drop is assumed)
        m_dot (kg/s): Mass flow
        h_channel (m): Channel height
        w_channel_margin (m): The amount of margin around the chamber for structural reasons. Important because it also radiates heat
        emmisivity (-) Emmisivity of the chamber material (for calculation of radation losses)
        fp (-  ): [description]
    """
    ## Function to minimize, obscuring the constant values, and exposing only the the variable to optimize
    def min_fun_P_total(x):
        # x[0] is T_wall_superheat [K]
        # x[1] is w_channel # [m]
        T_wall = float(x[0]) + T_chamber # [K]
        print("T_wall: {}".format(T_wall))
        w_channel = float(x[1]) # [m]
        print("w_channel: {}".format(w_channel))
        return total_power_single_channel(T_wall=T_wall, w_channel=w_channel, Nu_func=Nu_func, T_inlet=T_inlet, T_chamber=T_chamber, \
            T_ref=T_ref, p_ref=p_ref, m_dot=m_dot, h_channel=h_channel, w_channel_margin=w_channel_margin, emmisivity=emmisivity, fp=fp)

    ## Simply guess that the optimum is in the middle
    T_wall_superheat_guess = (T_wall_superheat_max + T_wall_superheat_min) / 2
    w_channel_guess = (w_channel_max + w_channel_min) / 2
    x_guess = np.array([T_wall_superheat_guess, w_channel_guess])
    x_bounds = np.array([[T_wall_superheat_min, T_wall_superheat_max],[w_channel_min,w_channel_max]])

    x_result = minimize(fun=min_fun_P_total, x0=x_guess, bounds=x_bounds)
    return x_result


def minimize_two_phase_Nusselt(T_wall_superheat_min, T_wall_superheat_max, w_channel_min, w_channel_max, Nu_func, T_inlet, T_chamber, T_ref, p_ref, m_dot, h_channel, w_channel_margin, fp):
    """Attempts to find the lowest total heater area by varying wall superheat and channel width

    Args:
        T_wall_superheat (K): Superheat of wall in relation to T_chamber (min and max for optim. algorithm)
        w_channel (m): Channel width (for optimization algorithm)
        Nu_func (-): Nusselt function
        T_inlet (K): Chamber inlet temperature
        T_chamber (K): Chamber outlet temperature (same as T_c in IRT)
        T_ref (K): Reference temperature for the Nusselt relation and flow similary parameters
        p_ref (Pa): Reference pressure for the Nusselt relation and flow similary parameters (same as inlet pressure as no pressure drop is assumed)
        m_dot (kg/s): Mass flow
        h_channel (m): Channel height
        w_channel_margin (m): The amount of margin around the chamber for structural reasons. Important because it also radiates heat
        emmisivity (-) Emmisivity of the chamber material (for calculation of radation losses)
        fp (-  ): [description]
    """

    def min_func_area(x):
        # x[0] is T_wall_superheat [K]
        # x[1] is w_channel # [m]
        T_wall = float(x[0]) + T_chamber # [K]
        print("T_wall: {}".format(T_wall))
        w_channel = float(x[1]) # [m]
        A_chip = 1 # [m^2]
        return A_chip

    ## Simply guess that the optimum is in the middle
    T_wall_superheat_guess = (T_wall_superheat_max + T_wall_superheat_min) / 2
    w_channel_guess = (w_channel_max + w_channel_min) / 2
    x_guess = np.array([T_wall_superheat_guess, w_channel_guess])
    x_bounds = np.array([[T_wall_superheat_min, T_wall_superheat_max],[w_channel_min,w_channel_max]])

    x_result = minimize(fun=min_func_area, x0=x_guess, bounds=x_bounds)
    return x_result

def two_phase_single_channel(T_wall, w_channel, Nu_func_gas, Nu_func_liquid, T_inlet, T_chamber, p_ref, m_dot, h_channel, fp: FluidProperties, print_info=True):
    """ Function that calculates the total power consumption of a specific chamber, in order to optimize the chamber

    Args:
        T_wall (K): Wall temperature
        w_channel (m): Channel width
        Nu_func_gas (-): Nusselt function for gas phase
        Nu_func_liquid (-) Nusselt function for liquid phase
        T_inlet (K): Chamber inlet temperature
        T_chamber (K): Chamber outlet temperature (same as T_c in IRT)
        p_ref (Pa): Reference pressure for the Nusselt relation and flow similary parameters (same as inlet pressure as no pressure drop is assumed)
        m_dot (kg/s): Mass flow
        h_channel (m): Channel height
        w_channel_margin (m): The amount of margin around the chamber for structural reasons. Important because it also radiates heat
        fp (-  ): [description]
        print_info(Bool): for debugging purposes
    """

    # Calculate saturation temperature, to determine where transition from gas to liquid occurs
    T_sat = fp.get_saturation_temperature(p=p_ref) # [K]
    # Sanity check on input
    assert(T_chamber > T_sat)
    assert(T_wall > T_chamber)
    

    # Calculate the two reference temperatures for the separated phases
    T_bulk_gas = (T_sat + T_chamber) / 2 # [K] Bulk temperature gas phase
    T_bulk_liquid_multi = (T_inlet + T_sat) / 2# [K] Bulk temperature of liquid and multi-phase flow
    # Calculate the density at these reference points
    rho_bulk_gas = fp.get_density(T=T_bulk_gas, p=p_ref) # [kg/m^3]
    rho_bulk_liquid_multi = fp.get_density(T=T_bulk_liquid_multi, p=p_ref) # [kg/m^3]

    # Channel geometry
    A_channel = w_channel * h_channel # [m^2] Area through which the fluid flows
    wetted_perimeter = wetter_perimeter_rectangular(w_channel=w_channel, h_channel=h_channel) # [m] Distance of channel cross-section in contact with fluid
    D_hydraulic = hydraulic_diameter_rectangular(w_channel=w_channel, h_channel=h_channel) # [m] Hydraulic diameter

    # Flow similarity parameters (for debugging and Nu calculatoin purposes)
    Re_bulk_gas = fp.get_Reynolds_from_mass_flow(m_dot=m_dot, p =p_ref, T=T_bulk_gas, L_ref=D_hydraulic, A=A_channel) # [-] Bulk Reynolds number in the gas phase
    Re_bulk_liquid_multi = fp.get_Reynolds_from_mass_flow(m_dot=m_dot, p=p_ref, T=T_bulk_liquid_multi, L_ref=D_hydraulic, A=A_channel) # [-] Bulk Reynolds number in the liquid/multi-phase
    Pr_bulk_gas = fp.get_Prandtl(T=T_bulk_gas, p=p_ref) # [-] Prandtl number in the gas phase
    Pr_bulk_liquid_multi = fp.get_Prandtl(T=T_bulk_liquid_multi, p=p_ref) # [-] Prandtl number in liquid/multi-phase
    Bo_sat= fp.get_Bond_number(p_sat=p_ref, L_ref=D_hydraulic) # [-] Bond number at saturation pressure (assumed to be p_ref)

    # Calculate Nusselt number in both sections
    args_gas = {'Re': Re_bulk_gas, 'Pr': Pr_bulk_gas, 'Bo': Bo_sat} # Arguments for Nusselt function (gas phase) 
    args_liquid_multi = {'Re': Re_bulk_liquid_multi, 'Pr': Pr_bulk_liquid_multi} # Arguments for Nusselt function (liquid phase)
    Nu_gas = Nu_func_gas(args=args_gas)
    Nu_liquid_multi = Nu_func_liquid(args=args_liquid_multi)
    # Calculate Stanton number in both sections
    St_gas = Stanton_from_Nusselt_and_velocity(Nu=Nu_gas, T_ref=T_bulk_gas, p_ref=p_ref, L_ref=D_hydraulic, m_dot=m_dot, A=A_channel, fp=fp) # [-] Stanton number in gas phase
    St_liquid_multi = Stanton_from_Nusselt_and_velocity(Nu_liquid_multi, T_ref=T_bulk_liquid_multi, p_ref=p_ref, L_ref=D_hydraulic, m_dot=m_dot, A=A_channel, fp=fp) # [-] Stanton number in liquid phase
    # Calculate velocity for convection parameter (bulk temp used as reference for phase)
    u_bulk_gas = velocity_from_mass_flow(A=A_channel, m_dot=m_dot, rho=rho_bulk_gas) # [m/s] Velocity at the gas bulk reference state
    u_bulk_liquid_multi = velocity_from_mass_flow(A=A_channel, m_dot=m_dot, rho=rho_bulk_liquid_multi) # [m/s] Velocity at the liquid/multi-phase bulk reference state
    # Convective parameter
    h_conv_gas = h_conv_from_Stanton(Stanton=St_gas, u=u_bulk_gas, T_ref=T_bulk_gas, p_ref=p_ref, fp=fp) # [W/(m^2*K)] Convective heat transfer coefficient at bulk gas state
    h_conv_liquid_multi = h_conv_from_Stanton(Stanton=St_liquid_multi, u=u_bulk_liquid_multi, T_ref=T_bulk_liquid_multi, p_ref=p_ref, fp=fp) # [W/(m^2*K)] Convective heat transfer coefficient at bulk liquid/multi-phase state
    # Required specific enthalpy change for heating the separate sections
    h_outlet = fp.get_enthalpy(T=T_chamber, p=p_ref) # [J/kg] Specific enthalpy at the outlet
    h_sat_gas = fp.get_saturation_enthalpy_gas(p=p_ref) # [J/kg] Specific enthalpy of saturated gas
    h_inlet = fp.get_enthalpy(T=T_inlet, p=p_ref) # [J/kg]
    # Required specific enthalpy increases
    delta_h_gas = h_outlet - h_sat_gas # [J/kg] Enthalpy increase needed to go from saturated gas to outlet enthalpy
    delta_h_liquid_multi = h_sat_gas - h_inlet # [J/k] Enthalpy increase needed to go from inlet enthalpy to saturated gas
    # Required power for those enthalpy changes at the given mass flow
    Q_dot_gas = required_power(m_dot=m_dot,delta_h=delta_h_gas) # [W]
    Q_dot_liquid_multi = required_power(m_dot=m_dot, delta_h=delta_h_liquid_multi) # [W]
    # Required heater area to achieve the required power
    A_heater_gas = required_heater_area(Q_dot=Q_dot_gas, h_conv=h_conv_gas, T_wall=T_wall, T_ref=T_bulk_gas) # [m^2]
    A_heater_liquid_multi = required_heater_area(Q_dot=Q_dot_liquid_multi, h_conv=h_conv_liquid_multi, T_wall=T_wall, T_ref=T_bulk_liquid_multi) # [m^2]
    # Required length to achieve desired area
    L_channel_gas = A_heater_gas/wetted_perimeter # [m] Length of channel after gas is saturated
    L_channel_liquid_multi = A_heater_liquid_multi/wetted_perimeter # [m] Length of channel after heater
    L_channel = L_channel_gas + L_channel_liquid_multi # [m]
    L_hydrodynamic_entrance = D_hydraulic * Re_bulk_liquid_multi * 0.04 # [m] Hydrodynamic entrance to estimate if the flow is fully developed


    assert(h_outlet > h_sat_gas)
    assert(h_sat_gas > h_inlet)

    if(print_info):
        print("\n--- SPECIFIC ENTHALPY AT DIFFERENT STATIONS ---")
        print("h_outlet: {:4.3f} J/kg".format(h_outlet))
        print("h_sat_gas: {:4.3f} J/kg".format(h_sat_gas))
        print("h_inlet: {:4.3f} J/kg".format(h_inlet))

        print("\n --- REQUIRED POWER ---")
        print("Q_dot_gas: {:2.5f} W".format(Q_dot_gas))
        print("Q_dot_liquid_multi: {:2.5f} W".format(Q_dot_liquid_multi))

        print("\n --- BULK GAS PHASE PARAMETERS --- ")
        print("u: {:3.2f} m/s".format(u_bulk_gas))
        print("Nu: {}".format(Nu_gas))
        print("Re: {}".format(Re_bulk_gas))
        print("Pr: {}".format(Pr_bulk_gas))
        print("St: {}".format(St_gas))
        print("Bo_sat: {}".format(Bo_sat))
        
        print("\n --- BULK LIQUID/MULTI-PHASE PARAMETERS ---")
        print("u: {:3.4f} m/s".format(u_bulk_liquid_multi))
        print("Nu: {}".format(Nu_liquid_multi))
        print("Re: {}".format(Re_bulk_liquid_multi))
        print("Pr: {}".format(Pr_bulk_liquid_multi))
        print("St: {}".format(St_liquid_multi))

        print("\n --- CHARACTERISTIC GEOMETRIC VALUES --- ")
        print("Hydrodynamic entance length: {:3.3f} micron".format(L_hydrodynamic_entrance*1e6))
        print("Hydraulic diameter: {:3.3f} micron".format(D_hydraulic*1e6))
        print("L/D: {:4.2f} ".format(L_channel/D_hydraulic))
        print("L/X_T {:4.2f}".format(L_channel/L_hydrodynamic_entrance))

        print("\n --- RESULTING GEOMETRY ---")
        print("Total length: {:3.3f} mm".format(L_channel*1e3))
        print("Length (liquid/multi): {:3.3f} mm".format(L_channel_liquid_multi*1e3))
        print("Length (gas): {:3.4f} mm".format(L_channel_gas*1e3))
        print("Relative length (gas) {:3.3f} \%".format(L_channel_gas/L_channel*100))

        ## Return a dictionary with results and interesting intermediate values
    res = { "L_channel": L_channel, # [m] Total length of channel
            "D_hydraulic": D_hydraulic, # [m] Hydraulic diameter of channel
            "Nu_liquid_multi": Nu_liquid_multi, # [-] Nusselt number of liquid/multi-phase flow
            "Pr_bulk_liquid_multi": Pr_bulk_liquid_multi, # [-] Prandlt number of liquid/multi-phase flow
            "Re_bulk_liquid_multi": Re_bulk_liquid_multi, # [-] Reynolds number of liquid/multi-phase flow
            "St_liquid_multi": St_liquid_multi, # [-] Stanton number of liquid/multi-phase flow
            "h_conv_liquid_multi": h_conv_liquid_multi, # [W/(m^2*K)] Heat transfer coefficient
            "A_heater_liquid_multi": A_heater_liquid_multi, # [m^2] Required heater area for liquid/multi-phase flow
            "L_channel_liquid_multi": L_channel_liquid_multi, # [m] Length of channel to get required heater area
            "u_bulk_liquid_multi": u_bulk_liquid_multi, # [m/s] Bulk flow velocity of liquid/multi-phase flow
            "rho_bulk_liquid_multi": rho_bulk_liquid_multi, # [kg/m^3] Bulk density of liquid/multi-phase flow
            "T_bulk_liquid_multi": T_bulk_liquid_multi, # [K] Bulk temperature of liquid/multi-phase flow
            "delta_h_liquid_multi": delta_h_liquid_multi, # [J/kg] Enthalpy change from inlet to saturated gas
            "Q_dot_liquid_multi": Q_dot_liquid_multi, # [W] Power required for enthalpy change
            ## Same thing but for gas values
            "Nu_gas": Nu_gas, # [-]
            "Pr_bulk_gas": Pr_bulk_gas, # [-]
            "Re_bulk_gas": Re_bulk_gas, # [-]
            "St_gas": St_gas, # [-]
            "h_conv_gas": h_conv_gas, # [W/(m^2*K)]
            "A_heater_gas": A_heater_gas, # [m^2]
            "L_channel_gas": L_channel_gas, # [m]
            "u_bulk_gas": u_bulk_gas, # [m/s]
            "rho_bulk_gas": rho_bulk_gas, # [kg/m^3]
            "T_bulk_gas": T_bulk_gas, # [K]
            "delta_h_gas": delta_h_gas, # [J/kg]
            "Q_dot_gas": Q_dot_gas, # [W]
            }
    return res

    
    