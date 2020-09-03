"""In this file the complete models from chamber inlet to nozzle outlet which require only 0-dimensional geometry inputs are considered
"""
from scipy.optimize import root_scalar
from basic.IRT_corrections import hydraulic_diameter

from basic.chamber import T_wall_from_heat_flow, h_conv_from_Stanton, wetter_perimeter_rectangular, convective_heat_flow, mass_flow, velocity_from_mass_flow, ideal_enthalpy_change
from thermo.convection import Nusselt_Dittus_Boelter, Stanton_from_Nusselt_and_velocity
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


def DB_Rajeev(T_inlet, T_wall, p_inlet, u, L_channel, w_channel, h_channel, w_throat, h_throat, throat_roc, AR_exit, p_back, divergence_half_angle, fp: FluidProperties, heating=True):
    
    # First, iterations must be performed to find outlet temperature that results from the provided wall temperature. 
    # This is because the Dittus Boelter relation works with the bulk temperature, which depends on outlet temperature.

    # Use scipy root_scalar to converge to a result where outlet temperature matches with the bulk temperature (which must match)
    assert( T_inlet < T_wall) # To avoid surprises when incorrect inputs are given
    x = lambda T_guess : T_guess - DB_T_chamber_outlet(T_inlet=T_inlet,T_outlet_guess=T_guess,T_wall=T_wall,p_inlet=p_inlet, L_channel=L_channel, u=u, w_channel=w_channel, h_channel=h_channel, fp=fp)
    root_result = root_scalar(x, bracket=[T_inlet, T_wall], xtol=0.001) # [K] Converge to outlet temperature with 0.01K precision
    if root_result.converged:
        T_outlet= root_result.root # [K] Result if converged
    else:
        raise Exception("Root finding of outlet temperature did not converge")
    # Capture engine performance from Rajeev
    ep = IRTc.Rajeev_complete(p_chamber=p_inlet, T_chamber=T_outlet,w_throat=w_throat,h_throat=h_throat, throat_roc=throat_roc, AR_exit=AR_exit, p_back=p_back, divergence_half_angle=divergence_half_angle, fp=fp, is_cold_flow=False)

    # Some results to verify if models are in agreement (probably not as this stage!)
    A_channel = w_channel * h_channel # [m^2] Channel cross-sectional area
    rho_inlet = fp.get_density(T=T_inlet,p=p_inlet)
    m_dot_input = mass_flow(A=A_channel,u=u, rho=rho_inlet)
    print("Mass flow input: {:3.5f} mg/s".format(m_dot_input*1e6))
    print("Ideal mass flow: {:3.5f} mg/s".format(ep['m_dot_ideal']*1e6))
    print("Real mass flow {:3.5f} mg/s".format(ep['m_dot_real']*1e6))

def DB_T_chamber_outlet(T_inlet, T_outlet_guess, T_wall, p_inlet, L_channel,u, w_channel, h_channel, fp: FluidProperties):
    """ WARNING: do not use function in isolation, but as part of iteration that converges T_outlet_guess towards actual T_outlet
    Intermediate function in the model uses Dittus Boelter relation . It returns the chamber outlet temperature, based on a guess of this temperature,
    which means that the result is only correct if T_outlet_guess and T_outlet match within an acceptable margin of error. This function helps to iterate and converge to the correct value

    Args:
        T_inlet (K): Inlet temperature
        T_outlet_guess (K): Guess of the outlet temperature
        T_wall (K): Wall temperature
        p_inlet (Pa): Temperature at inlet (no pressure drop assumed over the rest of the chamber)
        L_channel (m): Channel length
        u (m/s): Flow velcoity
        w_channel (m): Channel width
        h_channel (m): Channel height
        fp (FluidProperties): Object to access properties of fluid with

    Returns:
        T_outlet (K): Outlet temperature
    """

    A_channel = w_channel * h_channel # [m^2] Cross-sectional diameter of chamber
    wetted_perimeter = wetter_perimeter_rectangular(w_channel=w_channel, h_channel=h_channel) # [m] Perimeter of channel
    D_hydr_chamber = IRTc.hydraulic_diameter(A=A_channel, wetted_perimeter=wetted_perimeter) # [m] Hydraulic diameter of channel
    Nu = conv.Nusselt_Dittus_Boelter(T_inlet, T_outlet_guess, p=p_inlet, D_hydraulic=D_hydr_chamber, L_channel=L_channel, u=u, fp=fp,supressExceptions=True) # [-] Nusselt number of of the flow
    # It is important to use the exact same reference temperature and pressurefor Stanton, as for the case in which Nusselt was determined.
    # Pressure is in all cases assumed to be pressre at chamber inlet calculations
    # Similarly, the reference length must be the same so Nu and Re are in agreement and lengths cancel out. Which is the hydraulic diameter in this case
    T_bulk = (T_inlet+T_outlet_guess)/2 # [K] For Dittus Boelter this is the bulk temperature
    St = conv.Stanton_from_Nusselt_and_velocity(Nu=Nu,T_ref=T_bulk, p_ref=p_inlet, u=u, L_ref=D_hydr_chamber, fp=fp) # [-]
    h_conv = h_conv_from_Stanton(Stanton=St, u=u, T_ref=T_bulk, p_ref=p_inlet, fp=fp) # [W/(kg*m^2)]

    # Now the enthalpy increase of the fluid due to the heat flow into the fluid must be determined
    A_wall = L_channel*w_channel # [m^2] Chamber wall, through which the fluid is heated
    Q_wall = convective_heat_flow(heat_transfer_coefficient=h_conv, T_wall=T_wall, T_ref=T_bulk, A_wall=A_wall) # [W] Heat flow towards wall (is negative when heat flows away from wall to the flow)
    
    rho_inlet = fp.get_density(T=T_inlet,p=p_inlet) # [kg/m^3]
    m_dot = mass_flow(A=A_channel,u=u,rho=rho_inlet) # [kg/s]
    # Q_wall is negative as heat flows away from T_wall to T_bulk, so enthalpy increase is -Q_wall
    delta_h = -Q_wall/m_dot # [J/kg] Increase in specific enthalpy of the flow at the channel outlet. 

    print("Delta h: {:3.5f}".format(delta_h))

    # The initial enthalpy must be known to determine the thermodynamic state at the end
    h_inlet = fp.get_enthalpy(T_inlet,p_inlet) # [J/kg]
    print("h_inlet {:3.5f}".format(h_inlet))
    h_outlet = h_inlet + delta_h # [J/kg]
    T_outlet = fp.get_temperature_from_enthalpy_and_pressure(h=h_outlet,p=p_inlet) # [K] The temperature at the outlet, based on the guess.
    # Note, this is NOT the correct temperature unless T_outlet_guess and T_outlet match within an acceptable margin of error 
    # Because the reference temperature T_bulk depends on T_outlet_guess
    print (T_outlet)

    return T_outlet # [K] Outlet temperature (based on guess of outlet temperature)




