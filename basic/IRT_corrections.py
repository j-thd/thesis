"""File for applying corrections factors to basic IRT calculations. 

Corrections may stem from different sources, and the paper/book will be listed for each equation/calculation
"""

import math
import basic.IRT as IRT
from thermo.prop import FluidProperties

def divergence_loss_conical_2D(alpha):
    """Divergence loss for 2d conical nozzles according eq 2.17 from Makhan 2018

    Args:
        alpha (rad): nozzle divergence half angle
    """

    return math.sin(alpha)/alpha

## Commented out as it may not be necessary, and is not verified.
# def divergence_loss_conical_3D(alpha):
#     """Divergence loss for 3D conical nozzles according eq.28 from Makhan 2018

#     Args:
#         alpha (alpha): nozzle divergence half angle
#     """

#     return (1 + math.cos(alpha))/2

def viscous_loss(area_ratio, reynolds_throat_wall):
    """Viscous loss acoording to eq. 2.19 from Makhan

    Args:
        area_ratio (-): Exit area/throat area
        reynolds_throat_wall (-): Reynolds number at the throat wall
    """

    return 17.6*math.exp(0.0032*area_ratio)/math.sqrt(reynolds_throat_wall)

def Reynolds_throat_wall_cold(reynolds_throat):
    """Estimation of Reynolds at throat wall for cold flow from Spisz1965 eq. 18 and 19

    Args:
        reynolds_throat (-): Reynolds number at the throat
    """

    return reynolds_throat * 0.857**(5/3)

def Reynolds_throat_wall_hot(reynolds_throat):
    """Estimation of Reynolds at throat wall for cold flow from Spisz1965 eq. 18 and 19

    Args:
        reynolds_throat (-): Reynolds number at the throat
    """

    return reynolds_throat * 1.388**(5/3)


def throat_boundary_loss(gamma, reynolds_throat, throat_radius, throat_roc):
    """Throat boundary loss according to eq. 2.20, 2.21 from Makhan2018
        Made for axi-symmetric nozzles with no sharp throat

    Args:
        gamma (-): Specific heat ratio
        reynolds_throat (-): Reynolds number at the throat
        throat_radius (m): Throat radius (or half of hydraulic diameter)
        throat_roc (m): Radius of curvature at the thraot
    """

    Re_mod = reynolds_throat * math.sqrt(throat_radius/throat_roc) # [-] Modified reynolds number at the throat (2.21 - Makhan2018)

    #Intermediate variables for readability
    a = ((gamma+1)/2)**(3/4)
    b = (72-32*math.sqrt(6)) / (3 * (gamma+1))
    c = 4*math.sqrt(6)/3
    d = b + c
    e = 2*math.sqrt(2)*(gamma-1)*(gamma+2)
    f = e / ( 3 * math.sqrt(gamma+1) )

    return 1 - a*d/math.sqrt(Re_mod) + f / Re_mod

def reynolds(m_dot, A, D_hydraulic, viscosity):
    """Returns the reynolds number based on the local conditions

    Args:
        m_dot (kg/s): mass flow
        A (m^2): cross-sectional area
        D_hydraulic (m): Hydraulic diameter
        viscosity (Pa*s): Dynamic viscosity

    Returns:
        [-]: Reynolds numbers (dimensionless flow characteristic)
    """
    return m_dot*D_hydraulic/(A*viscosity)

def hydraulic_diameter(A, wetted_perimeter):
    """ Returns the hydraulic diameter, often used to replace the regular diameter when a channel is circular and a characteristic length is needed

    Args:
        A (m^2): cross-sectional area of channel
        wetted_perimeter (m): length of boundary of the rea in contact with fluid (at this time all parts of the code assume full contact)

    Returns:
        [m]: Hydraulic diameter
    """

    return 4 * A / wetted_perimeter

def Rajeev_complete(p_chamber, T_chamber, w_throat, h_throat, throat_roc, AR_exit, p_back, divergence_half_angle, fp: FluidProperties, is_cold_flow):
    """ Function that implements all corrections proposed by Makhan2018

    Args:
        p_chamber (Pa): Chamber pressure
        T_chamber (K): Chamber temperature
        w_throat (m): Throat width
        h_throat (m): Throat heigh (or channel depth)
        throat_roc (m): Throat radius of curvature
        AR_exit (-): Area ratio of nozzle exit area divided by throat area
        p_back (Pa): Back pressure
        divergence_half_angle (rad): Divergence half angle of nozzle 
        fp (FluidProperties): Object to access fluid properties
        is_cold_flow (bool): Reynolds number is adjusted depending on whether the chamber is heated or cooled

    Raises:
        ValueError: Is raised for hot flow, since no verification is done yet on that equation
    """
    
    # Get the (assumed to be) constant fluid properties
    gamma = fp.get_specific_heat_ratio(T=T_chamber,p=p_chamber) # [-] Specific heat ratio
    R = fp.get_specific_gas_constant() # [J/kg] Specific gas constant
    # Report calculated values for verification and comparison purposes
    print("Gamma: {:1.4f}".format(gamma))
    print("R: {:3.2f} J/kg\n".format(R))

    # Calculate basic peformance parameters
    A_throat = w_throat * h_throat # [m] Throat area
    
    ## IDEAL PERFORMANCE
    # First get ideal performance, and check if the nozzle is properly expanded.
    ep = IRT.get_engine_performance(p_chamber=p_chamber,T_chamber=T_chamber,A_throat=A_throat, AR_exit=AR_exit, p_back=p_back, gamma=gamma, R=R)

    # Report ideal performance
    print("Thrust: {:.2f} mN".format(ep['thrust']*1e3))
    print("Isp_ideal: {:.1f} s".format(ep['thrust']/ep['m_dot']/9.80655))
    print("Mass flow: {:.3f} mg/s".format(ep['m_dot']*1e6))

    m_dot_ideal = ep['m_dot'] # [kg/s] Ideal mass flow
    F_ideal = ep['thrust'] # [N] Ideal thrust

    ## CALCULATING THE CORRECTION FACTORS

    # Calculate the divergence loss and report it
    CF_divergence_loss = divergence_loss_conical_2D(alpha=divergence_half_angle)
    print("\n -- DIVERGENCE LOSS for {:2.2f} deg divergence half-angle".format(math.degrees(divergence_half_angle)))
    print("  Divergence loss (2D concical): {:.5f} ".format(CF_divergence_loss))

    # Calculate the viscous loss

    # To determine the Reynolds number at the throat, the hydraulic diameter at the throat and nozzle conditions must be determined
    # Get hydraulic diameter of the nozzle from the wetted perimeter and nozzle area
    wetted_perimeter_throat = 2 * (w_throat+h_throat) # [m] Wetted perimeter throat
    Dh_throat = hydraulic_diameter(A=A_throat, wetted_perimeter=wetted_perimeter_throat) # [m] Hydraulic diameter at throat
    p_throat = p_chamber/IRT.pressure_ratio(M=1, gamma=gamma) # [Pa] pressure in throat
    T_throat = T_chamber/IRT.temperature_ratio(M=1, gamma=gamma) # [K] Temperature in throat
    viscosity_throat = fp.get_viscosity(T=T_throat,p=p_throat)
    # Throat reynolds based on ideal mass flow?
    Re_throat = reynolds(m_dot=m_dot_ideal,A=A_throat,D_hydraulic=Dh_throat,viscosity=viscosity_throat)
    if is_cold_flow:
        Re_throat_wall = Reynolds_throat_wall_cold(reynolds_throat=Re_throat)
    else:
        Re_throat_wall = Reynolds_throat_wall_hot(reynolds_throat=Re_throat)
    print("\n-- THROAT CONDITIONS --")
    print("  p = {:2.4f} bar,     T = {:4.2f} K".format(p_throat*1e-5, T_throat))
    print("  mu = {:2.4f} [microPa*s]  Dh = {:3.4f} [microm]".format(viscosity_throat*1e6,Dh_throat*1e6))
    print(" Reynolds: {:6.6f} ".format(Re_throat))

    CF_viscous_loss = viscous_loss(area_ratio=AR_exit,reynolds_throat_wall=Re_throat_wall)
    print(" CF_viscous_loss: {:1.5f}".format(CF_viscous_loss))

    # Calculating throat boundary layer loss, which causes a reduction in effective throat area/mass flow
    Cd_throat_boundary_loss = throat_boundary_loss(gamma=gamma,reynolds_throat=Re_throat,throat_radius=0.5*Dh_throat,throat_roc=throat_roc)
    print("\n-- DISCHARGE FACTOR --")
    print("  Throat boundary layer: {:1.4f}".format(Cd_throat_boundary_loss))


    ## APPLYING THE CORRECTION FACTORS
    # Now all these loss factors must be combined into a new "real" thrust 
    # The divergence loss only applies to the jet/momentum thrust and not the pressure, so jet thrust is needed
    # This is equal to the exit velocity times corrected mass flow. The returned exit velocity does not include pressure terms!

    # First we must know the corrected mass flow
    m_dot_real = ep['m_dot'] * Cd_throat_boundary_loss # [kg/s]
    # Secondly, we must know the pressure thrust to add to the jet thrust again
    F_pressure = IRT.pressure_thrust(p_chamber=p_chamber, p_back=p_back, A_throat=A_throat,AR=AR_exit, gamma=gamma)
    F_divergence = m_dot_real * ep['u_exit'] * CF_divergence_loss +  F_pressure # [N] Thrust decreased by divergence loss, pressure term must be added again, since divergence only applies to jet thrust
    # This jet thrust is then again corrected by viscous losses, which are subtracted from the current CF
    CF_jet_divergence = F_divergence / (p_chamber * A_throat) # [-] Thrust coefficient after taking into account discharge factor and divergence loss
    CF_real_final = CF_jet_divergence - CF_viscous_loss # [-] The final thrust coefficient, also taking into account viscous loss
    F_real = CF_real_final * p_chamber * A_throat # [N] Real thrust, after taking into account of all the three proposed correction factors

    # Report "real" results
    print("\n === CORRECTED PERFORMANCE PARAMETERS === ")
    print("  Real mass flow: {:3.4f} mg/s".format(m_dot_real*1e6))
    print("  CF with divergence loss {:1.5f}".format(CF_jet_divergence))
    print("  Real CF: {:1.5f}".format(CF_real_final))
    print("  Real thrust: {:2.4f} mN".format(F_real*1e3))
    print("  Real Isp: {:3.2f}".format(F_real/m_dot_real/9.80655))

    return {    'm_dot_real': m_dot_real,
                'm_dot_ideal': ep['m_dot'],
                'F_real': CF_real_final}


