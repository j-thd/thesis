# File for basic IRT functions
# For source of equation check appendix A of literature study
from physical_constants import g0
from sys import float_info
import math
from scipy import optimize
from thermo.prop import FluidProperties


def vdk(gamma):
    """Calculates the vandenkerckhove function
    
    Arguments:
        gamma {-} -- Specific heat ratio
    
    Returns:
        Gamma(gamma) {-}
    """
    return math.sqrt(gamma)*(2/(gamma+1))**((gamma+1)/(2*(gamma-1)))

def mass_flow(p_chamber, A_throat, R, T_chamber, gamma):
    """Calculates mass flow based on IRT
    
    Arguments:
        p_chamber {Pa} -- Chamber pressure
        A_throat {m^2} -- Throat area of nozzle
        R {J/(kg*K)} -- Specific gas constant of propellant
        T_chamber {K} -- Chamber temperature
        gamma {-} -- Specific heat ratio
    
    Returns:
        Mass flow {kg/s} -- 
    """
    return p_chamber*A_throat/math.sqrt(R*T_chamber)*vdk(gamma)

def area_ratio(M, gamma):
    """Returns the area ratio for a specific Mach number
    
    Arguments:
        M {-} -- Mach number
        gamma {-} -- Specific heat ratio
    
    Returns:
        Area ratio {-} -- Local area divivded by throat area
    """

    exponent = (gamma+1)/(2*(gamma-1)) # Exponent in formula
    base = 2/(gamma+1)*( 1 + (gamma-1)/2*M**2)
    return 1/M*base**exponent

def temperature_ratio(M, gamma):
    """Returns the temperature ratio w.r.t to total/chamber temperature
    
    Arguments:
        M {-} -- Mach number
        gamma {-} -- Specific heat ratio
    
    Returns:
        Temperature Ratio {-} -- Total/chamber temperature divided by local temperature
    """
    return 1 + (gamma-1)/2*M**2

def pressure_ratio(M, gamma):
    """Returns the pressure ratio w.r.t total/chamber pressure by calculating the temperature ratio by raising it to the right power according to isentropic relations
    
    Arguments:
        M {-} -- Mach number
        gamma {-} -- Specific heat ratio
    
    Returns:
        Pressure ratio {-} -- Total/chamber pressure divided by local pressure
    """
    return temperature_ratio(M=M,gamma=gamma)**(gamma/(gamma-1))

def density_ratio(M,gamma):
    """Returns the density ratio w.r.t total/chamber density by calculating the temperature ratio by raising it to the right power according to isentropic relations
    
    Arguments:
        M {-} -- Mach number
        gamma {-} -- Specific heat ratio
    
    Returns:
        density ratio {-} -- Total/chamber density divided by local density
    """

    return temperature_ratio(M=M,gamma=gamma)**(1/(gamma-1))

def is_throat_sonic(p_chamber, p_back, exit_area_ratio, gamma):
    """Returns whether the throat has reached sonic conditions at all, to verify that chocked flow occurs in the nozzle
    
    Arguments:
        p_chamber {Pa} -- Chamber/total pressure
        p_exit {Pa} -- Exit pressure
        exit_area_ratio {-} -- Exit area divided by throat area
        gamma {-} -- Specific heat ratio
    
    Returns:
        whether flow is chocked/throat is sonic (True/False) -- [description]
    """
    # Get the subsonic Mach number at the exit to determine pressure ratio at exit
    M_exit = Mach_from_area_ratio(AR=exit_area_ratio, gamma=gamma,supersonic=False)
    pressure_ratio_exit = pressure_ratio(M=M_exit,gamma=gamma)

    assert p_back >=0 # No negative pressures
    assert p_chamber > 0 # Everything is useless if this is zero anyway.
    if p_back == 0:
        return True
    return (p_chamber/p_back) > pressure_ratio_exit

def exit_pressure_ratio_shockwave(exit_area_ratio, gamma):
    """Returns the pressure ratio for which a shock is present exactly at the exit of the throat.
    Useful as this is the limiting pressure under which the divergent part of the nozzle should be fully supersonic
    
    Arguments:
        exit_area_ratio {-} -- Exit divided by throat area
        gamma {-} -- Specific heat ratio
    
    Returns:
        pressure ratio{-} -- Chamber/total pressure divided by pressure for which shock occurs precisely at exit
    """
    M_exit = Mach_from_area_ratio(AR=exit_area_ratio, gamma=gamma)
    PR_exit = pressure_ratio(M=M_exit, gamma=gamma) # Total pressure divived by pressure before shockwave
    PR_shockwave = pressure_ratio_shockwave(M=M_exit, gamma=gamma) # Pressure after shockwave divided before shockwave
    # This means PR_shockwave is divived instead of multiplied to get PR w.r.t to total pressure
    return PR_exit/PR_shockwave

def pressure_ratio_shockwave(M, gamma):
    """Returns the pressure drop over a shockwave (Eq. 8.65 - Anderson2006)
    
    Arguments:
        M {-} -- Mach number before shockwave
        gamma {-} -- Specific heat ratio
    
    Returns:
        pressure ratio {-} -- Pressure AFTER shockwave divided by pressure BEFORE shockwave (NO TOTAL PRESSURE INVOLVED, CONTRARY WITH OTHER PRESSURE RATIOS IN CODE)
    """

    return 1 + (2*gamma)/(gamma+1)*(M**2-1)


def nozzle_status(p_chamber, p_back, AR_exit, gamma, isentropic_tolerance=0.01):
    """Returns the status of the nozzle depending on exit pressure to inform designer of (im)proper pressure ratios 
    Isentropic expansion is tested within a certain tolerance margin

    ADDITION: it now also check if no shock occurs in the divergent part of the nozzle. 
    This value should lie somewhere between the pressure where the nozzle turns sonic,
    and where it is overexpanded (but expanded outside the nozzle through shockwaves)
    
    Arguments:
        p_chamber {Pa} -- Chamber/total pressure
        p_exit {Pa} -- Exit pressure
        exit_area_ratio {-} -- Exit area divided by throat area
        gamma {-} -- Specific heat ratio
    
    Keyword Arguments:
        isentropic_tolerance {-} -- How much exit pressure can vary from isentropic expansion in order to report it as isentropically expanded

    Returns:
        string -- ('subsonic'/'underexpanded'/'overexpanded'/'isentropically expanded')
    """
    # First check if throat is sonic, if not return subsonic flow
    if not is_throat_sonic(p_chamber=p_chamber, p_back=p_back, exit_area_ratio=AR_exit, gamma=gamma):
        return "subsonic"
    # Else, it is calculated what the exit pressure should be for isentropic expansion towards the nozzle
    else:
        # Find the Mach number for isentropic expansion to exit to find the pressure ratio
        M_exit= Mach_from_area_ratio(AR=AR_exit, gamma=gamma)
        PR_exit = pressure_ratio(M=M_exit, gamma=gamma)
        p_exit_isentropic = p_chamber / PR_exit # This is what the exit pressure should be for isentropic expansion
        #Find the pressure ratio, for which a normal shockwave is at the nozzle exit
        PR_exit_shockwave = exit_pressure_ratio_shockwave(exit_area_ratio=AR_exit, gamma=gamma)
        p_exit_shockwave = p_chamber / PR_exit_shockwave
        # Low and high bound to consider flow isentropic
        p_exit_isentropic_low = p_exit_isentropic*(1-isentropic_tolerance)
        p_exit_isentropic_high = p_exit_isentropic*(1+isentropic_tolerance)
        # Now assess whether it over-, under- or isentropically expanded
        # If the exit pressure is higher than ambient pressure it is underexpanded
        if p_back >= p_exit_shockwave:
            return "shock in nozzle"
        elif p_back > p_exit_isentropic_high and p_back < p_exit_shockwave:
            return 'overexpanded'
        elif p_back <= p_exit_isentropic_high and p_back >= p_exit_isentropic_low:
            return f'isentropically expanded (margin: {isentropic_tolerance})'.format()
        elif p_back >=0 and p_back < p_exit_isentropic_low:
            return 'underexpanded'
        
def Mach_from_area_ratio(AR, gamma, supersonic=True,precision=float_info.epsilon):
    """Returns the Mach number for a specific area ratio. Defaults to supersonic value.
    For computational simplicity raises an error for really high area ratios, so that the root-finding
    can be simplified to an interval search, and to avoid problem with asymptote as Mach number goes to
    zero.
    
    Arguments:
        AR {-} -- Local area divided by throat area
        gamma {-} -- Specific heat ratio
    
    Keyword Arguments:
        supersonic {bool} -- Returns supersonic Mach number by default, and subsonic when set to False (default: {True})
    """
    # Check if area_ratio is not too ridiculously high in superonsic case
    if AR > 1e7 and supersonic == True:
        raise ValueError("area_ratio too high. Check function documentation to solve problem.")
    # First the bounds of the interval search must be determined
    if supersonic:
        bounds = (1,100)
    else:
        bounds = (1e-7,1)

    # Find root of area_ratio function subtracted by desired area ratio
    f = lambda M: area_ratio(M=M, gamma=gamma)-AR
    return optimize.bisect(f=f,a=bounds[0],b=bounds[1])

def exit_velocity(AR_exit, T_chamber, R, gamma):
    """Returns exit velocity, assuming flow in nozzle is supersonic.
    Must be verified by nozzle status function!
    
    Arguments:
        AR_exit {-} -- Exit area divided throat area
        T_chamber {K} -- Total/chamber temperature
        R {J/kg*K} -- Specific gas constant
        gamma {-} -- Specific heat ratio
    
    Returns:
        Exit velocity (for supersonic nozzle flow) {m/s}
    """
    # By determining the exit Mach number the pressure/temperature ratio can be found
    M_exit = Mach_from_area_ratio(AR=AR_exit, gamma=gamma)
    # Temperature is slightly easier to implement, more readable
    TR_exit = temperature_ratio(M=M_exit, gamma=gamma)

    return (2*R*gamma*T_chamber/(gamma-1)*(1-1/TR_exit))**0.5 

def throat_velocity(T_chamber, R, gamma):
    """Return throat velocity by returning speed of sound in throat

    Args:
        T_chamber (K): Total/chamber temperature
        R (J/(kg*K)): Specific gas constant
        gamma (-): Specific heat ratio

    Returns:
        u_throat (m/s): Throat velocity (assuming it is speed of sound)
    """
    TR_throat = temperature_ratio(M=1,gamma=gamma) # [-] Temperature ratio in the throat
    T_throat = T_chamber / TR_throat # [K] Throat temperature

    return (gamma*R*T_throat)**0.5 # [m/s] Speed of sound calculation in throat

def pressure_thrust(p_chamber, p_back, A_throat, AR, gamma):
    """Returns the pressure thrust according to IRT

    Args:
        p_chamber (Pa): Chamber/total pressure
        p_back (Pa): Back pressure
        A_throat (m^2): Throat area
        AR (-): Exit area ratio
        gamma (-): Ratio of specific heats

    Returns:
        Pressure thrust (N): The pressure differential at the exit, multiplied by the exit area
    """
    M_exit = Mach_from_area_ratio(AR=AR,gamma=gamma)
    PR_exit = pressure_ratio(M=M_exit,gamma=gamma)
    p_exit = p_chamber/PR_exit
    A_exit = A_throat*AR

    return (p_exit-p_back)*A_exit

def effective_ISP(thrust, m_dot):
    """Returns specific impulse in seconds, based on standard gravity

    Args:
        thrust (N): Thrust
        m_dot (kg/s): Mass flow

    Returns:
        I_sp (s): Effective specific impulse based on standard gravity
    """
    return thrust/(m_dot*g0)

def thrust(p_chamber, T_chamber, A_throat, AR_exit, p_back, gamma, R):
    """Returns the thrust according to IRT
    
    Arguments:
        p_chamber {Pa} -- Chamber/total pressure
        T_chamber {K} -- Chamber/total temperature
        A_throat {m^2} -- Throat area of nozzle
        AR {-} -- Exit area divided by throat area
        p_back {Pa} -- Pressure of environment
        gamma {-} -- Specific heat ratio
        R {-} -- Gas constant of propellant
    
    Returns:
        thrust F {N} -- Returns the thrust according to Ideal Rocket Theory
    """
    m_dot = mass_flow(p_chamber=p_chamber, A_throat=A_throat, R=R, T_chamber=T_chamber, gamma=gamma)
    u_exit = exit_velocity(AR_exit=AR_exit,T_chamber=T_chamber,R=R,gamma=gamma)
    # Find the exit pressure to determine the pressure force
    F_pressure = pressure_thrust(p_chamber=p_chamber,p_back=p_back,A_throat=A_throat, AR=AR_exit, gamma=gamma)


    # Finally, some checks are done and errors are thrown is the nozzle is not supersonic until at least the exit
    NS = nozzle_status(p_chamber=p_chamber, p_back=p_back, AR_exit=AR_exit, gamma=gamma)
    if(NS == "subsonic"):
        print("Calculations for thrust assume that M=1 at the throat. The throat is subsonic, so the assumption is invalid.")
    if(NS == "shock in nozzle"):
        print("The throat reaches sonic conditions, but pressure is so low that normal shock occurs in nozzle. Therefore M<1 at the exit and calculations assuming supersonic conditions there are invalid.")
    
    # print("Jet thrust: {:.2f} mN".format(m_dot*u_exit*1e3))
    # print("Jet Isp {:3.0f} s".format(u_exit/9.81))
    # Otherwise, it's fine and returns the thrust
    return m_dot*u_exit + F_pressure

def get_engine_performance(p_chamber, T_chamber, A_throat, AR_exit, p_back, gamma, R):
    """Returns a dictionary of most relevant engine performance. Raises an error if supersonic assumptions are broken
    
    Arguments:
        p_chamber {Pa} -- Chamber/total pressure
        T_chamber {K} -- Chamber/total temperature
        A_throat {m^2} -- Throat area
        AR_exit {-} -- Exit area divided throat area
        p_back {Pa} -- Ambient/environmental pressure, to determine pressure force of thrust
        gamma {-} -- Specific heat ratio
        R {-} -- Specific gas constant
    
    Returns:
        dictionary with following parameters: # A dictionary allows for easy extension of variable without breaking previous code
            thrust {N} -- Total thrust
            m_dot {kg/s} - Mass flow of engine
            u_exit {u_exit} - Exit velocity (NOT effective exit velocity, so no pressure terms included)
            nozzle_status - Status of nozzle (just an fyi function)
    """
    # First check if assumptions about sonic conditions in throat and supersonic conditions at exit hold
    NS = nozzle_status(p_chamber=p_chamber, p_back=p_back, AR_exit=AR_exit, gamma=gamma)
    if(NS == "subsonic"):
        print("Calculations for thrust assume that M=1 at the throat. The throat is subsonic, so the assumption is invalid.")
    if(NS == "shock in nozzle"):
        print("The throat reaches sonic conditions, but pressure is so low that normal shock occurs in nozzle. Therefore M<1 at the exit and calculations assuming supersonic conditions there are invalid.")

    # If this is not a problem, it is fine to return other values, as they will be correct
    m_dot = mass_flow(p_chamber=p_chamber, A_throat=A_throat,R=R, T_chamber=T_chamber, gamma=gamma) # [kg/s]
    u_exit = exit_velocity(AR_exit=AR_exit,T_chamber=T_chamber, R=R, gamma=gamma) # [m/s]
    F = thrust(p_chamber=p_chamber, T_chamber=T_chamber, A_throat=A_throat, AR_exit=AR_exit, p_back=p_back, gamma=gamma, R=R) # [N]
    Isp = effective_ISP(thrust=F, m_dot=m_dot) # [s]
    
    # Add information about throat for correction purposes (calculating Reynolds numbers etc..)
    T_throat = T_chamber / temperature_ratio(M=1,gamma=gamma)  # [K] Temperature in nozzle
    p_throat = p_chamber / pressure_ratio(M=1,gamma=gamma) # [Pa] Pressure in nozzle
    u_throat = throat_velocity(T_chamber=T_chamber,R=R,gamma=gamma) # [m/s] Velocity in throat


    # Add information about exit to check for condensation in nozzle exit
    M_exit = Mach_from_area_ratio(AR=AR_exit, gamma=gamma) # [-] Mach number at the xit
    T_exit = T_chamber / temperature_ratio(M=M_exit,gamma=gamma) # [K] Temperature at exit
    p_exit = p_chamber / pressure_ratio(M=M_exit, gamma=gamma) # [Pa] Pressure at exit

    return {'thrust': F,
            'Isp': Isp,
            'm_dot': m_dot,
            'u_exit': u_exit,
            'nozzle_status': NS,
            'T_throat': T_throat,
            'p_throat': p_throat,
            'u_throat': u_throat,
            'T_exit': T_exit,
            'p_exit': p_exit}

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
    A_high = 1e-4 # [K]

    R = fp.get_specific_gas_constant() # [J/(kg*K)]
    gamma = fp.get_specific_heat_ratio(T_chamber,p_chamber) # [-] Specific heat ratio
    # If x is zero, then the desired thrust is found. Gamma is changed depending on temperature as well
    x = lambda A_t: F_desired - get_engine_performance(p_chamber=p_chamber, T_chamber=T_chamber, A_throat=A_t, \
                        AR_exit=AR_exit,p_back=p_back, gamma=gamma,R=R)['thrust']

    root_result = optimize.root_scalar(x, bracket=[A_low,A_high], xtol=1e-12
    )

    if root_result.converged:
        A_throat = root_result.root # [m^2]
    else:
        raise Exception("No solution found for given temperature")

    # Now the chamber temperature is known, return the unknown parameters

    ep = get_engine_performance(p_chamber=p_chamber, T_chamber=T_chamber, A_throat=A_throat, \
                        AR_exit=AR_exit,p_back=p_back, gamma=fp.get_specific_heat_ratio(T_chamber,p_chamber),R=R)
    # Add throat area to dictionary
    ep['A_throat'] = A_throat # This is usually an input for get_engine_performance, but now added to make it one complete solution
    #ep['Isp'] = IRT.effective_ISP(thrust=F_desired,m_dot=ep['m_dot']) # [s] Effective specific impulse

    return ep