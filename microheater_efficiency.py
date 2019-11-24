import numpy as np
import matplotlib.pyplot as plt

from constants import stefan_boltzmann
import thermo.water as water

def calc_P_delta_h(T_in,T_out,m_dot,p_chamber):
    """Calculate power required to raise temperature of flow with mass flow m_dot
    
    Arguments:
        T_in {K} -- Inlet temperature
        T_out {K} -- Outlet temperature
        m_dot {kg/s} -- Mass flow
        p_chamber {Pa} -- Chamber pressure
    
    Returns:
        {W} - Required power
    """
    # Check if outlet state is gaseous
    outlet_phase = water.get_phase(T=T_out,p=p_chamber)
    if not (outlet_phase == "gas" or outlet_phase == "supercritical_gas"):
        print(water.get_phase(T=T_out,p=p_chamber))
        raise ValueError("Assumed outlet temperature is not not in gas phase")
    
    Delta_h = water.get_enthalpy(T=T_out, p=p_chamber) - water.get_enthalpy(T=T_in, p=p_chamber)

    return m_dot*Delta_h

def calc_P_rad(T_mh_range,emissivity,A_mh):
    """Calculation radation loss from micro-heater surface
    
    Arguments:
        T_mh_range {K} -- Micro-heater temperature -- Needs to a python object, because ^4 makes it very large
        emissivity {-} -- Emissivity of material
        A_mh {m^2} -- Micro-heater area
    
    Returns:
        W -- Power radiated from micro-heater surface
    """
    
    #check if it is a numpy array with a python objectas dtype
    if not isinstance(T_mh_range, np.ndarray):
        raise TypeError("T_mh_range must a numpy array (even if 1 value)")

    if not (T_mh_range.dtype == object):
        raise TypeError("T_mh_range must be a numpy array with dtype=object")

    return stefan_boltzmann*emissivity*A_mh*T_mh_range**4

def calc_P_cond(A_chip, T2, T1, L , thermal_conductivity):
    """Calculation conduction loss 
    
    Arguments:
        A_chip {m^2} -- Chip surface area assumed to be in contact with other surface
        T2 {K} -- Temperature of chip
        T1 {K} -- Temperature of other object on other side of chip
        L {m} -- Thickness of material
        k {W/mK} -- Thermal conductivity of material
    """

    return -thermal_conductivity*A_chip*(T2-T1)/L

def calc_microheater_efficiency(T_mh_range,emissivity,A_mh,T_in,T_margin,m_dot,p_chamber):
    """Calculates micro-heater efficiency based on assumption that:

    * Propellant is water
    * Radiation is only form of heat loss
    * Outlet temperature is raised to a certain margin below micro-heater temperature
    * Constant chamber pressure 
    
    Arguments:
        T_mh_range {np.array(dtype=object) [K]} -- Numpy array for which efficiencies will be calculated NOTE: must be converted to dtype=object because it will be ^4 and huge.
        emissivity {-} -- Emmissivity from Boltzmann's law
        A_mh {m^2} -- Assumed micro-heater area
        T_in {K} -- Assumed inlet temperature
        T_margin {K} -- Assumed margin below micro-heater temperature to which outlet temperature is raised
        m_dot {kg/s} -- Assumed mass flow of micro-heater
        p_chamber {Pa} -- Assumed chamber pressure (default: 1 bar)
    """
    # T_mh_range must be numpy array
    if not isinstance(T_mh_range, np.ndarray):
        raise TypeError("T_mh_range must a numpy array (even if 1 value)")

    if not (T_mh_range.dtype == object):
        raise TypeError("T_mh_range must be a numpy array with dtype=object")


    T_out = T_mh_range - T_margin # [K] Assume exit temperature is below micro-heater temperature by T_margin
   
    
    # If in range, calculate power provided to fluid
    P_Delta_h = np.zeros_like(T_out) # Make array of known size to store enthalpy difference for iteration
    for i, T_now in enumerate(T_out):
         # Get power required to raise temperature of flwo
        P_Delta_h[i] = calc_P_delta_h(T_in=T_in, T_out=T_now, m_dot=m_dot, p_chamber=p_chamber) # pylint: disable=unsupported-assignment-operation
    
    # Then calculate power lost through radiation
    P_rad = calc_P_rad(T_mh_range=T_mh_range,emissivity=emissivity,A_mh=A_mh)
    # Return efficiency
    return P_Delta_h / (P_Delta_h + P_rad)

def calc_microheater_efficiency_two(T_mh_range,emissivity,A_mh,T_in,T_margin,m_dot,p_chamber,thermal_conductivity,T_room,L_glass):
    """Calculates micro-heater efficiency based on assumption that:

    * Propellant is water
    * Radiation is only form of heat loss
    * Outlet temperature is raised to a certain margin below micro-heater temperature
    * Constant chamber pressure 
    * Heat is lost through conduction through glass
    
    Arguments:
        T_mh_range {np.array(dtype=object) [K]} -- Numpy array for which efficiencies will be calculated NOTE: must be converted to dtype=object because it will be ^4 and huge.
        emissivity {-} -- Emmissivity from Boltzmann's law
        A_mh {m^2} -- Assumed micro-heater area
        T_in {K} -- Assumed inlet temperature
        T_margin {K} -- Assumed margin below micro-heater temperature to which outlet temperature is raised
        m_dot {kg/s} -- Assumed mass flow of micro-heater
        p_chamber {Pa} -- Assumed chamber pressure (default: 1 bar)
    """
    # T_mh_range must be numpy array
    if not isinstance(T_mh_range, np.ndarray):
        raise TypeError("T_mh_range must a numpy array (even if 1 value)")

    if not (T_mh_range.dtype == object):
        raise TypeError("T_mh_range must be a numpy array with dtype=object")


    T_out = T_mh_range - T_margin # [K] Assume exit temperature is below micro-heater temperature by T_margin
   
    
    # If in range, calculate power provided to fluid
    P_Delta_h = np.zeros_like(T_out) # Make array of known size to store enthalpy difference for iteration
    for i, T_now in enumerate(T_out):
         # Get power required to raise temperature of flwo
        P_Delta_h[i] = calc_P_delta_h(T_in=T_in, T_out=T_now, m_dot=m_dot, p_chamber=p_chamber) # pylint: disable=unsupported-assignment-operation
    
    # Then calculate power lost through radiation
    P_rad = calc_P_rad(T_mh_range=T_mh_range,emissivity=emissivity,A_mh=A_mh)
    # Also add power lost through conduction
    P_cond = -calc_P_cond(A_chip=A_mh, T2=T_mh_range,T1=T_room,L=L_glass,thermal_conductivity=thermal_conductivity)

    # Add losses
    P_loss = P_cond+P_rad
    # Return efficiency
    return P_Delta_h / (P_Delta_h + P_loss)   
    

def run():
    emissivity = 1 # [-] Emmisivity of micro-heater, just 1 for illustrative purposes
    T_in = 297 # [K] Room temperature at inlet
    T_margin = 0 # [K] No margin
    m_dot = 0.5e-6 # [kg/s] 1 mg/s mass flow assumed
    A_mh = 6e-3*10e-3 # [m^2]

    p_chamber = 1e5 # [Pa] 1 bar chamber pressure assumed

    T_mh_range_1 = np.arange(400,1201, 10,dtype=object)
    efficiency_1 = calc_microheater_efficiency(T_mh_range=T_mh_range_1, emissivity=emissivity, A_mh=A_mh, T_in=T_in,T_margin=T_margin,m_dot=m_dot,p_chamber=p_chamber)
    # Again but with a margin
    T_margin_2 = 25 # [K]
    efficiency_2 = calc_microheater_efficiency(T_mh_range=T_mh_range_1, emissivity=emissivity, A_mh=A_mh, T_in=T_in,T_margin=T_margin_2,m_dot=m_dot,p_chamber=p_chamber)
    # Margin of 25 K now
    T_margin_3 = 100
    T_mh_range_3 = np.arange(480,1201, 10,dtype=object)
    efficiency_3 = calc_microheater_efficiency(T_mh_range=T_mh_range_3, emissivity=emissivity, A_mh=A_mh, T_in=T_in,T_margin=T_margin_3,m_dot=m_dot,p_chamber=p_chamber)
    # Margin of 50 K now
    T_margin_4 = 250
    T_mh_range_4 = np.arange(630,1201, 10,dtype=object)
    efficiency_4 = calc_microheater_efficiency(T_mh_range=T_mh_range_4, emissivity=emissivity, A_mh=A_mh, T_in=T_in,T_margin=T_margin_4,m_dot=m_dot,p_chamber=p_chamber)
    #Plot this shit
    plt.close()
    plt.plot(T_mh_range_1,efficiency_1, label="0 K margin")
    plt.plot(T_mh_range_1,efficiency_2, label="25 K margin")
    plt.plot(T_mh_range_3,efficiency_3, label="100 K margin")
    plt.plot(T_mh_range_4,efficiency_4, label="250 K margin")
    plt.xlabel("Micro-heater temperature $T_{mh} [K]$")
    plt.ylabel("Micro-heater efficiency $\\mu_{mh} [-]$")
    plt.title("Micro-heater efficiency calculation\n based only on radiative losses and incomplete propellant heating")
    plt.legend(title="$A_{{chip}}={Amh}$mm$^2$ \n $\\varepsilon={em}$\n$\dot{{m}}={m_dot}$ mg/s\n$p_c={p}$ bar\n$T_{{mh}} - T_{{out}}:$".format(Amh=A_mh*1e6, em=emissivity,m_dot=m_dot*1e6,p=p_chamber*1e-5))
    plt.grid()
    plt.show()


def run2():
    #Calcuate heating efficiency of heaters in literature
    m_dot = 0.83e-6 # [kg/s] mass flow
    T_in = 24+273.15 # [K] Inlet temperature
    T_out = 426.65 # [K] Outlet temperature
    p = 5.15e5 # [Pa] Pressure
    P_total = 8.19 # [W] Total electrical power

    # Specific enthalpy at inlet and outlet
    h_in = water.get_enthalpy(T=T_in, p=p) # [J/kg] Inlet
    h_out = water.get_enthalpy(T=T_out, p=p) # [J/kg] Outlet
    
    P_delta_h = m_dot*(h_out-h_in) # [W] Power raising enthalpy
    efficiency = P_delta_h/P_total # [-]

    print("Exit phase: {}".format(water.get_phase(T=T_out,p=p)))

    print("P_delta_h: {:1.2f} W".format(P_delta_h))
    print("Micro-heater efficiency: {:1.2f} ".format(efficiency))

    print("T_in = {:3.0f} K \t\t T_out =  {:3.0f} K".format(T_in,T_out))
    print("Mass flow = {:2.2f} mg/s".format(m_dot*1e6))
    
def run3():
    emissivity = 1 # [-] Emmisivity of micro-heater, just 1 for illustrative purposes
    T_in = 297 # [K] Room temperature at inlet
    T_margin = 0 # [K] No margin
    m_dot = 0.5e-6 # [kg/s] 1 mg/s mass flow assumed
    A_mh = 6e-3*10e-3 # [m^2]
    thermal_conductivity = 1.5 # [W/mK] Thermal conductivity of glass
    t_glass = 1e-3 # Thickness of glass assumed to be 1mm
    T_room = 300 # Room temperature (given at 24C but 300K is close enough)
    p_chamber = 1e5 # [Pa] 1 bar chamber pressure assumed

    T_mh_range_1 = np.arange(400,1201, 10,dtype=object)
    efficiency_1 = calc_microheater_efficiency_two(T_mh_range=T_mh_range_1, emissivity=emissivity, A_mh=A_mh, T_in=T_in,T_margin=T_margin,m_dot=m_dot,p_chamber=p_chamber,thermal_conductivity=thermal_conductivity,T_room=T_room,L_glass=t_glass)
    # Again but with a margin
    T_margin_2 = 0 # [K]
    efficiency_2 = calc_microheater_efficiency(T_mh_range=T_mh_range_1, emissivity=emissivity, A_mh=A_mh, T_in=T_in,T_margin=T_margin_2,m_dot=m_dot,p_chamber=p_chamber)
    # Margin of 25 K now
   
    #Plot this shit
    plt.close()
    plt.plot(T_mh_range_1,efficiency_2, label="Radiation only")
    plt.plot(T_mh_range_1,efficiency_1, label="Radiation and conduction")
    plt.xlabel("Micro-heater temperature $T_{{mh}} [K]$")
    plt.ylabel("Micro-heater efficiency $\\mu_{{mh}} [-]$")
    plt.title("Micro-heater efficiency calculation\n based on radiative and/or conductive losses (0K margin)")
    plt.legend(title="$A_{{chip}}={Amh}$mm$^2$ \n $\\varepsilon={em}$\n$\dot{{m}}={m_dot}$ mg/s\n$p_c={p}$ bar".format(Amh=A_mh*1e6, em=emissivity,m_dot=m_dot*1e6,p=p_chamber*1e-5))
    plt.grid()
    plt.show()

if __name__ == "__main__":
    print("Running microheater efficiency script")
    run3()