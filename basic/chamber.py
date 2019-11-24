import thermo.water as water
from constants import stefan_boltzmann

def ideal_enthalpy_change(T_inlet, p_inlet, T_outlet, p_outlet):
    """Returns specific enthalpy change based on simple chamber inlet and outlet conditions.
    This should give the power the micro-heater must transfer in ideal conditions with no heat losses.
    In addition returns a warning if the final state is not gaseous.
    
    Arguments:
        T_inlet {K} -- Inlet temperature
        p_inlet {Pa} -- Inlet pressure
        T_outlet {K} -- Outlet temperature
        p_outlet {Pa} -- Outlet pressure
    
    Returns:
        delta_h {J/(kg*K)} -- 
    """

    h_inlet = water.get_enthalpy(T=T_inlet, p=p_inlet)
    h_outlet = water.get_enthalpy(T=T_outlet, p=p_outlet)

    outlet_phase = water.get_phase(T=T_outlet,p=p_outlet)
    if not outlet_phase == 'gas':
        print("Warning: Phase at chamber exit is not gaseous but {}".format(outlet_phase))
        
    return h_outlet-h_inlet

def ideal_power_consumption(mass_flow, T_inlet, p_inlet, T_outlet, p_outlet):
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

    return mass_flow * ideal_enthalpy_change(T_inlet=T_inlet, p_inlet=p_inlet, T_outlet=T_outlet, p_outlet=p_outlet)

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

    return emmisivity * A * stefan_boltzmann * T**4

    