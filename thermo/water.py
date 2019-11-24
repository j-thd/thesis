from CoolProp.CoolProp import PropsSI, PhaseSI

def get_enthalpy(T, p):
    """Returns SPECIFIC enthalpy for water through CoolProp library
    
    Arguments:
        T {K} -- Temperature
        P {Pa} -- Pressure
    
    Returns:
        h {J/kg}
    """
    # Documentation says "HEOS::Water" is the most accurate, but it is also slow.
    # "IF97::Water" is faster but slightly less accurate
    h = PropsSI('H','T',T,'P',p,"HEOS::Water")
    return h

def get_phase(T,p):
    """Returns the phase of water through CoolProp library

    Arguments:
        T {K} -- Temperature
        P {Pa} -- Pressure
    
    Returns:
        Returns phase according to CoolProp
        ('liquid'/'gas'/'supercritical_gas'/'supercritical'/'supercritical_liquid')
    """ 

    return PhaseSI('T',T,'P',p,"HEOS::Water")

# def get_temperature(p,h):
#     """Return temperature from pressure and enthalpy through CoolProp libray
    
#     Arguments:
#         p {Pa} -- pressure
#         h {J/kg} -- [description]
#     """

#     # Documentation says "HEOS::Water" is the most accurate, but it is also slow.
#     # "IF97::Water" is faster but slightly less accurate
#     T = PropsSI('T','P',p,'H',h,"HEOS::Water")

#     return T

def get_specific_gas_constant():

    molar_mass = PropsSI("MOLAR_MASS","HEOS::Water") # Returned in kg/mol
    gas_constant = PropsSI("GAS_CONSTANT","HEOS::Water") # Returned in J/(mol*K)

    return gas_constant/molar_mass