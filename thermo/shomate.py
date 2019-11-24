# Implements Shomate equation and stores parameters for it


def getEnthalpy(parameters, T):
    # Inputs
    #   T [K] - Temperature of substance
    #  
    # Output:
    #   H [kJ/kg] - Specific enthalpy
    # Enthalpy given w.r.t to reference state, but unclear how that state is
    # defined so far
    t = T/1000
    H = (parameters['A']*t
         + parameters['B']*t**2/2
         + parameters['C']*t**3/3
         + parameters['D']*t**4/4
         - parameters['E']/t
         + parameters['F']
         - parameters['H'])
    return H


# Naming is given by <common-name>_<phase(L/G)>_(1/2/3) (for different ranges)
water_L = {'A': -203.6060,  # Parameters as given by NIST Webbook
           'B': 1523.290,
           'C': -3196.413,
           'D': 2474.455,
           'E': 3.855326,
           'F': -256.5478,
           'G': -488.7163,
           'H': -285.8304,
           'validity_range': (298, 500)}  # Range in which parameters are valid

water_G_1 = {'A': 30.09200,
             'B': 6.832514,
             'C': 6.793435,
             'D': -2.534480,
             'E': 0.082139,
             'F': -250.8810,
             'G': 223.3967,
             'H': 241.8264,
             'validity_range': (500, 1700)}
