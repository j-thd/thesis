# File to get quick and dirty results from the regular good old IRT

import basic.IRT as IRT
from thermo.prop import FluidProperties

# Design paramaters
fp = FluidProperties("water") # Object to access fluid properties with
p_c = 5e5 # [bar] Chamber pressure
T_c = 473 # [K] Chamber temperature
h_channel = 100e-6 # [m] Channel/nozzle depth
w_throat = 45e-6 # [m] Throat width
AR_exit = 16.971 # [-] Exit area ratio
p_back = 30 # [Pa] Atmospheric pressire

# Calculate throat area, and propellant properties
A_throat = h_channel*w_throat # [m^2] Thrpat area
gamma = 1.33 #fp.get_specific_heat_ratio(T=T_c, p=p_c) # [-] Get gamma at specified gas constant
R = fp.get_specific_gas_constant() # [J/(kg*K)] Specific gas constant

ep = IRT.get_engine_performance(p_chamber=p_c, T_chamber=T_c, A_throat=A_throat, AR_exit=AR_exit, p_back=30, gamma=gamma, R=R)
print(gamma)
print(ep)
print(ep['thrust']/ep['m_dot']/9.807)
print(ep['u_exit']/9.807)