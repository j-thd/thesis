# File to get quick and dirty results from the regular good old IRT

import basic.IRT as IRT
from thermo.prop import FluidProperties
import thrusters.thruster_data

td = thrusters.thruster_data.Cen2010_1

# Design paramaters
fp = FluidProperties(td['propellant']) # Object to access fluid properties with
p_c = td['p_inlet'] # [bar] Chamber pressure
T_c = td['T_chamber_guess'] # [K] Chamber temperature
h_channel = td['h_channel'] # [m] Channel/nozzle depth
w_throat = td['w_throat'] # [m] Throat width
AR_exit = td['AR_exit'] # [-] Exit area ratio
p_back = td['p_back'] # [Pa] Atmospheric pressire

# Calculate throat area, and propellant properties
A_throat = h_channel*w_throat # [m^2] Thrpat area
gamma = fp.get_specific_heat_ratio(T=T_c, p=p_c) # [-] Get gamma at specified gas constant
R = fp.get_specific_gas_constant() # [J/(kg*K)] Specific gas constant

ep = IRT.get_engine_performance(p_chamber=p_c, T_chamber=T_c, A_throat=A_throat, AR_exit=AR_exit, p_back=p_back, gamma=gamma, R=R)
print(gamma)
print(ep)
print(ep['thrust']/ep['m_dot']/9.807)
print(ep['u_exit']/9.807)