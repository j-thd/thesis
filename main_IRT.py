# File to get quick and dirty results from the regular good old IRT

import basic.IRT as IRT
from physical_constants import g0
from thermo.prop import FluidProperties
import thrusters.thruster_data

td = thrusters.thruster_data.Cen2010_6

# Design paramaters
fp = FluidProperties(td['propellant']) # Object to access fluid properties with
p_c = td['p_inlet'] # [bar] Chamber pressure
T_c = td['T_chamber_guess'] # [K] Chamber temperature
h_channel = td['h_channel'] # [m] Channel/nozzle depth
w_throat = td['w_throat'] # [m] Throat width
AR_exit = td['AR_exit'] # [-] Exit area ratio
p_back = td['p_back'] # [Pa] Atmospheric pressire

print("Chamber temperature: {:3.2f} K".format(T_c))

# Calculate throat area, and propellant properties
A_throat = h_channel*w_throat # [m^2] Thrpat area
gamma = fp.get_specific_heat_ratio(T=T_c, p=p_c) # [-] Get gamma at specified gas constant
R = fp.get_specific_gas_constant() # [J/(kg*K)] Specific gas constant

ep = IRT.get_engine_performance(p_chamber=p_c, T_chamber=T_c, A_throat=A_throat, AR_exit=AR_exit, p_back=p_back, gamma=gamma, R=R)
print("\n --- IRT predictions --- ")
print("Isp: {:3.2f} s".format(ep['Isp']))
print("Mass flow: {:2.3f} mg/s".format(ep['m_dot']*1e6))
print("Thrust: {:2.3f} mN".format(ep['thrust']*1e3))

zeta_CF = td['F']/ep['thrust']
Isp_real = td['F']/td['m_dot']/g0
zeta_Isp = Isp_real/ep['Isp']
discharge_factor = td['m_dot']/ep['m_dot']

print('\n --- Experimental values ---')
print("Isp: {:3.2f} s (zeta_Isp = {:1.3f})".format(Isp_real, zeta_Isp))
print("Mass flow: {:2.3f} mg/s (Cd = {:1.3f})".format(td['m_dot']*1e6, discharge_factor))
print("Thrust: {:2.3f} mN (zeta_CF ={:1.3f})".format(td['F'], zeta_CF))
