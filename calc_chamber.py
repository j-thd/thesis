""" Kladblok voor berekeningen aan VLM prestaties
"""

import basic.chamber as chamber
import basic.IRT as IRT
import thermo.water as water
from constants import atm

# Do calculation on Silva2017 engine - Thruster 5

p_chamber = 2.05e5 # Not sure how it was measured
p_inlet = p_chamber # Assuming it is same through the chamber
p_back = 1.05e5 # Was it vacuum?
T_chamber = 273.15+24 # Measured chamber temperature (through resistance measurement)
T_inlet = 273.15 + 24# For now assuming room temperature
gamma = 1.4

R = 297#water.get_specific_gas_constant()
throat_width = 130e-6#45e-6 # Measured value
exit_width = 130*8.25e-6
channel_depth = 500e-6

# Calculate areas
A_throat = throat_width*channel_depth
A_exit = exit_width*channel_depth
AR_exit = A_exit/A_throat

# Calculate engine peformance with IRT
try:
    ep = IRT.get_engine_performance(p_chamber=p_chamber, T_chamber=T_chamber, A_throat=A_throat, AR_exit=AR_exit, p_back=p_back,gamma=gamma,R=R)
except ValueError as ve:
    print("Errors raised.")

# Calculate ideal power consumption
P_mh = chamber.ideal_power_consumption(mass_flow=ep['m_dot'],T_inlet=T_inlet,p_inlet=p_inlet,T_outlet=T_chamber,p_outlet=p_chamber)
print("Thrust: {:.2f} mN".format(ep['thrust']*1e3))
print("Isp: {:.1f} s".format(ep['u_exit']/9.80655))
print("Mass flow: {:.2f} mg/s".format(ep['m_dot']*1e6))
print("Power {:1f} W".format(P_mh))

# Add chip area
A_chip = 1291287398127491824
k =1293710947104  