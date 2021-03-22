'''Validate different Nusselt number relations with da Silva's thruster
There is no clear wall temperature, so, instead it must be checked against a range of temperature and the worst solution is assumed to be the true error.

Chamber temperature is raised by 1 degree Kelvin to fix that the outlet is otherwise not saturated despite observations
'''
import numpy as np
import matplotlib.pyplot as plt
import models.zero_D as zD
import thrusters.thruster_data
import thermo.convection
from thermo.prop import FluidProperties

td = thrusters.thruster_data.td_Silva_5 # Dictionary with design/measured values

# First set of Nusselt number relations (Dittus-Boelter for both phases)
Nu_func_gas_1 = thermo.convection.Nu_DB # [-] Function to calculate Nusselt number
Nu_func_liquid_1 = thermo.convection.Nu_DB
# Fully developed laminar flow for liquid phase. DB for the rest
Nu_func_gas_2 = thermo.convection.Nu_DB
Nu_func_liquid_2 = thermo.convection.Nu_laminar_developed_constant_wall_temp_square
# Third set (Kandlikar for x=0 for liquid phase)
Nu_func_gas_3 = thermo.convection.Nu_DB
Nu_func_liquid_3 = thermo.convection.Nu_Kandlikar_NDB_Re_low_sat_gas_constant_wall_temp_square_water


#bla
# For da Silva's thrusters the numbers have to fudged a bit, because the reported temperatures
# seem inconsitent with saturation temperatures and/or reported wall temperatures

w_channel = td['w_channel']             # [m] Channel width
T_inlet = td['T_inlet']                 # [K] Inlet temperature
T_chamber = td['T_chamber']+1            # [K] Chamber temperature
p_inlet = td['p_inlet']                 # [Pa] Inlet pressure
m_dot = td['m_dot']                     # [kg/s] Mass flow (through all channels if multiple)
channel_amount = td['channel_amount']   # [-] Amount of channels
h_channel = td['h_channel']             # [m] Channel height/depth
fp = FluidProperties(td['propellant'])  # Object from which fluid properties can be accessed



# Calculate mass flow for one single channel
m_dot_channel = m_dot/channel_amount    # [kg/s] Mass flow through one single channel

## Wall temperature is unknown so a range is chosen instead:
T_wall = np.linspace(start=T_chamber+1, stop=700, num=250) # [K] 

it = np.nditer(T_wall, flags=['c_index'])

L_channel_1 = np.zeros_like(T_wall) # [m]
L_channel_2 = np.zeros_like(T_wall) # [m]
L_channel_3 = np.zeros_like(T_wall) # [m]

for T in it:
    ## First set of Nusselt relations
    res_1 = zD.two_phase_single_channel( T_wall=T, w_channel=w_channel, Nu_func_gas=Nu_func_gas_1, Nu_func_liquid=Nu_func_liquid_1,\
    T_inlet=T_inlet, T_chamber=T_chamber, p_ref=p_inlet, m_dot=m_dot_channel,\
        h_channel=h_channel, fp=fp,print_info=False)
    # Store results
    L_channel_1[it.index] = res_1['L_channel']
    
    ## Second set of Nusselt relations
    res_2 = zD.two_phase_single_channel( T_wall=T, w_channel=w_channel, Nu_func_gas=Nu_func_gas_2, Nu_func_liquid=Nu_func_liquid_2,\
    T_inlet=T_inlet, T_chamber=T_chamber, p_ref=p_inlet, m_dot=m_dot_channel,\
        h_channel=h_channel, fp=fp,print_info=False)
    # Store results
    L_channel_2[it.index] = res_2['L_channel']

    ## Third set of Nusselt relations
    res_3 = zD.two_phase_single_channel( T_wall=T, w_channel=w_channel, Nu_func_gas=Nu_func_gas_3, Nu_func_liquid=Nu_func_liquid_3,\
    T_inlet=T_inlet, T_chamber=T_chamber, p_ref=p_inlet, m_dot=m_dot_channel,\
        h_channel=h_channel, fp=fp,print_info=False)
    # Store results
    L_channel_3[it.index] = res_3['L_channel'] 

#print(T_wall)

plt.plot(T_wall, L_channel_1*1e3, label="DB - DB")
plt.plot(T_wall, L_channel_2*1e3, label="FD Laminar - DB")
plt.plot(T_wall, L_channel_3*1e3, label="Kandlikar Re<100 - DB")
plt.xlabel("Wall temperature [K]")
plt.ylabel("Length channel [mm] ")
plt.hlines(td['L_channel']*1e3, xmin=T_wall[0], xmax=T_wall[-1], linestyle='dashed', label="Real length")
plt.title("0D Two-phase model applied to Silva #5")
plt.tight_layout(pad=1.0)
plt.legend()
plt.grid()
plt.show()