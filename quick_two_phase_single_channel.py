# File to compare predictions of two-phase single channel calculations with real thrusters

import models.zero_D as zD
import thrusters.thruster_data
import thermo.convection
from thermo.prop import FluidProperties

td = thrusters.thruster_data.td_Silva_5 # Dictionary with design/measured values

Nu_func_gas = thermo.convection.Nu_Kandlikar_NDB_Re_low_sat_gas_constant_wall_temp_square_water # [-] Function to calculate Nusselt number
Nu_func_liquid = thermo.convection.Nu_DB

#bla
# For da Silva's thrusters the numbers have to fudged a bit, because the reported temperatures
# seem inconsitent with saturation temperatures and/or reported wall temperatures
T_wall = td['T_wall'] +50               # [K] Wall temperature
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

zD.two_phase_single_channel( T_wall=T_wall, w_channel=w_channel, Nu_func_gas=Nu_func_gas, Nu_func_liquid=Nu_func_liquid,\
    T_inlet=T_inlet, T_chamber=T_chamber, p_ref=p_inlet, m_dot=m_dot_channel,\
        h_channel=h_channel, fp=fp)

