# File to validate calculations for DaSilva's thruster #5

import numpy as np
import matplotlib.pyplot as plt

import models.one_D as oneD
import thrusters.thruster_data
import thermo.convection
import thermo.two_phase as tp
import basic.chamber
from thermo.prop import FluidProperties

td = thrusters.thruster_data.td_Silva_5 # Dictionary with design/measured values

# Fidelity of simulatoin
steps_per_section = 200 # [-] Amount of subdivision in each section 
steps_l = steps_per_section
steps_tp = steps_per_section
steps_g = steps_per_section

# Functions to calculate Nusselt numbers.
Nusselt_relations_1 = {
    'Nu_func_gas': thermo.convection.Nu_DB, # [-] Function to calculate Nusselt number (gas phase)
    'Nu_func_liquid': thermo.convection.Nu_DB,  # [-] Function to caculate Nusselt number (liquid phase)
    'Nu_func_two_phase': tp.Nu_Kandlikar_NBD_dryout, # [-] Function to calculate Nusselt number (two-phase)
    'Nu_func_le': thermo.convection.Nu_DB, # [-] Function to calculate Nusselt in two-phase, AS IF the flow was entirely liquid (two-phase, le)
    'Nu_func_dryout': thermo.two_phase.Nu_DB_two_phase, #thermo.two_phase.Nu_DB_two_phase # [-] Function to calculate the Nusselt number after dry-out. It is up to Nu_func_two_phase to decide if/how to apply it'
}

Nusselt_relations_2 = {
    'Nu_func_gas': thermo.convection.Nu_laminar_developed_constant_wall_temp_square, # [-] Function to calculate Nusselt number (gas phase)
    'Nu_func_liquid': thermo.convection.Nu_laminar_developed_constant_wall_temp_square,  # [-] Function to caculate Nusselt number (liquid phase)
    'Nu_func_two_phase': tp.Nu_Kandlikar_NBD_dryout, # [-] Function to calculate Nusselt number (two-phase)
    'Nu_func_le': thermo.convection.Nu_laminar_developed_constant_wall_temp_square, # [-] Function to calculate Nusselt in two-phase, AS IF the flow was entirely liquid (two-phase, le)
    'Nu_func_dryout': thermo.convection.Nu_laminar_developed_constant_wall_temp_square, #thermo.two_phase.Nu_DB_two_phase # [-] Function to calculate the Nusselt number after dry-out. It is up to Nu_func_two_phase to decide if/how to apply it'
}




# For da Silva's thrusters the numbers have to fudged a bit, because the reported temperatures
# seem inconsitent with saturation temperatures and/or reported wall temperatures
w_channel = td['w_channel']             # [m] Channel width
T_inlet = td['T_inlet']                 # [K] Inlet temperature
T_chamber = td['T_chamber'] +0.5          # [K] Chamber temperature <--- NOTE: this is where the number is increased by 0.5K, because it's below saturatoin temperature otherwise
p_inlet = td['p_inlet']                 # [Pa] Inlet pressure
m_dot = td['m_dot']                     # [kg/s] Mass flow (through all channels if multiple)
channel_amount = td['channel_amount']   # [-] Amount of channels
h_channel = td['h_channel']             # [m] Channel height/depth
fp = FluidProperties(td['propellant'])  # Object from which fluid properties can be accessed

# Calculate mass flow for one single channel
m_dot_channel = m_dot/channel_amount    # [kg/s] Mass flow through one single channel

# Wall temperature is unknown, so a range is taken
T_wall = np.linspace(start=T_chamber+0.01, stop=600, num=5000) # [K]

# Geometric values
wetted_perimeter = basic.chamber.wetted_perimeter_rectangular(w_channel=w_channel, h_channel=h_channel) # [m] Wetted perimeter of channel
A_channel = w_channel*h_channel # [m^2] Cross-sectional through which fluid flows
D_hydraulic = basic.chamber.hydraulic_diameter_rectangular(w_channel=w_channel, h_channel=h_channel) # [m] Hydraulic diameter
# Preparation functions calculations many intermediate values that are known before geometry is known
prepared_values = oneD.full_homogenous_preparation(
    T_inlet=T_inlet,
    T_outlet=T_chamber,
    m_dot=m_dot,
    p_ref=p_inlet,
    steps_l=steps_l,
    steps_tp=steps_tp,
    steps_g=steps_g,
    fp=fp)

 
# Storing the length results in here, one for each set of Nusselt relations
L_1 = np.zeros_like(T_wall) # [m] Total channel length
L_2 = np.zeros_like(L_1)

# Loop to calculate the channel length with each wall temperature
it_T = np.nditer(T_wall, flags=['c_index']) # [K] Wall temperature

for T in it_T:
    results_1 = oneD.full_homogenous_calculation(
        prepared_values=prepared_values,
        Nusselt_relations=Nusselt_relations_1,
        A_channel=A_channel,
        wetted_perimeter=wetted_perimeter,
        D_hydraulic=D_hydraulic,
        m_dot=m_dot,
        T_wall=T, # <--- Iterated variable
        p_ref=p_inlet,
        fp=fp
        )
    
    results_2 = oneD.full_homogenous_calculation(
        prepared_values=prepared_values,
        Nusselt_relations=Nusselt_relations_2,
        A_channel=A_channel,
        wetted_perimeter=wetted_perimeter,
        D_hydraulic=D_hydraulic,
        m_dot=m_dot,
        T_wall=T, # <--- Iterated variable
        p_ref=p_inlet,
        fp=fp
        )

    L_1[it_T.index] = results_1['L_total']
    L_2[it_T.index] = results_2['L_total']

print("Thruster name:")
print(td['name'])

plt.figure()
plt.title("1D Two-phase model applied to Silva #5")
plt.plot(T_wall,L_1*1e3, label="Turbulent")
plt.plot(T_wall,L_2*1e3, label="Laminar")
plt.xlabel("Wall temperature - $T_{{wall}}$ [K]")
plt.ylabel("Channel length - $L$ [m]")
plt.hlines(td['L_channel']*1e3, xmin=T_wall[0], xmax=T_wall[-1], linestyle='dashed', color='red', label="Real length")
plt.ylim([0,td['L_channel']*2e3])
plt.grid()
plt.legend()
plt.tight_layout(pad=1.0)
plt.show()