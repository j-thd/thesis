# File to quickly calculate heat losses of some thrusters, based on their geometry
import math
import numpy as np
import matplotlib.pyplot as plt

import thrusters.thruster_data
import basic.chamber as chamber

td = thrusters.thruster_data.Cen2010_6 # [-]

T_wall = np.linspace(start=400, stop=1000, num=50) # [K]
channel_amount = td['channel_amount'] # [-] Amount of channels
w_channel = td['w_channel'] # [m] Channel width
w_channel_spacing = td['w_channel_spacing'] # [m] Spacing between channels
w_throat = td['w_throat'] # [m] Throat width
convergent_half_angle = td['convergent_half_angle']  # [rad] Half angle of convergent nozzle
divergent_half_angle = td['divergent_half_angle'] # [rad] Half angle of divergent nozzle
AR_exit = td['AR_exit'] # [-] Area ratio of nozzle exit
l_channel = td['L_channel']
l_inlet = 4e-3 #  [m] Inlet manifold length
w_outer_margin = 0.5e-3 # [m] Outer margin around chip 
emmisivity = 1 # [-] Emmisivity of the chip


l_outlet = chamber.outlet_length(
    w_channel=w_channel,
    w_channel_spacing=w_channel_spacing,
    channel_amount=channel_amount,
    convergent_half_angle=convergent_half_angle,
    w_throat=w_throat,
    divergent_half_angle=divergent_half_angle,
    AR_exit=AR_exit
)

l_chip = chamber.total_chip_length(l_inlet=l_inlet, l_channel=l_channel, l_outlet=l_outlet) # [m] Total length of chip
w_chip = chamber.total_chip_width(channel_amount=channel_amount, w_channel=w_channel, w_channel_spacing=w_channel_spacing, w_outer_margin=w_outer_margin) # [m] Total width of chip
A_chip = l_chip*w_chip # [m^2] Area of the chip

# Must use special function to calculate radation loss in numpy arrays as T**4 easily overflows
P_rad_top_wall = chamber.radiation_loss_numpy(T_np_array=T_wall, A=A_chip, emmisivity=emmisivity) # [W] Top-side heat loss
P_substrate_loss = chamber.basic_substrate_heat_loss(T_top=T_wall[1],kappa=1.2, emissivity=0.1, thickness=0.1e-3, A_substrate=A_chip)
print(P_rad_top_wall[1]) 
print(P_substrate_loss)

# Report losses and geometry
print("\n --- CHIP GEOMETRY ---")
print("Chip length: {:1.3f} mm" .format(l_chip*1e3))
print("Chip width: {:1.3f} mm" .format(w_chip*1e3))
print("Chip area: {:2.3f} mm^2".format(A_chip*1e6))

# # Plot heat loss as function of wall temperature
# plt.figure()
# plt.plot(T_wall, P_rad_top_wall)
# plt.xlabel("Wall temperature $T_{{wall}}$ [K]")
# plt.ylabel("Radiation loss $P_{{rad}}$ [W] ")
# plt.title("One-sided radation loss (top wall)")
# plt.grid()
# plt.show()