'''
    File that ties together IRT and chamber optimization based on Nusselt numbers to figure out the best power consumption for a given thrust
'''

from basic.chamber import hydraulic_diameter_rectangular, radiation_loss
import numpy as np
import matplotlib.pyplot as plt

from thermo.prop import FluidProperties
import models.zero_D as zD
import thermo.convection

## Design parameters
fp = FluidProperties("water") # Object to get properties of propellant
F = 2.5e-3 # [N] Desired thrust
p_inlet = 5e5 # [Pa] Inlet pressure, no pressure drop assumed over chamber
T_inlet = 300 # [K] Inlet temperature
p_back = 0 # [Pa] Into vacuum
h_channel = 100e-6 # [m] Channel depths
w_channel = 3e-3 # [m] Channel width
AR_exit = 10 # [-] Exit area ratio
# Chamber temperatures ranges to consider
T_chamber_low = 600 # [K]
T_chamber_high = 1000 # [K]
# Superheat at wall
T_wall_superheat = 100 # [K] How much hotter is wall than fluid?

## Model parameters 
Nu_func = thermo.convection.Nu_DB # Function used to calculate heat transfer
emmissivity = 0.7 # [-]

# A range of mass flows and chamber temperatures for the given thrust must be made, for which the chamber must then be optimized
T_chamber = np.linspace(start=T_chamber_low, stop=T_chamber_high)
#T_wall = T_chamber + T_wall_superheat # [K] Wall temperature?

# Every element corresponded to a solution for the given design parameters and chamber temperature
m_dot = np.zeros_like(T_chamber) # [kg/s] Mass flow
Isp = np.zeros_like(m_dot) # [s] Specific impulse
A_heater = np.zeros_like(m_dot) # [m^2] Calculated heater area
P_radiation_loss = np.zeros_like(m_dot) # [W] Radiation loss
P_total = np.zeros_like(m_dot) # [W] Total power

# Get engine performance for given T and F
it_T = np.nditer(T_chamber, flags=['c_index']) # Iterator that keeps track of index
for T in it_T: # Iteration is used because IRT code does not work well with numpy arrays yet
    T_c = float(T) # Convert this shit back to a python Float for safety purposes. If it remains an np object it could overflow
    ## ENGINE PERFORMANCE
    # Calculate the performance according to ideal theory
    ep = zD.engine_performance_from_F_and_T(F_desired=F, p_chamber=p_inlet, T_chamber=T_c, AR_exit=AR_exit, p_back=p_back, fp=fp)
    m_dot[it_T.index] = ep['m_dot'] # [kg/s] Mass flow
    Isp[it_T.index] = ep['Isp'] # [s]

    ## CHAMBER PERFORMANCE
    # Calculate the reference temperature for the Nusselt function used
    T_bulk = (T_inlet + T_c) / 2 # [K] Bulk temperature is the reference used for Dittus Boelter relation
    T_wall = T_c + T_wall_superheat # [K] Wall temperature is a bit higher than chamber temperature
    A_channel = h_channel * w_channel # [m^2] Channel diameter
    Dh_channel = hydraulic_diameter_rectangular(w_channel=w_channel, h_channel=h_channel) # [m] Hydraulic diameter is the reference length
    cp = zD.chamber_performance_from_Nu(Nu_func=Nu_func, T_inlet=T_inlet, T_chamber=T_c, T_ref=T_bulk, T_wall=T_wall, p_ref=p_inlet, m_dot=m_dot[it_T.index], A_channel=A_channel, L_ref=Dh_channel, fp=fp)
    A_heater[it_T.index] = cp['A_heater'] # [m^2] Calculated heater area
    P_radiation_loss[it_T.index] = radiation_loss(T=T_wall, A=A_heater[it_T.index], emmisivity=emmissivity) # [W]
    P_total[it_T.index] = cp['Q_dot'] + P_radiation_loss[it_T.index] # [W]

    

# 
plt.plot(T_chamber, P_total)
plt.title("Thruster performance for $F = {:2.2f} mN$".format(F*1e3))
plt.grid()
plt.xlabel("$T_c$ [K]")
plt.ylabel("Power [W] ")

plt.show()