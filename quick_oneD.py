# File to do some quick calculations and tests with the one-dimensional

import numpy as np
import matplotlib.pyplot as plt

import models.one_D as oneD
import thrusters.thruster_data
import thermo.convection
import thermo.two_phase as tp
import basic.chamber
from thermo.prop import FluidProperties

td = thrusters.thruster_data.td_Silva_5 # Dictionary with design/measured values

# Functions to calculate Nusselt numbers
Nu_func_gas = thermo.convection.Nu_DB # [-] Function to calculate Nusselt number (gas phase)
Nu_func_liquid = thermo.convection.Nu_DB  # [-] Function to caculate Nusselt number (liquid phase)
Nu_func_two_phase = tp.Nu_Kandlikar_NBD_dryout # [-] Function to calculate Nusselt number (two-phase)
Nu_func_le = thermo.convection.Nu_DB # [-] Function to calculate Nusselt in two-phase, AS IF the flow was entirely liquid (two-phase, le)
Nu_func_dryout = thermo.two_phase.Nu_DB_two_phase #thermo.two_phase.Nu_DB_two_phase # [-] Function to calculate the Nusselt number after dry-out. It is up to Nu_func_two_phase to decide if/how to apply it


#bla
# For da Silva's thrusters the numbers have to fudged a bit, because the reported temperatures
# seem inconsitent with saturation temperatures and/or reported wall temperatures
T_wall = td['T_wall'] + 100                 # [K] Wall temperature
w_channel = td['w_channel']             # [m] Channel width
T_inlet = td['T_inlet']                 # [K] Inlet temperature
T_chamber = td['T_chamber']+100         # [K] Chamber temperature
p_inlet = td['p_inlet']                 # [Pa] Inlet pressure
m_dot = td['m_dot']                     # [kg/s] Mass flow (through all channels if multiple)
channel_amount = td['channel_amount']   # [-] Amount of channels
h_channel = td['h_channel']             # [m] Channel height/depth
fp = FluidProperties(td['propellant'])  # Object from which fluid properties can be accessed

# Calculate mass flow for one single channel
m_dot_channel = m_dot/channel_amount    # [kg/s] Mass flow through one single channel

# Geometric values
wetted_perimeter = basic.chamber.wetter_perimeter_rectangular(w_channel=w_channel, h_channel=h_channel) # [m] Wetted perimeter of channel
A_channel = w_channel*h_channel # [m^2] Cross-sectional through which fluid flows
D_hydr = basic.chamber.hydraulic_diameter_rectangular(w_channel=w_channel, h_channel=h_channel) # [m] Hydraulic diameter
# Preparation functions calculations many intermediate values that are known before geometry is known
p_l = oneD.prepare_single_phase_liquid(T_inlet=T_inlet, steps=50, p_ref=p_inlet, m_dot=m_dot, fp=fp) # {dict} Dictionary of prepared values
res_l = oneD.calc_channel_single_phase(\
    T = p_l['T'],
    Q_dot= p_l['Q_dot'],
    rho = p_l['rho'],
    Pr = p_l['Pr'],
    kappa = p_l['kappa'],
    mu = p_l['mu'],
    p_ref=p_inlet,
    m_dot=m_dot,
    T_wall=T_wall,
    D_hydr=D_hydr,
    wetted_perimeter=wetted_perimeter,
    A_channel=A_channel,
    Nu_func=Nu_func_liquid,
    fp=fp
    )

# Store results for plotting or adding new points
T = p_l['T'] # [K] Flow temperature 
L = res_l['L'] # [m] Channel length x at given temperature T
Nu = res_l['Nu'] # [-] Nusselt number at given temperature T
Re = res_l['Re'] # [-] Reynolds number at given temperature T
Pr = p_l['Pr'] # [-] Prandtl number at given temperature T
mu = p_l['mu'] # [Pa*s] Viscosity
rho = p_l['rho'] # [kg/m^3] Flow density
u = res_l['u'] # [m/s] Flow velocity
h_conv = res_l['h_conv'] # [W/(m*K)] Convective heat transfer parameter
x = np.zeros_like(T) # [-] Flow quality, zero before boiling
alpha = np.zeros_like(T) # [-] Void fraction, zero before boiling

## Intermediate liquid phase result
print("Liquid phase length: {:2.3f} mm".format(1e3*L[-1]))

## Now do the same for the transition from x=0 to x=1. The subscript tp is necessary for when the results are joined together, 
# and X_tp_l means that is at scalar value X at the start of the two-phase transition when it is still complete liquid phase (where X_l is the resulting X in the liquid phase flow)
# Example: rho_l is how the density changes as the temperature increases towards the saturation point. rho_tp_l is the density at this saturation point, rho_tp is how the density changes through-out the phase change
p_tp = oneD.prepare_homogenous_transition(p=p_inlet, m_dot=m_dot, steps=500, fp=fp)
x_tp = p_tp['x']
alpha_tp = p_tp['alpha']
T_sat = p_tp['T_sat']
rho_tp_l = p_tp['rho_l']
rho_tp_g = p_tp['rho_g']
rho_tp = p_tp['rho']
mu_tp_l = p_tp['mu_l']
#mu_tp_g = p_tp['mu_g']
mu_tp = p_tp['mu']
Pr_tp_l = p_tp['Pr_l']
Pr_tp = p_tp['Pr']
kappa_tp_l = p_tp['kappa_l']
kappa_tp = p_tp['kappa']
Q_dot_tp = p_tp['Q_dot']


res_tp = oneD.calc_homogenous_transition(
    p_sat=p_inlet,
    x=x_tp,
    alpha=alpha_tp,
    T_sat=T_sat,
    rho_l=rho_tp_l,
    rho_g=rho_tp_g,
    rho=rho_tp,
    m_dot=m_dot,
    mu_l=mu_tp_l,
    mu=mu_tp,
    Pr_l=Pr_tp_l,
    Pr=Pr_tp,
    kappa_l=kappa_tp_l,
    kappa=kappa_tp,
    Q_dot=Q_dot_tp,
    T_wall=T_wall,
    D_hydr=D_hydr,
    wetted_perimeter=wetted_perimeter,
    A_channel=A_channel,
    Nu_func_tp=Nu_func_two_phase,
    Nu_func_le=Nu_func_le,
    Nu_func_dryout=Nu_func_dryout,
    fp=fp
)
## Intermediate two-phase results
print("Two-phase length: {:2.3f} mm".format(1e3*res_tp['L'][-1]))

T_tp = np.ones_like(x_tp) * T_sat # [K]
L_tp = res_tp['L'] + L[-1] # [m] Add length of previous section
u_tp = res_tp['u']
rho_tp = res_tp['rho']
Nu_tp = res_tp['Nu']
mu_tp = p_tp['mu']
Re_tp = res_tp['Re']
h_conv_tp = res_tp['h_conv']



## Add the values to already existing results (liquid + two-phase)


T = np.hstack((T, T_tp)) # [K]
L = np.hstack((L, L_tp)) # [m]
x = np.hstack((x, x_tp)) # [-]
alpha = np.hstack((alpha, alpha_tp)) # [-]
rho = np.hstack((rho, rho_tp)) # [kg/m^3]
u = np.hstack((u, u_tp)) # [-]
Re = np.hstack((Re, Re_tp)) # [-]
Pr = np.hstack((Pr, Pr_tp)) # [-]
mu = np.hstack((mu, mu_tp)) # [Pa*s]
Nu = np.hstack((Nu, Nu_tp)) # [-]
h_conv = np.hstack((h_conv, h_conv_tp))

## Gas Phase

# Prepare intermediate (non-geometric-dependent) values in channel
p_g = oneD.prepare_single_phase_gas(T_outlet=T_chamber, steps=50, p_ref=p_inlet, m_dot=m_dot, fp=fp)
T_g = p_g['T'] # [K]
rho_g = p_g['rho'] # [kg/m^3]
Pr_g = p_g['Pr'] # [-]
kappa_g = p_g['kappa'] # [W/(m*K)]
mu_g = p_g['mu'] # [Pa*s]

# Calculate channel length and other interesting values in gas phase
res_g = oneD.calc_channel_single_phase(
    T=T_g,
    Q_dot=p_g['Q_dot'],
    rho=rho_g,
    Pr=Pr_g,
    kappa=kappa_g,
    mu=mu_g,
    p_ref=p_inlet,
    m_dot=m_dot,
    T_wall=T_wall,
    D_hydr=D_hydr,
    wetted_perimeter=wetted_perimeter,
    A_channel=A_channel,
    Nu_func=Nu_func_gas,
    fp=fp
)

## Store results and add them to previous sections of channel
print("Gas phase length: {:2.3f} mm".format(1e3*res_g['L'][-1]))
T_g = p_g['T']
x_g = np.ones_like(T_g)
alpha_g = np.ones_like(T_g) 
L_g = res_g['L'] + L[-1] # [m] Add length of previous section
u_g = res_g['u']
rho_g = p_g['rho']
Nu_g = res_g['Nu']
mu_g = p_g['mu']
Re_g = res_g['Re']
Pr_g = p_g['Pr']
h_conv_g = res_g['h_conv']

# Check Nusselt-number jump
print("Nu jump: {:3.3f}".format(Nu_g[0]/Nu_tp[-1]))
print("Kappa_l/kappa_g: {:3.3f}".format(kappa_tp_l/kappa_g[0]))
print("h_conv_tp/h_conv_g: {:3.9f}".format(h_conv_g[0]/h_conv_tp[-1]))

T = np.hstack((T, T_g))
x = np.hstack((x, x_g))
alpha = np.hstack((alpha, alpha_g))
L = np.hstack((L, L_g))
u = np.hstack((u, u_g))
rho = np.hstack((rho, rho_g))
Nu = np.hstack((Nu, Nu_g))
mu = np.hstack((mu, mu_g))
Re = np.hstack((Re, Re_g))
Pr = np.hstack((Pr, Pr_g))
h_conv = np.hstack((h_conv, h_conv_g))

print("Lengtish {:10.10f} mm".format(1e3*L[-2]))
print("Length: {:1.3f} mm".format(1e3*L[-1]))
print("Last Nu {:2.34}".format(Nu[-1]))

fig, axs = plt.subplots(3,3)
# Top left
axs[0][0].plot(L*1e3, T)
axs[0][0].set_ylabel("Temperature T [K]")
axs[0][0].grid()
# Middle left
axs[1][0].plot(L*1e3, Nu)
axs[1][0].set_ylabel("Nu [-]")
axs[1][0].grid() 
# Bottom left
axs[2][0].plot(L*1e3, rho)
axs[2][0].set_ylabel("Density $\\rho$ [kg$\\cdot$m$^{{-3}}$]")
axs[2][0].set_xlabel("Distance $x$ [mm]")
axs[2][0].grid()
# Top center
axs[0][1].plot(L*1e3, mu)
axs[0][1].set_ylabel("mu [Pa*s]")
axs[0][1].grid()
# Middle center
axs[1][1].plot(L*1e3, Pr)
axs[1][1].set_ylabel("Prandtl Pr [-]")
axs[1][1].grid()
# Bottom center
axs[2][1].plot(L*1e3, u)
axs[2][1].set_ylabel("Velocity $u$ [m/s] ")
axs[2][1].set_xlabel("Distance $x$ [mm]")
axs[2][1].grid()
## Top right
axs[0][2].plot(L*1e3, h_conv)
axs[0][2].set_ylabel("Convective heat flux parameter [W/m*K] ")
axs[0][2].grid()
## Middle right
axs[1][2].plot(L*1e3, alpha)
axs[1][2].set_ylabel("Void fraction [-]")
axs[1][2].grid()
## Bottom right
axs[2][2].plot(L*1e3, Re)
axs[2][2].set_ylabel("Reynolds $Re$ [-]")
axs[2][2].grid()
axs[2][2].set_xlabel("Distance $x$ [mm]")


plt.suptitle("1D-single-channel model ($\\dot{{m}}={:2.2f}$ mg/s; $p={:1.1f}$ bar)".format(m_dot_channel*1e6,p_inlet*1e-5))
plt.tight_layout(pad=0.5)
plt.show()

# ## PLOT SOME TRANSITION SHIT
# plt.figure()
# plt.plot(x_tp, alpha_tp)
# plt.xlabel("Quality $x$ [-]")
# plt.ylabel("Void fraction $\\alpha$ [-]")
# plt.grid()

# plt.figure()
# plt.plot(x_tp, u_tp)
# plt.xlabel("Quality $x$ [-]")
# plt.ylabel("Flow velocity $u$ [m/s]")
# plt.grid()
# plt.show()Nu_func_tp