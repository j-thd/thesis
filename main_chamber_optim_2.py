''' Same as the initial iteration, but now fixing chamber temperature, and varying wall temperature.
    To get an idea about how power losses change with varying wall temperatures
'''
from thermo.convection import Nu_DB, Stanton_from_Nusselt_func_and_velocity
import numpy as np
from thermo.prop import FluidProperties
from basic.chamber import h_conv_from_Stanton, hydraulic_diameter_rectangular, ideal_enthalpy_change, radiation_loss, required_heater_area, required_power, velocity_from_mass_flow
import matplotlib.pyplot as plt


# Given parameters
fp = FluidProperties("water") # Object to get water properties from
T_inlet = 300 # [K] Temperature at inlet (room temperature)
p_inlet = 5e5 # [Pa] Pressure at inlet
T_chamber = 600 # [K] Desired exit temperature
T_wall = T_chamber + np.array([0, 100, 200, 300, 400]) # [K] Given wall temperature
h_channel = 100e-6 # [m] Channel depth
m_dot = 1e-6 # [kg/s] Mass flow
# Radiation loss parameters
emmisivity = 1 # [-] Emmisivity in Stefan Boltzmann law
# Used Nusselt relation
Nu_func = Nu_DB


# Which channel width gives the least power usage?
w_channel_min = 0.5e-5 # [m] Minimum channel width
w_channel_max = 5e-5 # [m] Maximum channel width
w_channel = np.linspace(start=w_channel_min, stop=w_channel_max) # [m] Range of channel widths to evaluate

# Calculate some basic channel geometry
A_channel = w_channel*h_channel # [m^2] Cross-sectional area of channel
Dh_channel = hydraulic_diameter_rectangular(w_channel=w_channel, h_channel=h_channel) # [m] Hydraulic diameter of rectangular channels

# Calculate some basic flow parameter at the inlet
rho_inlet = fp.get_density(T=T_inlet, p=p_inlet) # [kg/m^3]
u_inlet = velocity_from_mass_flow(m_dot=m_dot, rho=rho_inlet, A=A_channel) # [m/s] Flow velocity inside channel



# Storing results
Re_bulk =np.zeros((T_wall.shape[0], w_channel.shape[0])) # [-] Bulk Reynolds number
Pr_bulk = np.zeros_like(Re_bulk)  # [-] Bulk Prandtl number
u_bulk = np.zeros_like(Re_bulk) # [-] Velocity under bulk condition
Nusselt = np.zeros_like(Re_bulk) # [-] Nusselt number
Stanton = np.zeros_like(Re_bulk) # [-] Stanton number
h_conv = np.zeros_like(Re_bulk) # [W/(m^2*K)] Convective heat transfer coefficient
A_heater = np.zeros_like(Re_bulk) # [m^2] Heater area required to raise temperature to T_chamber with temperature T_wall
L_channel = np.zeros_like(Re_bulk) # [m] Required channel length for A_heater and w_channel
P_radiation_loss = np.zeros_like(Re_bulk) # [W] Power lost through radation of A_heater with temperature T_wall
P_total = np.zeros_like(Re_bulk) # [W] Total power consumption

### Calculate the tihngs that vary both with T_wall and w_channel
it_T = np.nditer(T_wall,flags=['c_index']) # [K] To iterate over all T_chamber values

for T in it_T:
    ## Make sure all parameters are calculated at the same reference state (including velocity)
    # Calculate all the bulk conditions required to determine the Nusselt number relation
    T_bulk = (T_inlet + T_chamber) / 2 # [K] Bulk temperature as reference
    rho_bulk = fp.get_density(T=T_bulk, p=p_inlet) # [kg/m^3] Bulk density is required to get the correct velocity at the reference temperature
    u_bulk[it_T.index,:] = velocity_from_mass_flow(m_dot=m_dot, rho=rho_bulk, A=A_channel) # [m/s] The velocity at the bulk reference state
    Re_bulk[it_T.index,:] = fp.get_Reynolds_from_mass_flow(T=T_bulk,p=p_inlet, L_ref=Dh_channel, m_dot=m_dot, A=A_channel) # [-]
    Pr_bulk[it_T.index,:] = fp.get_Prandtl(T=T_bulk, p=p_inlet) # [-]
    # Now the convection parameters can be determined, taking care that the reference temperature and flow conditions are evaluated in the same location
    Nusselt[it_T.index,:] = Nu_func(args={'Re': Re_bulk[it_T.index,:], 'Pr': Pr_bulk[it_T.index,:]}) # [-] Not necessary, but useful for plots
    Stanton[it_T.index,:] = Stanton_from_Nusselt_func_and_velocity(Nu_func=Nu_func, m_dot=m_dot, A=A_channel, T_ref=T_bulk, p_ref=p_inlet, L_ref=Dh_channel, fp=fp) # [-]
    h_conv[it_T.index,:] = h_conv_from_Stanton(Stanton=Stanton[it_T.index,:], u=u_bulk[it_T.index,:], T_ref=T_bulk, p_ref=p_inlet, fp=fp) # [W/(m^2*K)]
    # Know we must know how much energy is convected, we can determine the required heater area and channel length
    delta_h = ideal_enthalpy_change(T_inlet=T_inlet, p_inlet=p_inlet, T_outlet=T_chamber, p_outlet=p_inlet, fp=fp) # [J/kg] Specific enthalpy change of the channel
    Q_dot = required_power(m_dot=m_dot, delta_h=delta_h) # [W] Power that must go into the flow to achieve delta_h 
    A_heater[it_T.index,:] = required_heater_area(Q_dot=Q_dot, h_conv=h_conv[it_T.index,:], T_wall=T, T_ref=T_bulk) # [m^2] Required area to deliver power Q_dot
    L_channel[it_T.index,:] = A_heater[it_T.index,:]/w_channel # [m] With this length , the required heater area is obtained for the given channel width
    # Estimated power loss through radation
    # Under the simple assumption that the heater only loses power through radiation, the power loss is calculated
    P_radiation_loss[it_T.index,:] = radiation_loss(T=float(T), A=A_heater[it_T.index,:], emmisivity=emmisivity) # [W] Radiation from single-sided heater surface
    P_total[it_T.index,:] = P_radiation_loss[it_T.index,:] + Q_dot

## PLOTTING RESULTS
fig0, axs0 = plt.subplots(3,1) # For all the flow similarity parameters
fig1, axs1 = plt.subplots(3,1) # For the resulting dimensions and power loss
# First plot the results that do not vary with T_chamber, and need only be plotted once


# Reset the iterator to plot all results with different lines for different chamber temperatures
it_T.reset()
for T in it_T:
    ## First set of subplots
    # Bulk velocity 
    axs0[0].plot(w_channel*1e6, u_bulk[it_T.index,:], label="{:3.0f} K".format(T))
    # Bulk Reynolds
    axs0[1].plot(w_channel*1e6, Re_bulk[it_T.index,:])
    # Bulk Prandtl
    axs0[2].plot(w_channel*1e6, Pr_bulk[it_T.index,:])
    ## Second set of subplots
    axs1[0].plot(w_channel*1e6, L_channel[it_T.index,:]*1e3, label="{:3.0f} K".format(T))
    axs1[1].plot(w_channel*1e6, P_radiation_loss[it_T.index,:], label="{:3.0f} K".format(T))
    axs1[2].plot(w_channel*1e6, P_total[it_T.index,:])
    
    

# Activate grid, and set labels
## First set of subplots
[a.grid() for a in axs0 ]  
axs0[0].set_ylabel("$u_{bulk}$ [m/s]")
axs0[1].set_ylabel("$Re_{bulk}$ [-]")
axs0[2].set_ylabel("$Pr_{bulk}$ [-]")
# Set x-label only under bottom axis
axs0[2].set_xlabel("Channel width $w_c$ [$\\mu$m]")
# Make layout tight, add title and legend
fig0.suptitle("Flow and heat transfer parameters for $T_{{chamber}} = {:4.0f}$ K & $\dot{{m}}= {:3.1f}$ mg/s".format(T_chamber, m_dot*1e6))
fig0.tight_layout()
axs0[0].legend(title="Wall temp. $T_w$:", fontsize="small")
## Second set of subplots
[a.grid() for a in axs1]
axs1[0].set_ylabel("$L_c$ [mm]")
axs1[1].set_ylabel("$P_{rad}$ [W]")
axs1[2].set_ylabel("$P_{t}$ [W]")
# Set x-label only under bottom axis
axs1[2].set_xlabel("Channel width $w_c$ [$\\mu$m]")
# Make layout tight, add title and legend
fig1.suptitle("Channel length and power consumption for $T_{{chamber}} = {:4.0f}$ K & $\dot{{m}}= {:3.1f}$ mg/s".format(T_chamber, m_dot*1e6))
fig1.tight_layout()
axs1[1].legend(title="Wall temp. $T_c$:", fontsize="small")



# Show results
plt.show()