from basic.IRT_corrections import  hydraulic_diameter
from basic.chamber import wetter_perimeter_rectangular
import numpy as np
from thermo.prop import FluidProperties
import models.zero_D as zD
import matplotlib.pyplot as plt


input_data = "NOPE" 

## Quick and dirty to swap between different input cases
if input_data == "Rajeev":
    fp = FluidProperties("nitrogen")
    # Design paramaters
    h_channel = 0.000081 # [m] Not necessary for calculations, but used as assumption to calculate throat width of a rectangular nozzle
    T_inlet = 300 # [K] Inlet temperature, used as assumption to calculate ideal power needed for heater
    p_chamber = 1.92e5 # [Pa] chamber pressure
    p_back = 0 # [Pa] expanded into vacuum
    AR_exit = np.array([30]) # [-] Range of area ratios to go over
    F = 0.46208e-3 # [mN] Desired thrust

    # Range of temperatures to go over
    T_chamber_min = 295.86 # [K]
    T_chamber_max = 1000 # [K]
    dT = (T_chamber_max - T_chamber_min)/50 # [K]
    T_chamber = np.arange(start=T_chamber_min, stop=T_chamber_max+dT/2, step=dT) # [K] Temperature range to plot (dT/2 added to include max in range)

else:
    fp = FluidProperties("HEOS::Water")
    # Design paramaters
    h_channel = 100e-6 # [m] Not necessary for calculations, but used as assumption to calculate throat width of a rectangular nozzle
    T_inlet = 300 # [K] Inlet temperature, used as assumption to calculate ideal power needed for heater
    p_chamber = 5e5 # [Pa] chamber pressure
    p_back = 0 # [Pa] expanded into vacuum
    AR_exit = np.array([5,10,25,50]) # [-] Range of area ratios to go over
    F = 2.5e-3 # [mN] Desired thrust

    # Range of temperatures to go over
    T_chamber_min = 800 # [K]
    T_chamber_max = 2000 # [K]
    dT = (T_chamber_max - T_chamber_min)/50 # [K]
    T_chamber = np.arange(start=T_chamber_min, stop=T_chamber_max+dT/2, step=dT) # [K] Temperature range to plot (dT/2 added to include max in range)

# Arrays to store results in
m_dot = np.zeros((AR_exit.shape[0],T_chamber.shape[0])) # [kg/s] Mass flow
A_throat = np.zeros_like(m_dot) # [m^2] Throat area
Isp = np.zeros_like(m_dot) # [s] Specific impulse
h_chamber = np.zeros_like(m_dot) # [J/kg] Hold enthalpy at chamber

# Results to check Reynolds number and condensation in throat
T_throat = np.zeros_like(m_dot) # [K] Temperature in throat
p_throat = np.zeros_like(m_dot) # [Pa] Pressure in throat
u_throat = np.zeros_like(m_dot) # [m/s] Flow velocity in throat
Re_throat = np.zeros_like(m_dot) # [-] Reynolds number in throat
Pr_throat = np.zeros_like(m_dot) # [-] Prandtl number in the throat
mu_throat = np.zeros_like(m_dot) # [Pa*s] Dynamic viscosity in the throat

# Results to check condensation in nozzle exit
T_exit = np.zeros_like(m_dot) # [K]
p_exit = np.zeros_like(m_dot) # [Pa]

it_AR = np.nditer(AR_exit, flags=['c_index'])

# Iterate over all possible combinations of area ratio and chamber temperature to find the mass flow and area ratio that finds according to basic IRT
for AR in it_AR:
    it_T = np.nditer(T_chamber, flags=['c_index'])
    for T in it_T:
        # Store mass flow, throat area and other results
        ep = zD.engine_performance_from_F_and_T(F_desired=F, p_chamber=p_chamber,T_chamber=float(T), AR_exit=float(AR),p_back=p_back, fp=fp)
        m_dot[it_AR.index][it_T.index] = ep['m_dot'] # [kg/s] Mass flow
        A_throat[it_AR.index][it_T.index] = ep['A_throat'] # [m^2] Throat area
        Isp[it_AR.index][it_T.index] = ep['Isp'] # [s] Specific impulse
        h_chamber[it_AR.index][it_T.index] = fp.get_enthalpy(T=float(T),p=p_chamber) # [J/kg] Enthalpy at chamber inlet
        # Throat stuff
        T_throat[it_AR.index][it_T.index] = ep['T_throat'] # [K] Throat temperature
        p_throat[it_AR.index][it_T.index] = ep['p_throat'] # [Pa] Throat pressure
        u_throat[it_AR.index][it_T.index] = ep['u_throat'] # [K] Throat velocity
        Pr_throat[it_AR.index][it_T.index] = fp.get_Prandtl(T=ep['T_throat'], p=ep['p_throat']) # [-] Prandtl number in throat
        mu_throat[it_AR.index][it_T.index] = fp.get_viscosity(T=ep['T_throat'], p=ep['p_throat']) # [Pa*s] Dynamic viscosity in the throat
        # Nozzle exit stuff
        T_exit[it_AR.index][it_T.index] = ep['T_exit'] # [K] Nozzle exit temperature
        p_exit[it_AR.index][it_T.index] = ep['p_exit'] # [Pa] Nozzle exit pressure

# Calculate ideal power
h_inlet = fp.get_enthalpy(T=T_inlet, p=p_chamber) # [J/kg] Enthalpy at chamber inlet
delta_h = h_chamber - h_inlet # [J/kg] Enthalpy change between inlet and chamber
Pt_ideal = m_dot * delta_h # [W] Total ideal power consumption (no heat losses)

## Plot the results

fig, axs = plt.subplots(5,2)


# Render the mass flow and throat width as a function of temperature and area ratios
w_throat = A_throat / h_channel # [m] Throat width, determined by assuming channel depth
# Find out the hydraulic diameter for a rectangular channel
wetted_perimeter  = wetter_perimeter_rectangular(w_channel=w_throat, h_channel=h_channel) # [m] Wetted perimeter
D_hydraulic_throat = hydraulic_diameter(A=A_throat, wetted_perimeter=wetted_perimeter) # [m] Hydraulic diameter
# Numpy and CoolProp don't play nicely together, so iterate over each element to calculate Reynolds's number
it_AR.reset()
for AR in it_AR:
    it_T = np.nditer(T_chamber, flags=['c_index'])
    for T in it_T:
        # Calculate Reynolds number at throat
        Re_throat[it_AR.index][it_T.index] = fp.get_Reynolds_from_velocity(T=T_throat[it_AR.index][it_T.index], p=p_throat[it_AR.index][it_T.index], \
             L_ref=D_hydraulic_throat[it_AR.index][it_T.index], u=u_throat[it_AR.index][it_T.index]) # [-] Reynolds number at throat

it_AR.reset()
for AR in it_AR:
    # Left side of plot
    axs[0][0].plot(T_chamber, m_dot[it_AR.index,:]*1e6, label="{:2.0f}".format(AR))
    axs[1][0].plot(T_chamber, w_throat[it_AR.index,:]*1e6)
    axs[2][0].plot(T_chamber, Isp[it_AR.index,:])
    axs[3][0].plot(T_chamber, Pt_ideal[it_AR.index,:])
    axs[4][0].plot(T_throat[it_AR.index,:], mu_throat[it_AR.index,:])

    # Right side of plot
    axs[0][1].plot(T_chamber, D_hydraulic_throat[it_AR.index,:]*1e6)
    axs[1][1].plot(T_chamber, Re_throat[it_AR.index,:])
    axs[2][1].plot(T_chamber, Pr_throat[it_AR.index,:])
    axs[3][1].plot(T_throat[it_AR.index,:], p_throat[it_AR.index,:]*1e-5)
    axs[4][1].plot(T_exit[it_AR.index,:], p_exit[it_AR.index,:]*1e-5)


# Plot saturation curve on the throat pressure plot to check for condensation
# Some calculations in here are to put proper bound on the plot
T_sat_max = fp.get_saturation_temperature(p=1.1*np.max(p_throat)) #fp.get_critical_temperature() # [K] Get saturation temperature at 1.1 times p_throat
T_sat = np.linspace(start=np.min(T_throat), stop=T_sat_max, num=50) # [K] Evenly spaced temperature from inlet temp to critical temp
p_sat = fp.get_saturation_pressure(T=T_sat) # [Pa] Saturation pressures matching the temps
# Plot it on the throat pressure curve
axs[3][1].plot(T_sat, p_sat*1e-5, label="Saturation curve")

# Again, but with bounds for exit pressures and temperature, to keep plot nice
T_sat_max_2 = fp.get_saturation_temperature(p=1.1*np.max(p_exit)) #fp.get_critical_temperature() # [K] Get saturation temperature at 1 bar
T_sat_2 = np.linspace(start=np.min(T_exit), stop=T_sat_max_2, num=50) # [K] Evenly spaced temperature from inlet temp to critical temp
p_sat_2 = fp.get_saturation_pressure(T=T_sat_2) # [Pa] Saturation pressures matching the temps
axs[4][1].plot(T_sat_2, p_sat_2*1e-5, label="Saturation curve")

# Left side labels
axs[3][0].set_xlabel("Chamber temperature $T_c$ [K]")
axs[0][0].set_ylabel("$\\dot{m}$ [mg/s]")
axs[1][0].set_ylabel("$w_t$ [$\\mu$m]")
axs[2][0].set_ylabel(" $I_{sp}$ [s]")
axs[3][0].set_ylabel("$P_{t,ideal}$ [W]")
axs[4][0].set_ylabel("$\\mu_{throat} [Pa*s]$")
axs[0][0].grid()
axs[1][0].grid()
axs[2][0].grid()
axs[3][0].grid()
axs[4][0].grid()

# Right side labels
axs[0][1].set_ylabel("$D_{h,throat}$ [$\\mu$m]")
axs[1][1].set_ylabel("$Re_{D_{h,throat}}$ [-]")
axs[2][1].set_ylabel("$Pr_{throat}$ [-]")
axs[3][1].set_ylabel("$p_{throat} $ [bar]")
#axs[4][1].set_ylabel("$T_{throat}$ [K]")
axs[1][1].grid()
axs[2][1].grid()
axs[3][1].grid()
axs[4][1].grid()
axs[0][1].grid()

## Print results
m_dot_min = np.min(m_dot)
m_dot_max = np.max(m_dot)
print("Minimum mass flow: {:2.2f} mg/s".format(m_dot_min*1e6))
print("Maximum mass flow: {:2.2f} mg/s".format(m_dot_max*1e6))

# Print results of first point for debugging purposes
print("\nResults of first element in each array")
print("------    AR = {:3.1f}    T_c = {:3.2f} K".format(AR_exit[0],T_chamber[0]))
print("Throat area: {} m^2".format(A_throat[0][0]))
print("D_h: {} m".format(D_hydraulic_throat[0][0]))
print("Throat width: {} m".format(w_throat[0][0]))
print("Throat height {} ".format(h_channel))
print("Wetted perimeter: {} m".format(wetted_perimeter[0][0]))
print("Mass flow: {} kg/s".format(m_dot[0][0]))
print("Throat velocity: {} m/s".format(u_throat[0][0]))
print("Gamma inlet {}".format(fp.get_specific_heat_ratio(T=T_chamber[0],p=p_chamber)))
print("Gamma throat {} ".format(fp.get_specific_heat_ratio(T=T_throat[0][0],p=p_throat[0][0])))
print("T_throat: {} K".format(T_throat[0][0]))
print("p_throat: {} Pa".format(p_throat[0][0]))
print("T_exit: {} K".format(T_exit[0][0]))
print("p_exit: {} Pa".format(p_exit[0][0]))

# Alternative calculatoin of Reynolds number
Re_throat_alt = m_dot[0][0]/A_throat[0][0]*D_hydraulic_throat[0][0]/mu_throat[0][0] # [-] Alternative calc. of Reynodls number m_dot/A instead of rho*u
print("ALT RE: {}".format(Re_throat_alt))
print("Re_throat {} ".format(Re_throat[0][0]))




fig.legend(title="Exit area ratio $\\frac{A_e}{A_t}$", loc='lower right')
fig.suptitle("Peformance vs. $T_c$ & $\\frac{{A_e}}{{A_t}}$ for $F={:2.2f}$ mN and $p_c={:2.0f}$ bar".format(F*1e3,p_chamber*1e-5))
#fig.tight_layout()
plt.show()

