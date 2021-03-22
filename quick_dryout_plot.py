## File to show results of crititical dryout-formula

import numpy as np
import matplotlib.pyplot as plt
import thermo.two_phase as tp


p = 5e5 # [Pa] Pressure in heater
D_hydr = 0.1e-3 # [m] Hydraulic diameter
G = np.linspace(start=20, stop=1000, num=1000) # [kg/m^2] Mass fluxes
m_dot = 1 # [kg/s] Mass flow, just a placeholder for this plot, as only mass flux matters
A_channel = m_dot / G # [m^2] Channel area

x_crit = tp.dryout_quality(p=p, m_dot=m_dot, A_channel=A_channel, D_hydr=D_hydr) # [-] Dry out quality as function of mass flux

## Plot to see if it replicates figure on page 633 of book Carey2008
plt.figure()
plt.plot(G,x_crit, label="{:1.1f} mm".format(D_hydr*1e3))

D_hydr = 0.2e-3 # [m] Hydraulic diameter
x_crit = tp.dryout_quality(p=p, m_dot=m_dot, A_channel=A_channel, D_hydr=D_hydr) # [-] Dry out quality as function of mass flux
plt.plot(G,x_crit, label="{:1.1f} mm".format(D_hydr*1e3))

D_hydr = 0.5e-3 # [m] Hydraulic diameter
x_crit = tp.dryout_quality(p=p, m_dot=m_dot, A_channel=A_channel, D_hydr=D_hydr) # [-] Dry out quality as function of mass flux
plt.plot(G,x_crit, label="{:1.1f} mm".format(D_hydr*1e3))

D_hydr = 1e-3 # [m] Hydraulic diameter
x_crit = tp.dryout_quality(p=p, m_dot=m_dot, A_channel=A_channel, D_hydr=D_hydr) # [-] Dry out quality as function of mass flux
plt.plot(G,x_crit, label="{:1.1f} mm".format(D_hydr*1e3))

plt.hlines(1,xmin=np.min(G),xmax=np.max(G), colors="gray", linestyles="--", label="$x=1$")

plt.xlabel("Mass flux  $G$ [kg$\\cdot$m$^{-2}$]")
plt.ylabel("Critical quality $x_{{crit}}$ [-]")
plt.title("Dry-out quality at $p={:1.0f}$ bar".format(p*1e-5))
plt.legend(title="Hydraulic diameter  $D_h$:")
plt.tight_layout(pad=0.5)
plt.grid()
plt.show()
