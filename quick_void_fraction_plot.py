## File to show the effect of the density and homogeneous assumption on void fraction and density

import numpy as np
import matplotlib.pyplot as plt

import thermo.two_phase as tp
from thermo.prop import FluidProperties

# Select propellant
fp = FluidProperties("water")

x = np.linspace(start=0, stop=1,num=1000)
p = 1e5*np.array((1,5,10)) # [Pa] Saturation pressure
rho_g = fp.get_vapour_density_at_psat(p_sat=p)
rho_l = fp.get_liquid_density_at_psat(p_sat=p)


p_iter = np.nditer(p, flags=['c_index'])
alpha = np.zeros((np.size(p),np.size(x)))
print(alpha.shape)

# Calculate void fraction for different pressures p
for p in p_iter:
    alpha[p_iter.index] = tp.homogenous_void_fraction(x=x, rho_g=rho_g[p_iter.index], rho_l=rho_l[p_iter.index]) # [-] Void fraction


# Plot it!
plt.figure()
p_iter.reset()

for p in p_iter:
    plt.plot(x, alpha[p_iter.index], label="$p=${:3.0f} bar $\\frac{{\\rho_l}}{{\\rho_g}}={:4.0f}$".format(p*1e-5, rho_l[p_iter.index]/rho_g[p_iter.index]))

plt.legend()
plt.xlabel('Vapour quality $x$ [-]')
plt.ylabel("Void fraction $\\alpha$ [-]")
plt.title("Void fraction  for water (homogeneous)")
plt.tight_layout(pad=0.5)
plt.grid()
plt.show()