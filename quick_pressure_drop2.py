import thrusters.thruster_data
from basic.chamber import hydraulic_diameter_rectangular, velocity_from_mass_flow
from thermo.prop import FluidProperties
# Load a specific thruster for quick pressure drop calculations
td = thrusters.thruster_data.td_Silva_5

fp = FluidProperties(td['propellant']) # Object to access fluid properties from
m_dot = td['m_dot']/td['channel_amount'] # [kg/s] Mass flow
h_channel = td['h_channel'] # [m] Channel height/depth
w_channel = td['w_channel'] # [m] Channel width
L_channel = td['L_channel'] # [m] Channel length
T_outlet = td['T_chamber']  # [K] Outlet of channel
p_inlet = td['p_inlet'] # [Pa]
p_outlet = p_inlet # [Pa] Outlet pressure (assumed roughly equal to inlet)
T_inlet = td['T_inlet'] # [K]

# Geometry calculations
D_h = hydraulic_diameter_rectangular(w_channel=w_channel, h_channel=h_channel) # [m] Hydraulic diameter
A_channel = w_channel*h_channel
# Inlet conditions
rho_inlet = fp.get_density(T=T_inlet, p=p_inlet) # [kg/m^3]
# Outlet conditions
rho_outlet = fp.get_density(T=T_outlet, p=p_outlet)
u_outlet = velocity_from_mass_flow(A=A_channel, m_dot=m_dot, rho=rho_outlet)
Re_outlet = fp.get_Reynolds_from_mass_flow(T=T_outlet, p =p_outlet, L_ref=D_h, m_dot=m_dot, A=A_channel)


print(" --- GEOMETRY ---")
print("Channel length: {:4.3f} mm".format(L_channel*1e3))
print("Hydraulic diameter: {:4.3f} um".format(D_h*1e6))
print("L/D: {:3.2f} ".format(L_channel/D_h))

print("\n --- INLET FLOW CONDITIONS --- ")
print("Density: {:3.2f} kg/m^3".format(rho_inlet))

print("\n --- OUTLET FLOW CONDITIONS ---")
print("Mass flow: {:1.3f} mg/s".format(m_dot*1e6))
print("Density: {:1.3f} kg/m^3".format(rho_outlet))
print("Velocity: {:3.3f} m/s".format(u_outlet))
print("Reynolds: {:3.4f} ".format(Re_outlet))


# Quick and dirty pressure drop estimates
f = 64/Re_outlet # [-] Friction factor
print("\n --- FRICTIONAL PRESSURE DROP --- ")
delta_P = f * L_channel/D_h * 0.5 * rho_outlet * u_outlet**2
print("Friction factor: {:1.4f}".format(f))
print("Pressure drop: {:1.5f} bar".format(delta_P*1e-5))
print("Delta P / P: {} ".format(delta_P/p_outlet))

mass_flux = m_dot/A_channel
delta_P_accelaration = mass_flux**2 * ( 1/rho_inlet - 1/rho_outlet) # [Pa] Pressure drop due to acceleration
print("\nAccelerational pressure drop: ")
print("Mass flux: {} kg/m^2".format(mass_flux))
print("Delta P: {:1.6f} bar".format(delta_P_accelaration*1e-5))