# File to validate calculate Cen's thruster for one-D calculations
import numpy as np
import matplotlib.pyplot as plt

import models.one_D as oneD
import thrusters.thruster_data
import thermo.convection
import thermo.two_phase as tp
import basic.chamber
from thermo.prop import FluidProperties

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



# Fidelity of simulatoin
steps_per_section = 50 # [-] Amount of subdivision in each section 
steps_l = steps_per_section
steps_tp = steps_per_section
steps_g = steps_per_section

def calc_and_plot_thruster(td, axs_to_plot):
    # For Cen the chamber temperature is unknown, so  a range is taken instead
    # seem inconsitent with saturation temperatures and/or reported wall temperatures
    T_wall = td['T_wall']                   # [K] Wall temperature
    w_channel = td['w_channel']             # [m] Channel width
    T_inlet = td['T_inlet']                 # [K] Inlet temperature
    p_inlet = td['p_inlet']                 # [Pa] Inlet pressure
    m_dot = td['m_dot']                     # [kg/s] Mass flow (through all channels if multiple)
    channel_amount = td['channel_amount']   # [-] Amount of channels
    h_channel = td['h_channel']             # [m] Channel height/depth
    fp = FluidProperties(td['propellant'])  # Object from which fluid properties can be accessed

    # Calculate mass flow for one single channel
    m_dot_channel = m_dot/channel_amount    # [kg/s] Mass flow through one single channel

    # Chamber temperature is unknown, so a range is taken
    T_sat = fp.get_saturation_temperature(p=p_inlet) # [K]
    T_chamber = np.linspace(start=T_sat+1, stop=T_wall-1, num=250) # [K] 

    # Geometric values
    wetted_perimeter = basic.chamber.wetter_perimeter_rectangular(w_channel=w_channel, h_channel=h_channel) # [m] Wetted perimeter of channel
    A_channel = w_channel*h_channel # [m^2] Cross-sectional through which fluid flows
    D_hydraulic = basic.chamber.hydraulic_diameter_rectangular(w_channel=w_channel, h_channel=h_channel) # [m] Hydraulic diameter
    # Preparation functions calculations many intermediate values that are known before geometry is known


    
    # Storing the length results in here, one for each set of Nusselt relations
    L_1 = np.zeros_like(T_chamber) # [m] Total channel length
    L_2 = np.zeros_like(L_1)

    # Loop to calculate the channel length with each wall temperature
    it_T = np.nditer(T_chamber, flags=['c_index']) # [K] Wall temperature
    T_chamber_guess = None # [K] Stores chamber temperature when actual length is first reached (in laminar case).
    for T in it_T:


        prepared_values = oneD.full_homogenous_preparation(
        T_inlet=T_inlet,
        T_outlet=T, # <---- Iterated variable
        m_dot=m_dot_channel,
        p_ref=p_inlet,
        steps_l=steps_l,
        steps_tp=steps_tp,
        steps_g=steps_g,
        fp=fp)

        results_1 = oneD.full_homogenous_calculation(
            prepared_values=prepared_values,
            Nusselt_relations=Nusselt_relations_1,
            A_channel=A_channel,
            wetted_perimeter=wetted_perimeter,
            D_hydraulic=D_hydraulic,
            m_dot=m_dot,
            T_wall=T_wall,
            p_ref=p_inlet,
            fp=fp
            )
        
        results_2 = oneD.full_homogenous_calculation(
            prepared_values=prepared_values,
            Nusselt_relations=Nusselt_relations_2,
            A_channel=A_channel,
            wetted_perimeter=wetted_perimeter,
            D_hydraulic=D_hydraulic,
            m_dot=m_dot_channel,
            T_wall=T_wall, # <--- Iterated variable
            p_ref=p_inlet,
            fp=fp
            )
        
        # First time the length crosses the actual length, store the value
        if (T_chamber_guess == None) and (results_2['L_total'] > td['L_channel']): # results_2 stores laminar results
            T_chamber_guess = T # [K]
        L_1[it_T.index] = results_1['L_total']
        L_2[it_T.index] = results_2['L_total']

    ## Print info about thruster, including, estimated chamber temperature for laminar relations
    print("Thruster name: {}".format(td['name']))
    print("Estimated chamber temperature: {:3.2f} K".format(T_chamber_guess))
    

    axs_to_plot.set_title("$\\dot{{m}}={:1.2f}$ mg/s, $p={:1.2f}$ bar".format(m_dot*1e6,p_inlet*1e-5))
    axs_to_plot.plot(T_chamber,L_1*1e3, label="Turbulent")
    axs_to_plot.plot(T_chamber,L_2*1e3, label="Laminar")
    axs_to_plot.hlines(td['L_channel']*1e3, xmin=T_chamber[0], xmax=T_chamber[-1], linestyle='dashed', color='red', label="Real length")
    axs_to_plot.grid()

        

def run():
    # Put results of each Cen2010 thruster in subplot
    fig, axs = plt.subplots(3,2)

    td_1 = thrusters.thruster_data.Cen2010_1 # Dictionary with design/measured values
    calc_and_plot_thruster(td_1, axs[0][0])

    td_2 = thrusters.thruster_data.Cen2010_2 # Dictionary with design/measured values
    calc_and_plot_thruster(td_2, axs[1][0])

    td_3 = thrusters.thruster_data.Cen2010_3 # Dictionary with design/measured values
    calc_and_plot_thruster(td_3, axs[2][0])

    td_4 = thrusters.thruster_data.Cen2010_4 # Dictionary with design/measured values
    calc_and_plot_thruster(td_4, axs[0][1])

    td_5 = thrusters.thruster_data.Cen2010_5 # Dictionary with design/measured values
    calc_and_plot_thruster(td_5, axs[1][1])

    td_6 = thrusters.thruster_data.Cen2010_6 # Dictionary with design/measured values
    calc_and_plot_thruster(td_6, axs[2][1])

    # Change y-limits on axes
    axs[0][0].set_ylim([0,8])
    axs[1][0].set_ylim([0,8])
    axs[2][0].set_ylim([0,8])
    axs[0][1].set_ylim([0,8])
    axs[1][1].set_ylim([0,8])
    axs[2][1].set_ylim([0,8])

    # Put x-labels on right places
    axs[2][0].set_xlabel("Chamber temperature - $T_c$ [K] ")
    axs[2][1].set_xlabel("Chamber temperature - $T_c$ [K] ")
    # Put y-labels on right places
    axs[0][0].set_ylabel("$L$ [mm]")
    axs[1][0].set_ylabel("Channel length - $L$ [mm]")
    axs[2][0].set_ylabel("$L$ [mm]")

    axs[1][1].legend(bbox_to_anchor=(1.0, 0.4))
    plt.suptitle("1D model applied to Cen2010 thruster data")
    plt.tight_layout(pad=0.5)
    plt.show() 

def subplot_results():
    pass


if __name__ == "__main__":
    print("Running validation script for Cen2010")
    run()